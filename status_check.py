import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
import urllib3
from pathlib import Path

# Disable HTTPS warning
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------

def fetch_plasmid(id_, session):
    """Fetch plasmid metadata → returns a list of rows (long format)."""
    url = (
        "https://bioregistration.onetakeda.com:9943"
        f"/bioreganon/api/v1/data/get-data/plasmid-name/{id_}"
    )

    try:
        r = session.get(url, verify=False, timeout=10)
        r.raise_for_status()
    except Exception:
        return []

    data = r.json()
    if not data:
        return []

    entry = data[0]
    lotId = entry.get("lotId")
    raw_biomass = entry.get("biomassName", None)

    # Normalize biomassName field
    if raw_biomass is None:
        biomass_list = []
    elif isinstance(raw_biomass, str):
        biomass_list = [raw_biomass]
    else:
        biomass_list = list(raw_biomass)

    return [{"plasmid_id": id_, "lotId": lotId, "biomassName": b} for b in biomass_list]


def parallel_fetch_plasmids(id_list, session, max_workers=20):
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(fetch_plasmid, pid, session): pid for pid in id_list}
        for fut in as_completed(futures):
            try:
                results.extend(fut.result())
            except Exception:
                pass
    return results


def fetch_biomass(biomassName, session):
    url = (
        "https://bioregistration.onetakeda.com:9943"
        f"/bioreganon/api/v1/data/get-data/biomass/{biomassName}"
    )

    try:
        r = session.get(url, verify=False, timeout=10)
        r.raise_for_status()
    except Exception:
        return None

    data = r.json()

    # Normalize API returning dict vs list
    if isinstance(data, dict):
        data = [data]
    if not data:
        return None

    entry = data[0]

    # Normalize attachments (may be None/str/list)
    attachments = entry.get("attachments")
    if attachments is None:
        attachments = []
    elif isinstance(attachments, str):
        attachments = [attachments]
    elif isinstance(attachments, list):
        attachments = attachments
    else:
        attachments = []

    return {
        "biomassName": biomassName,
        "biomass_lot_id": entry.get("biomass_lot_id"),
        # NOTE: your code stores attachments in biomass_description
        "biomass_description": attachments,
    }


def parallel_fetch_biomass(biomass_list, session, max_workers=20):
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(fetch_biomass, b, session): b for b in biomass_list}
        for fut in as_completed(futures):
            try:
                info = fut.result()
                if info:
                    results.append(info)
            except Exception:
                pass
    return results


def get_matching_biomasses(row, df_plasmid_full):
    ids = row["IDs"]
    date_str = str(row["date"])

    subset = df_plasmid_full[df_plasmid_full["plasmid_id"].isin(ids)]

    matches = subset[
        subset["biomass_description"]
        .astype(str)
        .str.contains(date_str, na=False)
    ]

    return pd.Series({
        "match_biomassName": matches["biomassName"].tolist(),
        "match_biomass_lotID": matches["biomass_lot_id"].tolist(),
        "match_biomass_description": matches["biomass_description"].tolist(),
    })


def dedupe_any_list(x):
    if x is None:
        return []
    if not isinstance(x, list):
        return [x]

    flat = []
    for item in x:
        if isinstance(item, list):
            flat.extend(item)
        else:
            flat.append(item)

    seen = set()
    deduped = []
    for item in flat:
        if item not in seen:
            seen.add(item)
            deduped.append(item)
    return deduped


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Load file names from a text file into a DataFrame"
    )
    parser.add_argument(
        "-input",
        required=True,
        help="Path to input text file (e.g., FileNames_Biortus.txt)"
    )
    args = parser.parse_args()

    input_file = args.input


    p = Path(input_file)



    with open(input_file, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    df = pd.DataFrame(lines[1:], columns=["File Name"])

    # Extract 8-digit date
    df["date"] = df["File Name"].str.extract(r'(\d{8})')

    # Extract all matching IDs: SECC-, SBVC-, or SMCC-
    df["IDs"] = df["File Name"].str.findall(r'(?:SECC|SBVC|SMCC)-\d+')

    # All unique plasmid IDs
    all_ids = sorted({pid for row in df["IDs"] for pid in row})

    # Create one session for this run
    with requests.Session() as session:
        plasmid_rows = parallel_fetch_plasmids(all_ids, session)
        df_plasmid = pd.DataFrame(plasmid_rows)

        all_biomass = df_plasmid["biomassName"].dropna().unique()
        biomass_rows = parallel_fetch_biomass(all_biomass, session)
        df_biomass = pd.DataFrame(biomass_rows)

    df_plasmid_full = df_plasmid.merge(df_biomass, on="biomassName", how="left")

    df_result = pd.concat(
        [df, df.apply(get_matching_biomasses, axis=1, df_plasmid_full=df_plasmid_full)],
        axis=1
    )

    df_result["match_biomassName"] = df_result["match_biomassName"].apply(dedupe_any_list)
    df_result["match_biomass_lotID"] = df_result["match_biomass_lotID"].apply(dedupe_any_list)
    df_result["match_biomass_description"] = df_result["match_biomass_description"].apply(dedupe_any_list)
    df_result["match_count"] = df_result["match_biomassName"].apply(len)

    out_path = p.with_name(f"{p.stem}_status.tsv")  # same folder as input file
    df_result.to_csv(out_path, sep="\t", index=False)
    print("Wrote:", out_path.resolve())


if __name__ == "__main__":
    main()
