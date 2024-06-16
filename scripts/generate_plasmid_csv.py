#!/usr/bin/env python3
import json
import argparse

# json_path = "/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_plasmid_data/GCF_929606205.1/data.json"

def extract_hits(data):
    hits = []
    for key, value in data.items():
        if isinstance(value, dict):
            if "plasmid" in value:  # Check if this dict represents a hit
                hit = value.copy()
                hit["hit_id"] = key
                hits.append(hit)
            else:
                hits.extend(extract_hits(value))
    return hits

def json_to_csv(json_file):
    with open(json_file, 'r') as file:
        parsed = json.load(file)
        # print(json.dumps(parsed, indent=4))

        hits = extract_hits(parsed["plasmidfinder"]["results"])

        if hits:
            # Specify the column order
            keys = ["plasmid", "identity", "HSP_length", "template_length", "position_in_ref",
                    "contig_name", "positions_in_contig", "accession", "coverage", "hit_id"]

            # Print column headers
            print(",".join(keys))

            # Print each hit as a CSV row
            for hit in hits:
                values = []
                for key in keys:
                    value = hit.get(key, "")
                    if key in ["contig_name", "hit_id"] and "," in value:
                        value = value.replace(",", "")  # Remove commas from contig_name and hit_id
                    values.append(str(value))
                print(",".join(values))


def main(json_file):
    json_to_csv(json_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PlasmidFinder JSON FILE")
    parser.add_argument("--json_path", type=str, help="Path to PlasmidFinder JSON output.")

    args = parser.parse_args()

    main(args.json_path)

