#!/usr/bin/env python3

# This script counts how many taxons there are per taxonomic level.

import json, os, csv, argparse

directory = os.path.abspath(os.path.dirname(__file__))

parser = argparse.ArgumentParser()
parser.add_argument("-gs", "--gold_standard",
                    nargs="+",
                    help="Input reads mapping json files from the cami")
parser.add_argument("-out", "--outdir",
                    help="Sets the directory to save the output files",
                    default=os.path.abspath(os.path.dirname(__file__)))
args = parser.parse_args()
gs_files = args.gold_standard
outdir = os.path.abspath(args.outdir)

results = []
ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

for sample_id in ["s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"]:
    gs_file = "".join([file for file in gs_files if sample_id in file])
    if gs_file == "":
        break
    # Loads the CAMI json file
    gs_result = {"superkingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    with open(gs_file) as jfile:
        print(f"Loading {gs_file}")
        gs_dict = json.load(jfile)
        # Gets only the selected taxonomy level
        print(f"Began analysis")
        for read, taxonomy in gs_dict.items():
            taxonomy = taxonomy[0]
            for rank in ranks:
                if taxonomy.get(rank, "none") != "none":
                    gs_result[rank].add(taxonomy[rank])
        gs_sum = {"sample": sample_id}
        for rank in ranks:
            gs_sum[rank] = len(gs_result[rank])
        results.append(gs_sum)

# Saves the output as a CSV file
print("Saving output")
with open(f"{outdir}/count_diversity_cami.csv", "w+") as outfile:
    fieldnames = ["sample", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    out = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=",")
    out.writeheader()
    for result in results:
        out.writerow(result)
