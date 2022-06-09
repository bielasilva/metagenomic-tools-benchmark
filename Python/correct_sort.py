#!/usr/bin/env python3

import json, os, csv, argparse

directory = os.path.abspath(os.path.dirname(__file__))

parser = argparse.ArgumentParser()
parser.add_argument("-in", "--infiles",
                    nargs="+",
                    help="Input files with reads and taxonomy IDs in JSON format from results_to_json.py. It can be called multiple times.\nSupported: CAMI, Kaiju, Kraken 2, Centrifuge, Metacache")
parser.add_argument("-gs", "--gold_standard",
                    nargs="+",
                    help="Input reads mapping files")
parser.add_argument("-t", "--taxonomy",
                    help="Taxonomy level. Currently accepts superkingdom, phylum, class, order, species. Default is species",
                    default="species")
parser.add_argument("-out", "--outdir",
                    help="Sets the directory to save the output files",
                    default=os.path.abspath(os.path.dirname(__file__)))
args = parser.parse_args()
tax_level = args.taxonomy
gs_files = args.gold_standard
files = args.infiles
outdir = os.path.abspath(args.outdir)

results = []
ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "NULL"]
low_tax_level = ranks[ranks.index(tax_level) + 1]

# Checks and divide the files per samples
for sample_id in ["s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"]:
    sample_files = [file for file in files if sample_id in file]
    gs_file = "".join([file for file in gs_files if sample_id in file])
    
    print(f"Working with: {sample_files} and {gs_file}")
    
    # Loads the CAMI json file
    gs = {}
    with open(gs_file) as jfile:
        print(f"Loading {gs_file}")
        gs_dict = json.load(jfile)
        # Gets only the selected taxonomy level
        for read, taxonomy in gs_dict.items():
            try:
                gs[read] = taxonomy[0][tax_level]
            except KeyError:
                gs[read] = taxonomy[0][low_tax_level]
        gs_size = len(gs)

    # Compares the results from the tools with the CAMI standard
    for file in sample_files:
        print(f"Working with {file}")
        
        # Deals with the different tools
        if any(st in file for st in ["kaiju", "kj"]):
            program = "kaiju"
        elif any(st in file for st in ["kraken", "krk"]):
            program = "kraken2"
        elif any(st in file for st in ["centrifuge", "cfg"]):
            program = "centrifuge"
        elif any(st in file for st in ["metacache", "mc_out"]):
            program = "metacache"
        else:
            print("Program not recognized or incompatible.")
        
        unclassified = 0
        correct = 0
        incorrect = 0

        with open(file, "r") as jfile:
            print(f"Loading {file}")
            reads_dict = json.load(jfile)
            for read,taxonomies in reads_dict.items():
                # Checks for the selected taxonomy level and assigns the read to each category
                for taxonomy in taxonomies:
                    if taxonomy in ["NA", "0"]:
                        unclassified += 1
                    elif taxonomy.get(tax_level, "none") == "none":
                        if taxonomy.get(low_tax_level) == gs[read]:
                            correct += 1
                        else:
                            incorrect += 1
                    elif taxonomy.get(tax_level) == gs[read]:
                        correct += 1
                    else:
                        incorrect += 1
            total = correct + incorrect + unclassified
            # Checks if total size is compatible with the CAMI standard
            if total != gs_size:
                print(f"Total sum of {file} is incorrect. Real is {gs_size} but got {total}")
            
        # Append to the results list to be saved later
        results.append([program, sample_id, correct, incorrect, unclassified])

# Saves the output as a CSV file
print("Saving output")
with open(f"{outdir}/statistics_result.csv", "w+") as outfile:
    out = csv.writer(outfile, delimiter=",")
    out.writerow(["program", "sample", "correct", "incorrect", "unclassified"])
    out.writerows(results)