#!/usr/bin/env python3

from Bio import Entrez
import csv, os, json, argparse
from collections import defaultdict

def main():
    infiles, outdir = get_user_input()
    programs_configs = {"standard": {"reads_col": 0, "taxid_col": 1, "delimiter": "\t", "n_header": 0                                       },
                        "kaiju": {"reads_col": 1, "taxid_col": 2, "delimiter": "\t", "n_header": 0},
                        "kraken2": {"reads_col": 1, "taxid_col": 2, "delimiter": "\t", "n_header": 0},
                        "centrifuge": {"reads_col": 0, "taxid_col": 2, "delimiter": "\t", "n_header": 1},
                        "metacache": {"reads_col": 0, "taxid_col": 1, "delimiter": "\t", "n_header": 11}}
    for infile in infiles:
        if any(st in infile for st in ["mapping", "reads"]):
            config = programs_configs["standard"]
        elif any(st in infile for st in ["kaiju", "kj"]):
            config = programs_configs["kaiju"]
        elif any(st in infile for st in ["kraken", "krk"]):
            config = programs_configs["kraken2"]
        elif any(st in infile for st in ["centrifuge", "cfg"]):
            config = programs_configs["centrifuge"]
        elif any(st in infile for st in ["metacache", "mc_out"]):
            config = programs_configs["metacache"]
        else:
            print("Program not recognized or incompatible.")
        infile_name = infile.split("/")[-1]
        reads_dict, taxIDs = get_infile_data(infile, config)
        sorted_tax_data = get_tax_data(taxIDs)
        sorted_reads_dict = sort_save_read_data(sorted_tax_data, reads_dict,outdir,infile_name)

def get_user_input():
    """Deals with the user commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--infiles",
                        nargs="+",
                        help = "Input files with reads and taxonomy IDs. It can be called multiple times.\nCurrently supported: Standard, Clark, Kraken 2, Centrifuge, Metacache",)
    parser.add_argument("-out", "--outdir",
                        help = "Sets the directory to save the output files",
                        default=os.path.abspath(os.path.dirname(__file__)))
    args = parser.parse_args()
    return args.infiles, args.outdir

def get_infile_data(infile, program_configs):  # Input is a CSV file from the programs
    """Deals with the program reads output"""
    print(f"\nGetting data from {infile}")
    reads_dict = defaultdict(list)
    taxIDs = set()
    # Reads the output file
    with open(infile, "r") as csv_in:
        in_reader = csv.reader(csv_in, delimiter=program_configs["delimiter"])
        # Only skips i there is a header
        for i in range(0, program_configs["n_header"]):
            next(in_reader)
        # Gets the readID and taxID
        for row in in_reader:
            if not row[0].startswith("#"):
                reads_dict[row[program_configs["reads_col"]]].append(row[program_configs["taxid_col"]]) # reads_dict = {readID_1: [taxID_1, taxID_1], [...], readID_n: [taxID_n]}
                if row[program_configs["taxid_col"]] not in ["0", "NA"]:
                    taxIDs.add(row[program_configs["taxid_col"]])  # taxIDs = [taxID_1, taxID_2, [...], taxID_n]
    taxIDs = list(taxIDs)
    print(f"Got data from {len(reads_dict)} reads and {len(taxIDs)} TaxIDs")
    return reads_dict, taxIDs

def get_tax_data(taxids_list): # Input is a list with TaxIDs
    """Deals with getting the data for each TaxID"""
    print("Getting TaxIDs data")
    Entrez.email = 'gabrielamorimsilva@gmail.com'
    Entrez.api_key = 'f19b2580b4e240476bdef13bab28f8bb7808'
    taxids_data = []
    # Creates sub lists of 10000 TaxIDs
    sub_taxid_list = [taxids_list[i:i + 10000] for i in range(0, len(taxids_list), 10000)]
    # Gets the data from NCBI
    for sub_list in sub_taxid_list:
        sub_taxid_str = ','.join(sub_list)
        search = Entrez.efetch(id=sub_taxid_str, db="taxonomy", retmode="xml")
        taxids_data.extend(Entrez.read(search)) # [{'TaxId': '1200747', [...], 'Rank': 'species', [...],'LineageEx': [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'}, [...]][...]]
    sorted_tax_data = sort_tax_data(taxids_data)
    print(f"Got data from {len(taxids_data)} TaxIDs")
    return sorted_tax_data

def sort_tax_data(taxids_data): # Input is the list with the data from Entrez
    """Sorts the data in the ranks for later comparison"""
    print("Sorting the data")
    sorted_tax_data = {} # {taxid_1: {species: taxid_s, genus: taxid_g, [...]}, taxid_2: {species: taxid_s, genus: taxid_g, [...]} [...]}
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    # Sorts the LineageEx taxIDs from the data
    for entry in taxids_data:
        if entry["Rank"] == "superkingdom":
            sorted_tax_data[entry["TaxId"]] = {"superkingdom": entry["TaxId"]}
        elif entry["TaxId"] in ["131567", "1"]:
            sorted_tax_data[entry["TaxId"]] = {}
        else:
            sorted_tax_data[entry["TaxId"]] = {}
            for item in entry["LineageEx"]:
                if item["Rank"] in ranks:
                    sorted_tax_data[entry["TaxId"]][item["Rank"]] = item["TaxId"]
            if entry["Rank"] in ranks:
                sorted_tax_data[entry["TaxId"]][entry["Rank"]] = entry["TaxId"]
            if entry.get("AkaTaxIds", False):
                for alt_taxid in entry["AkaTaxIds"]:
                    sorted_tax_data[alt_taxid] = sorted_tax_data[entry["TaxId"]].copy()
                    if entry["Rank"] in ranks:
                        sorted_tax_data[alt_taxid][entry["Rank"]] = alt_taxid
    return sorted_tax_data


def sort_save_read_data(sorted_tax_data, reads_dict, outdir, infile_name):
    """Substitutes the taxID by the taxonomic information """
    reads_tax_dict = defaultdict(list)
    for read, r_taxids in reads_dict.items():
        for r_taxid in r_taxids:
            if r_taxid not in ["NA", "0"]:
                reads_tax_dict[read].append(sorted_tax_data[r_taxid])
            elif r_taxid in ["NA", "0"]:
                reads_tax_dict[read].append(r_taxid)
    print("Saving output dictionary json file")
    with open(f"{outdir}/{infile_name}-dict.json", "w+") as json_file:
        json.dump(reads_tax_dict, json_file)
    return reads_dict

if __name__ == "__main__":
    main()
