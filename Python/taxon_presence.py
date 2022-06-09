#!/usr/bin/env python3

# This script gets the presence of all classified taxa between the tools and outputs a binary matrix

import json, os, argparse
import pandas as pd


def main():
    files, outdir = get_intake()
    programs_taxa = get_taxa(files)
    df_dict = get_df(programs_taxa)
    save_csv(df_dict, outdir)

def get_intake():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--infiles",
                        nargs = "+",
                        help = "Input json files from 'meta_tools_taxon_benchmark.py' It can be called multiple times.")
    parser.add_argument("-out", "--outdir",
                        help = "Sets the directory to save the output files",
                        default = os.path.abspath(os.path.dirname(__file__)))
    
    args = parser.parse_args()
    files = args.infiles
    outdir = os.path.abspath(args.outdir)

    return files, outdir


def get_taxa(files):
    ranks = ["phylum", "class", "order", "family", "genus", "species", "NULL"]

    programs_taxa = {}
    for file in files:
        if any(st in file for st in ["kaiju", "kj"]):
            program = "kaiju"
        elif any(st in file for st in ["kraken", "krk"]):
            program = "kraken2"
        elif any(st in file for st in ["centrifuge", "cfg"]):
            program = "centrifuge"
        elif any(st in file for st in ["metacache", "mc_out"]):
            program = "metacache"
        elif any(st in file for st in ["mapping"]):
            program = "cami"
        else:
            print("Program not recognized or incompatible.")
        with open(file) as jfile:
            print(f"Working with {file}")
            programs_taxa[program] = {}
            jdict = json.load(jfile)
            for tax_level in ranks[:-1]:
                programs_taxa[program][tax_level] = set([])
                for read, tax_info in jdict.items():
                    if tax_info[0] != "0":
                        programs_taxa[program][tax_level].add(tax_info[0].get(tax_level, "NULL"))
                
                programs_taxa[program][tax_level] = set([y for y in programs_taxa[program][tax_level] if y != "NULL"])

    return programs_taxa

def get_df(programs_taxa):
    df_dict = {}
    for tax_level in ["order", "family", "genus", "species"]:
        dict_csv = {}
        for program in programs_taxa:
            for taxon in programs_taxa[program][tax_level]:
                if dict_csv.get(taxon, None) is None:
                    dict_csv[taxon] = {}
                    dict_csv[taxon][program] = "1"
                elif dict_csv[taxon]:
                    dict_csv[taxon][program] = "1"
        df = pd.DataFrame.from_dict(dict_csv).fillna(0)
        df_dict[tax_level] = df
    
    return df_dict

def save_csv(df_dict, outdir):
    for tax_level in ["order", "family", "genus", "species"]:
        df_dict[tax_level].to_csv(os.path.join(outdir, f"taxa_matrix_{tax_level}.csv"))


if __name__ == '__main__':
    main()
