#! /usr/bin/env python3

from Bio import Entrez
import csv, argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cami_files", nargs="+", help="Input CAMI reads mapping files")

    cami_files = parser.parse_args().cami_files

    for file in cami_files:
        print(f"Working with {file}")
        
        accession_reads = {}  # {accession:[read, read], accession:[read, read]}
        
        with open(file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                read = row[0]
                accession = row[3].split(".")[0]
                if accession not in accession_reads:
                    accession_reads[accession] = [read]
                elif accession in accession_reads:
                    accession_reads[accession].append(read)
        
        print("Getting taxids")
        taxid_accessions = get_Entrez_data(list(accession_reads.keys())) # {taxid:[accession, accession], taxid:[accession, accession]}
        
        print("Saving output")
        with open(f"{file.split('.')[0]}_corrected.tsv", "w+") as f_out:
            writer = csv.writer(f_out, delimiter="\t")
            for taxid, accessions in taxid_accessions.items():
                for accession in accessions:
                    for read in accession_reads[accession]:
                        writer.writerow([read, taxid, accession])

def get_Entrez_data(accessions): 
    Entrez.email = 'EMAIL'
    Entrez.api_key = 'API_KEY'
    accessions_taxid = {}
    # Creates sub lists of 10000 accessions
    sub_accessions = [accessions[i:i + 10000] for i in range(0, len(accessions), 10000)]
    for sub_list in sub_accessions:
        search = Entrez.efetch(id=sub_list, db="nucleotide", retmode="xml")
        records = Entrez.read(search)
        for record in records:
            accession = record["GBSeq_locus"]
            for i in record["GBSeq_feature-table"][0]["GBFeature_quals"]:
                if i['GBQualifier_name'] == "db_xref":
                    taxid = i['GBQualifier_value'].split(":")[-1]
            if taxid not in accessions_taxid:
                accessions_taxid[taxid] = [accession]
            elif taxid in accessions_taxid:
                accessions_taxid[taxid].append(accession)
    return accessions_taxid

if __name__ == '__main__':
    main()
