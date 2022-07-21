#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:35:43 2022

@author: chrisbl
"""

import pyranges as pr

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib/Homo_sapiens.GRCh38.107.gtf

def extract_protein_coding_ids(ensembl_path):
    
    print("Loading local ensembl dataframe...")
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    print("Extracting IDs of protein coding transcripts...")
    protein_coding_ids = [ (gene_id, transcript_id) for gene_id,
                                     transcript_id,
                                     biotype in zip(gtf['gene_id'],
                                                    gtf['transcript_id'],
                                                    gtf['transcript_biotype']) if biotype == "protein_coding"]
    print("Assembling IDs...")
    transcript_dict = dict()
    transcript_list = set()
    protein_coding_gene_ids = set()
    for key, transcript_id in protein_coding_ids:
        transcript_dict[key] = []
        protein_coding_gene_ids.add(key)
    for key, transcript_id in protein_coding_ids:
        transcript_dict[key].append(transcript_id)
        transcript_list.add(transcript_id)
    print("Sorting IDs...")
    for key in transcript_dict.keys():
        transcript_dict[key] = set(transcript_dict[key])
        transcript_dict[key] = sorted(list(transcript_dict[key]))
    protein_coding_gene_ids = sorted(list(protein_coding_gene_ids))
    transcript_list = list(transcript_list)
    
    return transcript_dict, transcript_list, protein_coding_gene_ids
    
def main():
    transcript_dict, transcript_list, protein_coding_gene_ids = extract_protein_coding_ids("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/Homo_sapiens.GRCh38.107.gtf")
    print(transcript_list[0:20])
    count = 0
    flag = False
    for key in transcript_dict.keys():
        count += 1
        for ID in transcript_dict[key]:
            if ID == "E":
                print("KEY IS ", key, "::")
                print(transcript_dict[key])
        if flag:
            break
                
                

if __name__ == "__main__":
    main()