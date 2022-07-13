#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:35:43 2022

@author: chrisbl
"""

import pyranges as pr

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib/Homo_sapiens.GRCh38.107.gtf

def extract_protein_coding_ids(ensembl_path):
    
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    protein_coding_ids = [ (gene_id, transcript_id) for gene_id,
                                     transcript_id,
                                     biotype in zip(gtf['gene_id'],
                                                    gtf['transcript_id'],
                                                    gtf['transcript_biotype']) if biotype == "protein_coding"]
    transcript_dict = dict()
    transcript_list = list()
    protein_coding_gene_ids = set()
    for key, transcript_id in protein_coding_ids:
        transcript_dict[key] = []
        protein_coding_gene_ids.add(key)
    for key, transcript_id in protein_coding_ids:
        transcript_dict[key].append(transcript_id)
        transcript_list.append(transcript_id)
    for key in transcript_dict.keys():
        transcript_dict[key] = set(transcript_dict[key])
        transcript_dict[key] = sorted(list(transcript_dict[key]))
    protein_coding_gene_ids = sorted(list(protein_coding_gene_ids))
    
    return transcript_dict, transcript_list, protein_coding_gene_ids
    
def main():
    extract_protein_coding_ids()

if __name__ == "__main__":
    main()