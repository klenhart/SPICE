#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:35:43 2022

@author: chrisbl
"""

import pyranges as pr


tsl_dict = {"1":"1", "2":"2", "3":"3", "4":"4", "5":"5", "n":"6", "N":"6"}

def load_gtf(ensembl_path):
    print("Loading local ensembl dataframe...")
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    reduced = [ (gene_id, protein_id, transcript_id, biotype, tsl, tag) for gene_id,
                                      protein_id,
                                      transcript_id,
                                      biotype,
                                      tsl,
                                      tag in zip(gtf['gene_id'],
                                                 gtf['protein_id'],
                                                 gtf['transcript_id'],
                                                 gtf['transcript_biotype'],
                                                 gtf['transcript_support_level'],
                                                 gtf['tag'])]
    return reduced

def id_checker(gene_id, protein_id, transcript_id, biotype):
    flag_biotype = biotype in ("protein_coding", "protein_coding_LoF")
    flag_protein_id = str(protein_id).startswith("ENS") and (protein_id[3] == "P" or protein_id[6] == "P")
    return flag_biotype and flag_protein_id

def extract_protein_coding_ids(ensembl_path):
    gtf = load_gtf(ensembl_path)
    print("Extracting IDs of protein coding transcripts...")
    protein_coding_ids = [ (gene_id, protein_id, transcript_id, tsl_dict[str(tsl)[0]], tag) for gene_id, 
                          protein_id,
                          transcript_id,
                          biotype,
                          tsl,
                          tag in gtf if id_checker(gene_id, protein_id, transcript_id, biotype) ]
    protein_coding_ids_no_dups = list(set([i for i in protein_coding_ids]))
    protein_coding_ids_no_dups.sort()
    return protein_coding_ids_no_dups
    
def main():
    pass


if __name__ == "__main__":
    main()