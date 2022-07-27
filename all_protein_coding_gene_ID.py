#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:35:43 2022

@author: chrisbl
"""

import pyranges as pr

prefix = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/"
ensembl_path1 = "Homo_sapiens.GRCh38.107.gtf"
ensembl_path2 = "Tetraodon_nigroviridis.TETRAODON8.107.gtf"
    


# /share/project/zarnack/chrisbl/FAS/utility/protein_lib/Homo_sapiens.GRCh38.107.gtf

def extract_protein_coding_ids(ensembl_path):
    
    print("Loading local ensembl dataframe...")
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    print("Extracting IDs of protein coding transcripts...")
    protein_coding_ids = [ (gene_id, protein_id) for gene_id,
                                      protein_id,
                                      biotype in zip(gtf['gene_id'],
                                                    gtf['protein_id'],
                                                    gtf['transcript_biotype']
                                                    ) if biotype == "protein_coding" and str(type(protein_id)) != "<class 'float'>"]

    protein_coding_ids = list(set([i for i in protein_coding_ids]))
        
    print("Assembling IDs...")
    protein_coding_ids.sort()
    return protein_coding_ids
    
def main():
    prefix = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/"
    ensembl_path1 = "Homo_sapiens.GRCh38.107.gtf"
    ensembl_path2 = "Tetraodon_nigroviridis.TETRAODON8.107.gtf"
    protein_coding_ids = extract_protein_coding_ids(prefix + ensembl_path2)

                

if __name__ == "__main__":
    main()