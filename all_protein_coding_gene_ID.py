#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:12:20 2022

@author: chrisbl
"""

import pyensembl
import os

def load_ensembl_assembly(cache_dir, release_num):
    """
    

    Parameters
    ----------
    cache_dir : file path to cached human ensembl databases
    release_num : required release number

    Returns
    -------
    

    """
    os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir

    ENSEMBL = pyensembl.EnsemblRelease(release=release_num)
    
    print("Gathering protein coding genes...")
    protein_coding_gene_ids = [gene.id for gene in ENSEMBL.genes() if gene.biotype == "protein_coding"][0:3]
    
    print("Gathering transcripts...")
    # Get all sets of transcript_ids of protein coding genes.
    transcript_ids_list = [ [gene_id, ENSEMBL.transcript_ids_of_gene_id(gene_id)] for gene_id in protein_coding_gene_ids]
    
    print("Filtering all transcripts that are not protein coding...")
    # Remove all transcript IDs that belong to non-protein coding transcripts 
    # and assemble them in pairs with the gene_id
    transcript_list = list()
    transcript_dict = dict()
    for gene_id, transcript_ids in transcript_ids_list:
        transcript_dict[gene_id] = list()
        transcripts = [ENSEMBL.transcript_by_id(transcript_id) for transcript_id in transcript_ids]
        for transcript in transcripts:
            if transcript.biotype == "protein_coding":
                transcript_dict[gene_id].append([transcript.transcript_id])
                transcript_list.append(transcript.transcript_id)
    
    print("Sorting transcripts...")
    # Sort by 1. gene_id 2. transcript id
    transcript_list = sorted(transcript_list)
    for key in transcript_dict.keys():
        transcript_dict[key] = sorted(transcript_dict[key])

    return transcript_dict, transcript_list, protein_coding_gene_ids
    

if __name__ == "__main__":
    pass