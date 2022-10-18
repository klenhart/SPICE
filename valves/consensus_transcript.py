#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 10:50:12 2022

@author: chrisbl
"""

import requests, sys

from valves.fas_utility import make_request_data
from valves.fas_utility import tsv_to_tuple_list
from valves.fas_utility import chunks
from valves.fas_utility import build_fasta_from_ids_and_fasta
from valves.fas_utility import transcript_id_to_protein_id

def make_fasta_of_canonical_transcript_ids(fas_lib):
    all_ids_path = fas_lib.get_config("protein_coding_ids_path") 
    id_triples = tsv_to_tuple_list(all_ids_path)
    gene_chunks = chunks(id_triples, 1000)
    
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    
    canonical_path = fas_lib.get_config("canonical_path")
    isoforms_path = fas_lib.get_config("isoforms_path")
    
    
    transcript_id_list = []
    protein_id_list = []
    
    for chunk in gene_chunks:
        count = 0
        request = make_request_data(chunk, 0) 
        while count < 3:
            r = requests.post(server+ext, headers=headers, data=request)
            if r.ok:
                break
            elif count == 2:
                print("Failed to request ids.")
                sys.exit()
            else:
                count += 1
        decoded = r.json()
        for key in decoded.keys():
            canonical_id = decoded[key]["canonical_transcript"]
            canonical_id = canonical_id.split(".")[0]
            transcript_id_list.append(canonical_id)
    protein_id_list = transcript_id_to_protein_id(transcript_id_list, all_ids_path)
    build_fasta_from_ids_and_fasta(protein_id_list, isoforms_path, canonical_path)
    

def main():
    pass

