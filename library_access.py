#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 09:57:16 2022

@author: chrisbl
"""

# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

import os
import requests
import argparse
import sys

from ensembl_access import get_species_info
from ensembl_access import make_local_ensembl_name
from ensembl_access import get_release

#Test command line:
# python ensembl_access.py -s human -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib

def get_taxon_id(species): 
    server = "https://rest.ensembl.org"
    ext = "/info/genomes/taxonomy/" + species +"?"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    decoded = r.json()
    return decoded[0]["species_taxonomy_id"]

def make_rootpath(library_path, species, assembly_num):
    return library_path + species + "/release-" + str(assembly_num) + "/"

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
def ping_ensembl():
    import requests, sys
 
    server = "https://rest.ensembl.org"
    ext = "/info/ping?"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    decoded = r.json()
    return bool(decoded["ping"])

def make_request_data(id_list):
    request_data = '{ "ids" : ['
    for entry in id_list:
        request_data += '"' + entry + '", '
    request_data = request_data[:-2]
    request_data += ' ] }'
    return request_data

def extract_protein_coding_ids(ensembl_path):
    print("Loading local ensembl dataframe...")
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    print("Extracting IDs of protein coding transcripts...")
    protein_coding_ids = [ (gene_id, transcript_id, protein_id) for gene_id,
                                      protein_id,
                                      transcript_id,
                                      biotype in zip(gtf['gene_id'],
                                                    gtf['protein_id'],
                                                    gtf['transcript_id'],
                                                    gtf['transcript_biotype']
                                                    ) if biotype == "protein_coding" and str(type(protein_id)) != "<class 'float'>"]
    protein_coding_ids = list(set([i for i in protein_coding_ids]))
    print("Assembling IDs...")
    protein_coding_ids.sort()
    return protein_coding_ids

def assemble_protein_seqs(protein_coding_ids, assembly_num, species, library_path, root_path, taxon_id):
    """

    Parameters
    ----------
    list with gene and protein IDs as tuples
        [
            (gene_id_1, protein_id_1),
            (gene_id_1, protein_id_2),
            (gene_id_2, protein_id_3),
            .
            .
            .
            ...etc...
            
            ]

    Returns
    -------
    counts on how many request were found, not found and how many genes were assembled.
    """
    if not os.path.exists(root_path):
        os.makedirs(root_path)
    # not_found_path = root_path + "not_found.txt"
    # with open(not_found_path, 'w') as fp:
    #     pass
    isoforms_path = root_path + "isoforms.fasta"
    if not os.path.isfile(isoforms_path):
        with open(isoforms_path, "w") as fp:
            pass

    gene_ids = [gene_id for gene_id, protein_id in protein_coding_ids]
    count_genes = len(list(set(gene_ids)))
    
    protein_ids = [protein_id for gene_id, protein_id in protein_coding_ids]
    total_length = len(protein_ids)
    request_chunks = chunks(range(total_length), 50)
    ensembl_requests = []
    
    step = 0
    
    for i, chunk in enumerate(request_chunks):
        start = chunk[0]
        end = chunk[-1] + 1
        if end == total_length:
             ensembl_requests.append((step ,make_request_data(protein_ids[start:])))
        ensembl_requests.append((step, make_request_data(protein_ids[start:end])))
        step += 50

    server = "https://rest.ensembl.org/sequence/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}    
    
    header_dict = dict()
    for gene_id in gene_ids:
        header_dict[gene_id] = []
    
    for step, request in ensembl_requests:
        for x in range(3):
            r = requests.post(server, headers=headers, data=request)
            if r.ok:
                break
            elif x > 1:
                print("Failed to request sequenes 3 times. Checking if ensembl service is up. Step:", step)
                if not ping_ensembl():
                    print("Ensembl is currently down. Can't download sequences.")
                else:
                    print("Ensembl is up. Weird...")
                r.raise_for_status()
                sys.exit()
        decoded = r.json()
        id_seq_tuple_list = [(entry["query"], entry["seq"]) for entry in decoded]
        for i, id_seq_tuple in enumerate(id_seq_tuple_list):
            id_complement = i + step
            query_id, seq = id_seq_tuple
            gene_id, protein_id = protein_coding_ids[id_complement]
            if query_id != protein_id:
                print(protein_id,  "does not match", query_id, "Step was", step, "and index was", i)
                sys.exit()
            else:
                header = gene_id + "|" + protein_id + "|" + str(taxon_id)
                header_dict[gene_id].append(header)
                with open(isoforms_path, "a") as fasta:
                    fasta.write(">" + header + "\n" + seq + "\n")
    return header_dict, count_genes

def collect_data_paths(data_path):
    with open(data_path, "r") as f:
        data = f.read()
        data_paths = data.split("\n")
    return data_paths

def parser_setup():
    """
    Returns
    -------
    parent directory of library,
    output directory,
    species of expression data,
    path to textfile containing path to all expression files.    

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-l", "--library", type=str,
                        help="""Specify path of the parent folder of the library.""")

    parser.add_argument("-s", "--species", type=str,
                        help="Specify the species of the expression data.")
    
    parser.add_argument("-o", "--output", type=str,
                        help="Specify the path for the output.")
    
    parser.add_argument("-d", "--data", type=str,
                        help="Specify the path to a text file containing all paths to the expression data. One path per row.")

    args = parser.parse_args()
    
    library_path = args.library
    output_path = args.output
    species = args.species
    data_path = args.data

    return library_path, output_path, species, data_path

def main():
    """
    Returns
    -------
    Matches expression data with pre-calculated FAS Scores in library.
    """
    #Collect input
    library_path, output_path, species, data_path = parser_setup()
    
    #Collect input paths.
    data_paths = collect_data_paths(data_path)
    
    #Get ensembl stats for naming conventions
    species, url_name, assembly_default = get_species_info(species)
    release_num = get_release()
    
    #Make paths
    ensembl_path = make_local_ensembl_name(library_path,
                                           release_num,
                                           species,
                                           ".gtf",
                                           assembly_default,
                                           url_name)
    root_path = make_rootpath(library_path,
                              species,
                              release_num)
    
    #Collect IDs of protein coding transcripts.
    protein_coding_ids = extract_protein_coding_ids(ensembl_path)
    for path in data_paths:
        extract_ids_from_expression_data(path, protein_coding_ids)
        
    
    
    
if __name__ == "__main__":
    main()

