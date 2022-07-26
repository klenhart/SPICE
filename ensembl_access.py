# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

import os
import requests
import argparse
import sys

from ntree import make_rootpath

from all_protein_coding_gene_ID import extract_protein_coding_ids

from install_local_ensembl import get_release
from install_local_ensembl import install_local_ensembl
from install_local_ensembl import get_species_info
from install_local_ensembl import make_local_ensembl_name

from healthcheck import healthcheck


#Test command line:
# python ensembl_access.py -s human -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib

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

def assemble_protein_seqs(protein_coding_ids, assembly_num, species, library_path, root_path):
    """

    Parameters
    ----------
    list with gene, transcript and protein IDs as tuples
        [
            (gene_id_1, transcript_id_1, protein_id_1),
            (gene_id_1, transcript_id_2, protein_id_2),
            (gene_id_2, transcript_id_3, protein_id_3),
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
    not_found_path = root_path + "not_found.txt"
    with open(not_found_path, 'w') as fp:
        pass
    isoforms_path = root_path + "isoforms.fasta"
    if not os.path.isfile(isoforms_path):
        with open(isoforms_path, "w") as fp:
            pass

    gene_ids = [gene_id for gene_id, transcript_id, protein_id in protein_coding_ids]
    count_genes = len(list(set(gene_ids)))
    
    protein_ids = [protein_id for gene_id, transcript_id, protein_id in protein_coding_ids]
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
            gene_id, transcript_id, protein_id = protein_coding_ids[id_complement]
            if query_id != protein_id:
                print(protein_id,  "does not match", query_id, "Step was", step, "and index was", i)
                sys.exit()
            else:
                data = ">" + gene_id + "|" + transcript_id + "|" + protein_id + "|species:" + species + "|assembly:" + str(assembly_num) + "\n" + seq + "\n"
                with open(isoforms_path, "a") as fasta:
                    fasta.write(data)

    server = "https://rest.ensembl.org/sequence/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    request_data = '{ "ids" : ["ENSTNIP00000003408", "ENSTNIP00000003894", "ENSTNIP00000002623", "ENSTNIP00000000936", "ENSTNIP00000000289", "ENSTNIP00000001327", "ENSTNIP00000002520", "ENSTNIP00000000812" ] }'
    r = requests.post(server, headers=headers, data=request_data)

def parser_setup():
    """
    

    Returns
    -------
    Output an Cache folders for the run.

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--output", type=str,
                        help="""Specify location of the library. FAS_library folder can already exist in this folder.
                        If creating the library, this folder should contain the local ensembl assembly. If it does not,
                        you can download the local ensembl assembly using the -l (--local) argument. It will then be automatically
                        download into the FAS_library folder within the folder given in this argument.""")

    parser.add_argument("-s", "--species", type=str, default=None,
                        help="Specify the species.")
    
    parser.add_argument("-l", "--local", action="store_true",
                        help="Download and unpack a local ensembl assembly into the folder defined in --output.")
    
    parser.add_argument("-k", "--healthcheck", action='store_true', #todo
                        help="""Checks several parameters of the general state of the library
                        but also the library specified in via the species and assembly arguments.""")

    args = parser.parse_args()

    output = args.output
    species = args.species
    flag_install_local = args.local
    flag_healthcheck = args.healthcheck

    return output, species, flag_install_local, flag_healthcheck

def main():
    """
    Returns
    -------
    fasta that contains the transcript_id, the gene_id as a header and the
    consensus sequence like this:
        >trancript_id1 gene_id1
        consensus_sequence1
        >trancript_id2 gene_id1
        consensus_sequence2
        >trancript_id3 gene_id2
        consensus_sequence3
        .
        .
        .
        ...etc...
    """
    OUTPUT_DIR, species, flag_install_local, flag_healthcheck = parser_setup()  
    library_path = OUTPUT_DIR + "/FAS_library/"
    release_num = get_release()
    species, url_name, assembly_default = get_species_info(species)
    ensembl_path = make_local_ensembl_name(library_path, release_num, species, ".gtf", assembly_default, url_name)
  
    if flag_install_local:
        print("Local ensembl installation commencing...")
        install_local_ensembl(species, release_num, library_path, url_name, assembly_default)
    
    elif flag_healthcheck:
        print("Library healthcheck commencing...")
        healthcheck(library_path, species, release_num)
    else:
        print("Library generation commencing...")
        if not os.path.isfile(ensembl_path):
            print(ensembl_path, "does not exist. Maybe you have an old release of the local ensembl GTF. You can download a current one for your species by using the -l argument additional to your currently used arguments.")
            sys.exit()
        protein_coding_ids = extract_protein_coding_ids(ensembl_path)
        if not os.path.exists(library_path):
            os.makedirs(library_path)
        root_path = make_rootpath(library_path, species, release_num) 
        assemble_protein_seqs(protein_coding_ids, release_num, species, library_path, root_path)
        print("Saved isoforms as fasta in", root_path + "/isoforms.fasta")
        print(count_genes, "protein coding genes processed.")
        print(count_found, "protein sequences integrated into library assembled.")
        print(count_not_found, "protein sequences not found. IDs written into", root_path + "not_found.txt.")
        print("Library assembly complete.")

if __name__ == "__main__":
    main()

