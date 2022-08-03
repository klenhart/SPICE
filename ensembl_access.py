# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

import os
import requests
import argparse
import sys

from all_protein_coding_gene_ID import extract_protein_coding_ids

from install_local_ensembl import get_release
from install_local_ensembl import install_local_ensembl
from install_local_ensembl import get_species_info
from install_local_ensembl import make_local_ensembl_name
from FAS_handler import tsv_collection_maker

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
    
    for i in range(3):
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
        if r.ok:
            break
        elif i >= 2:
            print("Ensembl is currently down. Can't download sequences. Aborting...")
            r.raise_for_status()
            sys.exit() 
    decoded = r.json()
    return bool(decoded["ping"])

def make_request_data(id_list):
    id_list = [protein_id for gene_id, protein_id, transcript_id in id_list]
    request_data = '{ "ids" : ['
    for entry in id_list:
        request_data += '"' + entry + '", '
    request_data = request_data[:-2]
    request_data += ' ] }'
    return request_data

def make_folders_and_files(root_path):
    if not os.path.exists(root_path):
        os.makedirs(root_path)
    if not os.path.exists(root_path + "tsv_buffer"):
        os.makedirs(root_path + "tsv_buffer")
    if not os.path.exists(root_path + "FAS_buffer"):
        os.makedirs(root_path + "FAS_buffer")
    isoforms_path = root_path + "isoforms.fasta"
    if not os.path.isfile(isoforms_path):
        with open(isoforms_path, "w") as fp:
            pass
    phyloprofile_ids_path = root_path + "phyloprofile_ids.tsv"
    if not os.path.isfile(phyloprofile_ids_path):
        with open(phyloprofile_ids_path, "w") as fp:
            pass
    gene_ids_path = root_path + "gene_ids.txt"
    if not os.path.isfile(gene_ids_path):
        with open(gene_ids_path, "w") as fp:
            pass
    return gene_ids_path, isoforms_path, phyloprofile_ids_path

def extract_gene_ids(protein_coding_ids):
    gene_ids = sorted(list(set([gene_id for gene_id, protein_id, transcript_id in protein_coding_ids])))
    count = len(gene_ids)
    return count, gene_ids

def save_gene_ids_txt(gene_ids, gene_ids_path):
    gene_id_str = "\n".join(["gene"] + gene_ids)
    with open(gene_ids_path, "w") as file:
        file.write(gene_id_str)

def assemble_protein_seqs(protein_coding_ids, assembly_num, species, library_path, root_path, taxon_id):
    """

    Parameters
    ----------
    list with gene and protein IDs as tuples
        [
            (gene_id_1, protein_id_1, transcript_id_1),
            (gene_id_1, protein_id_2, transcript_id_2),
            (gene_id_2, protein_id_3, transcript_id_3),
            .
            .
            .
            ...etc...
            
            ]

    Returns
    -------
    counts on how many request were found, not found and how many genes were assembled.
    """
    gene_ids_path, isoforms_path, phyloprofile_ids_path = make_folders_and_files(root_path)

    count_genes, gene_ids = extract_gene_ids(protein_coding_ids)
    
    save_gene_ids_txt(gene_ids, gene_ids_path)

    request_chunks = list(chunks(protein_coding_ids, 50))
    ensembl_requests = []
    
    print(len(request_chunks))
    
    for i, chunk in enumerate(request_chunks):
        ensembl_requests.append(make_request_data(chunk))

    server = "https://rest.ensembl.org/sequence/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}    
    
    header_dict = dict()
    for gene_id in gene_ids:
        header_dict[gene_id] = []
    
    for i, request in enumerate(ensembl_requests):
        ids_list = list(request_chunks)[i]
        for x in range(3):
            r = requests.post(server, headers=headers, data=request)
            if r.ok:
                break
            elif x > 1:
                r.raise_for_status()
                sys.exit()
        decoded = r.json()
        id_seq_tuple_list = [(entry["query"], entry["seq"]) for entry in decoded]
        for j, id_seq_tuple in enumerate(id_seq_tuple_list):
            gene_id, protein_id, transcript_id = ids_list[j]
            query_id, seq = id_seq_tuple
            header = gene_id + "|" + protein_id + "|" + str(taxon_id)
            header_dict[gene_id].append(header)
            with open(phyloprofile_ids_path, "a") as file:
                file.write(header + "\t" + str(taxon_id) + "\n")
            with open(isoforms_path, "a") as fasta:
                fasta.write(">" + header + "\n" + seq + "\n")
    return header_dict, count_genes

def check_isoforms(isoforms_path):
    with open(isoforms_path, "r") as f:
        fasta = f.read()
    fasta_lines = fasta.split("\n")
    fasta_lines = fasta_lines[:-1]
    header = fasta_lines[::2]
    gene_ids = [ entry.split("|")[0][1:] for entry in header ]
    protein_ids = [ entry.split("|")[1] for entry in header ]
    
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    ext2 = "/lookup/id/"

    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    prefix_request = '{ "ids" : ["' 
    infix_request = '", "'
    suffix_request ='" ] }'
    for i, protein_id in enumerate(protein_ids): 
        
        gene_id = gene_ids[i]
        request = prefix_request + protein_id + infix_request + gene_id + suffix_request
        r = requests.post(server+ext, headers=headers, data=request)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        gene_r= decoded[gene_id]
        flag_biotype = gene_r["biotype"].startswith("protein_coding")
        
        protein_r = decoded[protein_id]

        transcript_id = protein_r["Parent"]
        r = requests.get(server+ext2+transcript_id+"?", headers={ "Content-Type" : "application/json"})
        
        transcript_r = r.json()
        flag_same_id = transcript_r["Parent"] == gene_id
        flag_biotype2 = transcript_r["biotype"].startswith("protein_coding")
        if not all([flag_same_id, flag_biotype2, flag_biotype]):
            print("Problem found!")
            print(gene_id, transcript_id, protein_id)

def ensembl_access(OUTPUT_DIR, species, flag_install_local):
    """
    Returns
    -------
    fasta called isoform.fasta that contains the the gene_id, protein_id, taxon_id as a header and the
    consensus sequence like this:
        >gene_id1|protein_id1|taxon_id
        consensus_sequence1
        >gene_id1|protein_id2|taxon_id
        consensus_sequence2
        >gene_id2|protein_id3|taxon_id
        consensus_sequence3
        .
        .
        .
        ...etc...
    """
    ping_ensembl()
    library_path = OUTPUT_DIR + "/FAS_library/"
    release_num = get_release()
    species, url_name, assembly_default = get_species_info(species)
    taxon_id = get_taxon_id(species)
    ensembl_path = make_local_ensembl_name(library_path, release_num, species, ".gtf", assembly_default, url_name)
  
    if flag_install_local:
        print("Local ensembl installation commencing...")
        install_local_ensembl(species, release_num, library_path, url_name, assembly_default)

    else:
        print("Library generation commencing...")
        if not os.path.isfile(ensembl_path):
            print(ensembl_path, "does not exist. Maybe you have an old release of the local ensembl GTF. You can download a current one for your species by using the -l argument additional to your currently used arguments.")
            sys.exit()
        protein_coding_ids = extract_protein_coding_ids(ensembl_path)
        if not os.path.exists(library_path):
            os.makedirs(library_path)
        root_path = make_rootpath(library_path, species, release_num) 
        header_dict, count_genes = assemble_protein_seqs(protein_coding_ids, release_num, species, library_path, root_path, taxon_id)
        tsv_collection_maker(header_dict, root_path)
        print(count_genes, "genes assembled.")
        print("Saved isoforms as fasta in", root_path + "/isoforms.fasta")
        print("Library assembly complete.")

def test():
    species = "human"
    release_num = 107
    species, url_name, assembly_default = ('homo_sapiens', 'Homo_sapiens', 'GRCh38')
    library_path = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/"
    ensembl_path = make_local_ensembl_name(library_path, release_num, species, ".gtf", assembly_default, url_name)
    return extract_protein_coding_ids(ensembl_path)
    
def main():
    species = "human"
    # release_num = 107
    # species, url_name, assembly_default = ('homo_sapiens', 'Homo_sapiens', 'GRCh38')
    # library_path = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/"
    # ensembl_path = make_local_ensembl_name(library_path, release_num, species, ".gtf", assembly_default, url_name)
    
if __name__ == "__main__":
    main()

