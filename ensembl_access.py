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


def assemble_protein_seqs(transcript_dict, assembly_num, species, library_path, root_path):
    """

    Parameters
    ----------
    dictionary with gene IDs as keys and transcript IDs in a list of one in another list.
        {
            gene_id_1 : [[transcript_id_1], [transcript_id_2]],
            gene_id_2 : [[transcript_id_3], [transcript_id_4]],
            .
            .
            .
            ...etc...
            
            }

    Returns
    -------
    library_dict with gene IDs as keys and transcript IDs combined with protein sequences in a list of two.
        {
            gene_id_1 : [[transcript_id_1, seq1], [transcript_id_2, seq2]],
            gene_id_2 : [[transcript_id_3, seq3], [transcript_id_4, seq4]],
            .
            .
            .
            ...etc...
            
            }
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
    
    count_not_found = 0
    count_found = 0
    count_genes = 0
    
    url_prefix = "https://rest.ensembl.org/sequence/id/"
    url_suffix = "?object_type=transcript;type=protein;species=" + species + ";"
    headers = { "Content-Type" : "application/json"}
    
    print(transcript_dict.keys(), " protein coding genes prepared for ensembl sequence requests...")
    
    for key in transcript_dict.keys():
        count_genes += 1
        for i, transcript_id in enumerate(transcript_dict[key]):
                for x in range(3):
                    r = requests.get(url_prefix + transcript_id + url_suffix, headers=headers)
                    if r.ok:
                        count_found += 1
                        seq = r.json()["seq"]
                        data = ">" + key + "|" + transcript_id + "|species:" + species + "|assembly:" + str(assembly_num) + "\n" + seq + "\n"
                        with open(isoforms_path, "a") as fasta:
                            fasta.write(data)
                        break
                    elif x == 2:
                        count_not_found += 1
                        with open(not_found_path, "a") as f:
                            f.write(key + "|" + transcript_id + "\n")
    return count_found, count_not_found, count_genes

def parser_setup():
    """
    

    Returns
    -------
    Output an Cache folders for the run.

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--output", type=str,
                        help=""""Specify location of the library. FAS_library folder can already exist in this folder.
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
        transcript_dict, transcript_id_list, gene_id_list = extract_protein_coding_ids(ensembl_path)
        if not os.path.exists(library_path):
            os.makedirs(library_path)
        root_path = make_rootpath(library_path, species, release_num) 
        count_found, count_not_found, count_genes = assemble_protein_seqs(transcript_dict, release_num, species, library_path, root_path)
        print("Saved isoforms as fasta in", root_path + "/isoforms.fasta")
        print(count_genes, "protein coding genes processed...")
        print(count_found, "protein sequences integrated into library assembled.")
        print(count_not_found, "protein sequences not found. IDs written into", root_path + "not_found.txt...")
        print("Library assembly complete.")

if __name__ == "__main__":
    main()
