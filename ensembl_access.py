# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

import requests
import argparse
import os
import sys
from all_protein_coding_gene_ID import extract_protein_coding_ids
from check_library import check_for_transcript




#Test command line:
# python ensembl_access.py -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/ -g /share/project/zarnack/chrisbl/FAS/utility/protein_lib/Homo_sapiens.GRCh38.107.gtf

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib

ENSEMBL_ASSMBLY = 107

def assemble_protein_seqs(transcript_dict, library_path):
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
    not_found_path = library_path + "not_found.txt"
    transcript_list_path = library_path + "transcript_list.txt"
    with open(not_found_path, 'w') as fp:
        pass
    if not os.path.isfile(transcript_list_path): 
        with open(transcript_list_path, 'w') as fp:
            pass
    
    count_not_found = 0
    count_found = 0
    count_already = 0
    
    url_prefix = "https://rest.ensembl.org/sequence/id/"
    url_suffix = "?object_type=transcript;type=protein;species=human;"
    headers = { "Content-Type" : "application/json"}
    
    print("Retrieving protein sequences...")
    for key in transcript_dict.keys():
        for i, transcript_id in enumerate(transcript_dict[key]):
            flag_already_loaded = check_for_transcript(key, transcript_id, library_path)
            if flag_already_loaded:
                count_already += 1
                continue
            else:
                for x in range(3):
                    r = requests.get(url_prefix + transcript_id + url_suffix, headers=headers)
                    if r.ok:
                        count_found += 1
                        seq = r.json()["seq"]
                        with open(transcript_list_path, "a") as f:
                            f.write(transcript_id + "\n")
                        with open(library_path + key + "/isoform.fasta", "a") as fasta:
                            fasta.write("> " + transcript_id + " " + key + " EnsemblAssembly 107\n" + seq + "\n")
                        break
                    elif x == 2:
                        count_not_found += 1
                        with open(not_found_path, "a") as f:
                            f.write(key + " " + transcript_id + "\n")
    return count_found, count_not_found, count_already

def parser_setup():
    """
    

    Returns
    -------
    Output an Cache folders for the run.

    """
    #Setting up parser:
    parser = argparse.ArgumentParser()
        
    parser.add_argument("-g", "--ensemblgtf", type=str,
                        help="specify the path of the ensembl human genome gtf.")
    
    parser.add_argument("-o", "--output", type=str,
                        help="specify output folder.")

    args = parser.parse_args()
    
    path = args.ensemblgtf
    output = args.output
    
    return output, path

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
    OUTPUT_DIR, ensembl_path = parser_setup()
        
    transcript_dict, transcript_id_list, gene_id_list = extract_protein_coding_ids(ensembl_path)    
    
    print("Generating subfolders in /library/[gene_id]...")
    library_path = OUTPUT_DIR + "/library/"
    if not os.path.exists(library_path):
        os.makedirs(library_path)
    for gene in gene_id_list:
        gene_path = library_path + gene
        if not os.path.exists(gene_path):
            os.makedirs(gene_path)
            file_path = gene_path + "/isoform.fasta"
            with open(file_path, 'w') as fp:
                pass

    count_found, count_not_found, count_already = assemble_protein_seqs(transcript_dict, library_path)    
    print("Saved transcript IDs in ", library_path)
    print("Saved isoforms as fasta in", library_path + "[gene_id]/isoforms.fasta")
    print(count_found, "protein sequences integrated into library assembled.")
    print(count_not_found, "protein sequences not found. IDs written into", library_path + "not_found.txt...")
    print(count_already, "protein sequences already found in library and were not download again.")
    print("Library assembly complete.")

if __name__ == "__main__":
    main()
