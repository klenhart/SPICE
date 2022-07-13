# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

import requests
import argparse
import os
import sys
import subprocess
from all_protein_coding_gene_ID import load_ensembl_assembly




#Test command line:
# python ensembl_access.py -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/ -c /share/project/zarnack/chrisbl/FAS/utility/ensembl_cache

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib

ENSEMBL_ASSMBLY = 106

def assemble_protein_seqs(transcript_dict):
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
    url_prefix = "https://rest.ensembl.org/sequence/id/"
    url_suffix = "?object_type=transcript;type=protein;species=human;"
    headers = { "Content-Type" : "application/json"}
    
    print("Retrieving protein sequences...")
    for key in transcript_dict.keys():
        for i, entry in enumerate(transcript_dict[key]):
            transcript_id = entry[0]
            sequence_retrieved = False
            attempt_count = 0
            while not sequence_retrieved:
                attempt_count += 1
                r = requests.get(url_prefix + transcript_id + url_suffix, headers=headers)
                if r.ok:
                    sequence_retrieved = True
                    seq = r.json()["seq"]
                    transcript_dict[key][i].append(seq)
                if attempt_count > 30:
                    r.raise_for_status()
                    sys.exit()
    return transcript_dict

    
def parser_setup():
    """
    

    Returns
    -------
    Output an Cache folders for the run.

    """
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-s", "--setup", action="store_true",
                        help="If the ensembl human assembly needs to be downloaded beforehand.")
    
    parser.add_argument("-c", "--cache", type=str,
                        help="specify location ensembl was cached in.")
    
    parser.add_argument("-o", "--output", type=str,
                        help="specify output folder.")

    args = parser.parse_args()
    
    cache = args.cache
    output = args.output
    install_pyensembl = args.setup
    
    return output, cache, install_pyensembl

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
    OUTPUT_DIR, CACHE_DIR, install_pyensembl = parser_setup()
    
    if install_pyensembl:
        subprocess.run(["export", "PYENSEMBL_CACHE_DIR=" + CACHE_DIR])
        subprocess.run(["pyensembl", "install", "--release 106", "--species human"])
        
    transcript_dict, transcript_id_list, gene_id_list = load_ensembl_assembly(CACHE_DIR,
                                                                              ENSEMBL_ASSMBLY)    
    library_dict = assemble_protein_seqs(transcript_dict)
    
    print("Generating subfolders in /library/[gene_id]...")
    newpath = OUTPUT_DIR + "/library/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    for gene in gene_id_list:
        gene_path = newpath + gene
        if not os.path.exists(gene_path):
            os.makedirs(gene_path)

    print("Saving transcript IDs in ", newpath)
    count = 0
    with open(newpath + "transcript_ids.txt", "w") as fasta:
        for transcript_id in transcript_id_list:
            fasta.write(transcript_id)
            fasta.write("\n")
            count += 1

    print("Saving isoforms as fasta in", newpath + "[gene_id]/isoforms.fasta")
    for key in library_dict.keys():
        with open(newpath + key + "/isoform.fasta", "w") as fasta:
            for entry in library_dict[key]:
                fasta.write("> " + entry[0] + " " + key + " Ensembl Assembly " + str(ENSEMBL_ASSMBLY) + "\n" + entry[1])
                fasta.write("\n")

    print(count, "protein sequences assembled.")
    print("Library assembly complete.")

if __name__ == "__main__":
    main()
