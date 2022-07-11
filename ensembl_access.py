# -*- coding: utf-8 -*-
"""
__author__ = "Christian Bluemel"

"""

from Bio.Seq import Seq

import requests
import pyensembl
import argparse
import os
import subprocess

#Test command line:
# python ensembl_access.py -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/ -c /share/project/zarnack/chrisbl/FAS/utility/ensembl_cache

# /share/project/zarnack/chrisbl/FAS/utility/protein_lib

#some test genes
# TEST = ['ENSG00000197580',
#         'ENSG00000133884',
#         'ENSG00000184313',
#         'ENSG00000122707',
#         'ENSG00000183696',
#         'ENSG00000133794',
#         'ENSG00000168137',
#         'ENSG00000198467',
#         'ENSG00000167460',
#         'ENSG00000164093']

ENSEMBL_ASSMBLY = 106

def load_ensembl_assembly(cache_dir, release_num):    
    os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir

    ENSEMBL = pyensembl.EnsemblRelease(release=release_num)
    return ENSEMBL
    
def assemble_ensembl_ids(ensembl):
    """
    
    
    Returns
    ------

"""
        
    print("Gathering protein coding genes...")
    # Get all gene_ids of protein coding genes.
    protein_coding_gene_ids = [gene.id for gene in ensembl.genes() if gene.biotype == "protein_coding"]
    
    print("Gathering transcripts...")
    # Get all sets of transcript_ids of protein coding genes.
    transcript_ids_list = [ [gene_id, ensembl.transcript_ids_of_gene_id(gene_id)] for gene_id in protein_coding_gene_ids]

    print("Filtering all transcripts that are not protein coding...")
    # Remove all transcript IDs that belong to non-protein coding transcripts 
    # and assemble them in pairs with the gene_id
    transcript_list = list()
    transcript_dict = dict()
    for gene_id, transcript_ids in transcript_ids_list:
        transcript_dict[gene_id] = list()
        transcripts = [ensembl.transcript_by_id(transcript_id) for transcript_id in transcript_ids]
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
    url_suffix = "?object_type=transcript;type=cds;species=human"
    headers = { "Content-Type" : "text/plain"}
    
    print("Retrieving transcript sequences and translating them to protein sequences...")
    for key in transcript_dict.keys():
        for i, entry in enumerate(transcript_dict[key]):
            transcript_id = entry[0]
            sequence_retrieved = False
            attempt_count = 0
            while not sequence_retrieved:
                attempt_count += 1
                r = requests.get(url_prefix + transcript_id + url_suffix, headers=headers)
                if r.text[0] != "<":
                    sequence_retrieved = True
                if attempt_count > 30:
                    print("Could not retrieve sequence",
                          transcript_id,
                          "of gene", key,  "after 30 attempts. Aborting...")
                    print("Brace for TranslationError!")
                    sequence_retrieved = True
            seq = Seq(r.text)
            seq = seq.translate()
            if seq[-1] == "*":
                seq = seq[0:-1]
            transcript_dict[key][i].append(str(seq))
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
        
    ENSEMBL = load_ensembl_assembly(CACHE_DIR, ENSEMBL_ASSMBLY)
    
    transcript_dict, transcript_id_list, gene_id_list = assemble_ensembl_ids(ENSEMBL)
    
    library_dict = assemble_protein_seqs(transcript_dict)
    
    #fasta_entries = [">" + entry[0] + " " + entry[1] + "Ensembl Assemlby" + str(ENSEMBL_ASSMBLY)  + "\n" + entry[2] for entry in library]
    
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
    
#     import requests, sys
 
# server = "https://rest.ensembl.org"
# ext = "/xrefs/symbol/homo_sapiens/BRCA2?external_db=HGNC"
 
# r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
# if not r.ok:
#   r.raise_for_status()
#   sys.exit()
 
# decoded = r.json()
# print(repr(decoded))