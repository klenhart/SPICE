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

def load_ensembl_assembly(cache_dir, release_num):
    os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir

    ENSEMBL = pyensembl.EnsemblRelease(release=release_num)
    return ENSEMBL
    
def assemble_ensembl_ids(ensembl):
    """
    
    
    Returns
    ------
    2xn matrix of transcripts that are protein coding where n is the number 
    of transcripts and 2 counts the transcript and the gene id:
        [
            [transcript_id_1,gene_id_1],
            [transcript_id_2,gene_id_1],
            [transcript_id_3,gene_id_2],
        .
        .
        .
        ...etc..."""
        
    print("Gathering protein coding genes...")
    # Get all gene_ids of protein coding genes.
    protein_coding_gene_ids = [gene.id for gene in ensembl.genes() if gene.biotype == "protein_coding"]
    
    print("Gathering transcripts...")
    # Get all sets of transcript_ids of protein coding genes.
    transcript_ids_list = [ensembl.transcript_ids_of_gene_id(gene_id) for gene_id in protein_coding_gene_ids]

    library = list()
    
    print("Filtering all transcripts that are not protein coding and attaching them to the gene IDs...")
    # Remove all transcript IDs that belong to non-protein coding transcripts 
    # and assemble them in pairs with the gene_id
    for transcript_ids in transcript_ids_list:
        transcripts = [ensembl.transcript_by_id(transcript_id) for transcript_id in transcript_ids]
        protein_coding_id_pairs = [[transcript.transcript_id,
                                    transcript.gene_id] for transcript in transcripts if transcript.biotype == "protein_coding"]
        for pair in protein_coding_id_pairs:
            library.append(pair)
    
    print("Sorting assembled headers...")
    # Sort by 1. gene_id 2. transcript id
    library = (sorted(library, key=lambda x: (x[1], x[0])))
    
    return library


def assemble_protein_seqs(id_library):
    """

    Parameters
    ----------
    id_library : 2xn matrix with transcript id and gene id in each row

    Returns
    -------
    library: 3xn matrix with transcript id, gene id and protein sequence in
    each row

    """
    url_prefix = "https://rest.ensembl.org/sequence/id/"
    url_suffix = "?object_type=transcript;type=cds;species=human"
    headers = { "Content-Type" : "text/plain"}
    
    print("Retrieving transcript sequences and translating them to protein sequences...")
    for i, entry in enumerate(id_library):
        sequence_retrieved = False
        attempt_count = 0
        while not sequence_retrieved:
            attempt_count += 1
            r = requests.get(url_prefix + entry[0] + url_suffix, headers=headers)
            if r.text[0] != "<":
                sequence_retrieved = True
            if attempt_count > 30:
                print("Could not retrieve sequence " + i + " after 30 attempts. Aborting...")
                print("Brace for TranslationError!")
                sequence_retrieved = True
        seq = Seq(r.text)
        seq = seq.translate()
        if seq[-1] == "*":
            seq = seq[0:-1]
        id_library[i].append(str(seq))
        print(seq)
    return id_library
    
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
        
    ENSEMBL = load_ensembl_assembly(CACHE_DIR, 106)
    
    library = assemble_ensembl_ids(ENSEMBL)
    
    library = assemble_protein_seqs(library)
    
    library = [">" + entry[0] + " " + entry[1] + "\n" + entry[2] for entry in library]
    
    print("Saving library as fasta in " + OUTPUT_DIR + 'library.fasta')
    count = 1
    with open(OUTPUT_DIR + 'library.fasta', 'w') as fasta:
        count += 1
        for entry in library:
            fasta.write(entry)
            fasta.write('\n')
    print(count, "protein sequences assembled.")
    print("Library assembly complete.")

if __name__ == "__main__":
    main()