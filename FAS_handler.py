#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 10:00:07 2022

@author: chrisbl
"""

import json
import itertools
import argparse
import os

# fas.run -s /share/project/zarnack/chrisbl/FAS_library/homo_sapiens/release-107/buffer.fa -q query.fa -a /share/project/zarnack/chrisbl/FAS_library -o /share/project/zarnack/chrisbl/FAS_library/homo_sapiens/release-107/

TEST = {"A" : ["A", "B", "C", "D"], "B" : ["1", "2", "3", "4"]}

# class FAS_Memory:
#     def __init__(self, gene_ids, root_path):
        

def tsv_collection_maker(header_dict, root_path):
    tsv_dict = dict()
    for gene_id in header_dict.keys():
        tsv_dict[gene_id] = ""
        pairs = itertools.product(header_dict[gene_id], header_dict[gene_id])
        pairs = [pair for pair in pairs if pair[0] < pair[1]]
        for header_1, header_2 in pairs:
            tsv_dict[gene_id] += header_1 + "\t" + header_2 + "\n"
    with open(root_path + "pairings_tsv.json", 'w') as fp:
        json.dump(tsv_dict, fp,  indent=4)

def tsv_access(gene_id, root_path):
    buffer_path = root_path + "tsv_buffer/"
    with open(root_path + "pairings_tsv.json", "r") as fp: 
        tsv_data = json.load(fp)[gene_id]
    with open(buffer_path + gene_id + ".tsv", "w") as fp:
        fp.write(tsv_data)

def tsv_remove(gene_id, root_path):
    buffer_path = root_path + "tsv_buffer/"
    os.remove(buffer_path + gene_id + ".tsv")

def parser_setup():
    """
    

    Returns
    -------
    Get run options if either to create a tsv file, to delete it again or to concatenate FAS results into a main file.
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--output", type=str,
                        help="""Specify directory of where the tsv should be deleted or be found.
                        This directory should also contain the pairings_tsv.json""")
    
    parser.add_argument("-g", "--gene", type=str,
                        help="""Ensembl gene ID of the gene that the tsv file should be created or deleted for.""")
    
    parser.add_argument("-r", "--remove", action="store_true",
                        help="Remove the tsv file of the given gene.")

    parser.add_argument("-m", "--maketsv", action="store_true",
                        help="Make the tsv file for the given gene.")

    args = parser.parse_args()

    root_path = args.output
    gene_id = args.gene
    flag_remove = args.local
    flag_maketsv = args.output

    return root_path, gene_id, flag_remove, flag_maketsv

def main():
    """
    Create tsv output for FAS or delete a tsv that was used by FAS already.
    """
    root_path, gene_id, flag_remove, flag_maketsv = parser_setup()
    if flag_maketsv:
        tsv_access(gene_id, root_path)
    elif flag_remove:
        tsv_remove(gene_id, root_path)
    
    
    
    

if __name__ == "__main__":
    main()