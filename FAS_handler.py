#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 10:00:07 2022

@author: chrisbl
"""

import json
import itertools

# fas.run -s /share/project/zarnack/chrisbl/FAS_library/homo_sapiens/release-107/buffer.fa -q query.fa -a /share/project/zarnack/chrisbl/FAS_library -o /share/project/zarnack/chrisbl/FAS_library/homo_sapiens/release-107/

TEST = {"A" : ["A", "B", "C", "D"], "B" : ["1", "2", "3", "4"]}

# class FAS_Memory:
#     def __init__(self, gene_ids, root_path):
        
        

def tsv_collection_maker(header_dict, root_path):
    tsv_dict = dict()
    for gene_id in header_dict.keys():
        tsv_dict[gene_id] = ""
        pairs = itertools.product(header_dict[gene_id], header_dict[gene_id])
        for x, y in pairs:
            tsv_dict[gene_id] += x + "\t" + y + "\n"
    with open(root_path + "pairings_tsv.json", 'w') as fp:
        json.dump(tsv_dict, fp,  indent=4)

def tsv_access(gene_id, root_path, output_name):
    with open(root_path + "pairings_tsv.json", "r") as fp: 
        tsv_data = json.load(fp)[gene_id]
    with open(root_path + output_name + ".tsv", "w") as fp:
        fp.write(tsv_data)
    

def main():
    """

    """
    tsv_collection_maker(TEST, "")
    tsv_access("A", "", "A")
    tsv_access("B", "", "B")
    # with open("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/test_annofas.json", "r") as f:
    #     fas = json.load(f)
    
    

if __name__ == "__main__":
    main()