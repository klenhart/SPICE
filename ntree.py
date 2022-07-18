#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:37:31 2022

@author: chrisbl
"""

import os
import sys

def cut_ids(id_list):
    pairs = [ (entry[-8:-6], entry[-8:-4], entry[-8:-2], entry[-8:]) for entry in id_list ]
    first_level = set()
    second_level = set()
    third_level = set()
    fourth_level = set()
    for first, second, third, fourth in pairs:
        first_level.add(first)
        second_level.add(second)
        third_level.add(third)
        fourth_level.add(fourth)
    first_level = list(first_level)
    second_level = list(second_level)
    third_level = list(third_level)
    fourth_level = list(fourth_level)
    
    return first_level, second_level, third_level, fourth_level

def generate_folder_tree(library_path, species, assembly_num, gene_id_list):
    if not os.path.exists(library_path):
        os.makedirs(library_path)
        
    species_path = library_path + species + "/"
    if not os.path.exists(species_path):
        os.makedirs(species_path)
        
    assembly_path = species_path + "/" + str(assembly_num)
    if not os.path.exists(assembly_path):
        os.makedirs(assembly_path)
    
    first_level, second_level, third_level, fourth_level = cut_ids(gene_id_list)
    for entry in first_level:
        pass
    for entry in second_level:
        pass
    for entry in third_level:
        pass
    for entry in fourth_level:
        pass
    
    
    

def main():
    generate_folder_tree(library_path, species, assembly_num, gene_id_list)

if __name__ == "__main__":
    main()
    
    

    # library_path = OUTPUT_DIR + "/library/"
    # if not os.path.exists(library_path):
    #     os.makedirs(library_path)
    # for gene in gene_id_list:
    #     gene_path = library_path + gene
    #     if not os.path.exists(gene_path):
    #         os.makedirs(gene_path)
    #         file_path = gene_path + "/isoform.fasta"
    #         with open(file_path, 'w') as fp:
    #             pass