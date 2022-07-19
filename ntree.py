#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:37:31 2022

@author: chrisbl
"""

import os

def make_rootpath(library_path, species, assembly_num):
    return library_path + species + "/assembly" + str(assembly_num) + "_GeneId000/"

def make_filepath(library_path, species, assembly_num, gene_id):
    rootpath = make_rootpath(library_path, species, assembly_num)
    return rootpath + gene_id[-8:-6] + "/" + gene_id[-6:-4] + "/" + gene_id[-4:-2] + "/" + gene_id[-2:] + "/isoform.fasta"
    

def append_to_leaf(library_path, assembly_num, gene_id, species, data):
    filepath = make_filepath(library_path, species, assembly_num, gene_id)
    with open(filepath, "a") as fasta:
        fasta.write(data)    

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
        
    assembly_path = species_path + "assembly" + str(assembly_num) + "_GeneId000/"
    if not os.path.exists(assembly_path):
        os.makedirs(assembly_path)
    
    first_level, second_level, third_level, fourth_level = cut_ids(gene_id_list)
    for entry in first_level:
        path = assembly_path + entry
        if not os.path.exists(path):
            os.makedirs(path)
    for entry in second_level:
        path = assembly_path + entry[0:2] + "/" + entry[2:]
        if not os.path.exists(path):
            os.makedirs(path)
    for entry in third_level:
        path = assembly_path + entry[0:2] + "/" + entry[2:4] + "/" + entry[4:]
        if not os.path.exists(path):
            os.makedirs(path)
    for entry in fourth_level:
        path = assembly_path + entry[0:2] + "/" + entry[2:4] + "/" + entry[4:6] + "/" + entry[6:] + "/"
        if not os.path.exists(path):
            os.makedirs(path)
            isoform_path = path + "isoform.fasta"        
            if not os.path.isfile(isoform_path):
                with open(isoform_path, 'w') as fp:
                    pass             

def main():
    # Just some test inputs...
    library_path = "/home/chrisbl/project/test/"
    species = "human"
    assembly_num = 107
    gene_id_list = ["ENST00000002165",
                    "ENST00000286031",
                    "ENST00000341376",
                    "ENST00000353205",
                    "ENST00000359326",
                    "ENST00000359637",
                    "ENST00000367429",
                    "ENST00000367770",
                    "ENST00000367771",
                    "ENST00000367772",
                    "ENST00000371582",
                    "ENST00000371584",
                    "ENST00000371588", #
                    "ENST00000373020",
                    "ENST00000373031",
                    "ENST00000374003",
                    "ENST00000374004",
                    "ENST00000374005",
                    "ENST00000399173",
                    "ENST00000413082",
                    "ENST00000413811",
                    "ENST00000423670",
                    "ENST00000451668",
                    "ENST00000457296",
                    "ENST00000472795",
                    "ENST00000496973",
                    "ENST00000505197",
                    "ENST00000509541",
                    "ENST00000513939",
                    "ENST00000514004",
                    "ENST00000612152",
                    "ENST00000614008",
                    "ENST00000616923",
                    "ENST00000630130",
                    "ENST00000643939",
                    "ENST00000650454",
                    "ENST00000683466"]
    generate_folder_tree(library_path, species, assembly_num, gene_id_list)
    gene_id = gene_id_list[12]
    seq = "ACGTCGATCGATATATCGATAGCTAGGGATATAGCTATATATCGATAGCGGC"
    data = ">ENSGsomethingsomething|" + gene_id + "|Evidence:3|species:" + species + "|assembly:" + str(assembly_num) + "\n" + seq + "\n"
    append_to_leaf(library_path, assembly_num, gene_id, species, data)

if __name__ == "__main__":
    main()
    
    