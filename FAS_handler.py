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

# /home/chrisbl/miniconda3/envs/FAS/bin/fas.doAnno
# /home/chrisbl/miniconda3/envs/FAS/bin/fas.run

TEST = {"A" : ["A", "B", "C", "D"], "B" : ["1", "2", "3", "4"]}

def bash_command_maker(root_path, python_path, FAS_handler_path, fas_path):
    gene_ids_path = root_path + "gene_ids.txt"
    isoforms_path = root_path + "isoforms.fasta"
    phyloprofile_path = root_path + "phyloprofile_ids.tsv"
    buffer_path = root_path + "tsv_buffer/"
    
    with open(gene_ids_path, "r") as txt:
        gene_ids = txt.read()
    gene_ids = gene_ids.split("\n")
    
    command_list = []
    for gene_id in gene_ids:
        access_command, remove_command = pairings_command_maker(gene_id,
                                                                python_path,
                                                                FAS_handler_path,
                                                                root_path)
        pairings_tsv_path = buffer_path + gene_id + ".tsv"
        FAS_command = FAS_command_maker(gene_id,
                                        isoforms_path,
                                        phyloprofile_path,
                                        pairings_tsv_path,
                                        fas_path,
                                        root_path)
        full_command = access_command + " && " + FAS_command + " ; " + remove_command
        command_list.append(full_command)
    command_str = command_list.join("\n")
    with open(root_path + "commands.txt", "w") as commands:
        commands.write(command_str)
        

def pairings_command_maker(gene_id, python_path, FAS_handler_path, root_path):
    options_access = []
    options_remove = []
    options_access.append(python_path)
    options_remove.append(python_path)
    options_access.append(FAS_handler_path)
    options_remove.append(FAS_handler_path)
    options_access.append("-m")
    options_remove.append("-r")
    options_access.append("-g " + gene_id)
    options_remove.append("-g " + gene_id)
    options_access.append("-o " + root_path)
    options_remove.append("-o " + root_path)
    access_command = options_access.join(" ")
    remove_command = options_remove.join(" ")
    return access_command, remove_command


def FAS_command_maker(gene_id, isoforms_path, phyloprofile_path, pairings_tsv_path, fas_path, root_path):
    options = []
    options.append(fas_path)
    options.append("--seed " + isoforms_path)
    options.append("--query " + isoforms_path)
    options.append("--annotation_dir "+ root_path)
    options.append("--out_dir " + root_path + "FAS_buffer/")
    options.append("--bidirectional")
    options.append("--pairwise " + pairings_tsv_path)
    options.append("--out_name " + gene_id)
    options.append("--tsv")
    options.append("--phyloprofile " + phyloprofile_path)
    options.append("--domain")
    options.append("--empty_as_1")
    fas_command = options.join(" ")
    return fas_command

    

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
                        help="""Specify the directory that contains the pairings_tsv.json and should contain the output.""")
    
    parser.add_argument("-p", "--python", type=str, default=None,
                        help="""Specify the location of the pythonversion to use. This is only necessary when using the --bash option.""")

    parser.add_argument("-h", "--handler", type=str, default=None,
                        help="""Specify the location FAS_handler.py. This is only necessary when using the --bash option.""")
    
    parser.add_argument("-f", "--fas", type=str, default=None,
                        help="""Specify the location of fas.run. This is only necessary when using the --bash option.""")
    
    parser.add_argument("-g", "--gene", type=str, default=None
                        help="""Ensembl gene ID of the gene that the tsv file should be created or deleted for.""")
    
    parser.add_argument("-r", "--remove", action="store_true",
                        help="Remove the tsv file of the given gene.")

    parser.add_argument("-m", "--maketsv", action="store_true",
                        help="Make the tsv file for the given gene.")
    
    parser.add_argument("-b", "--bash", action="store_true",
                        help="Make one FAS bash command for each gene in a text file.")
    

    args = parser.parse_args()

    root_path = args.output
    gene_id = args.gene
    python_path = args.python
    FAS_handler_path = args.handler
    fas_path = args.fas
    flag_remove = args.local
    flag_maketsv = args.output
    flag_bash = args.bash

    return root_path, python_path, FAS_handler_path, fas_path, gene_id, flag_remove, flag_maketsv, flag_bash

def main():
    """
    Create tsv output for FAS or delete a tsv that was used by FAS already.
    """
    root_path, python_path, FAS_handler_path, fas_path, gene_id, flag_remove, flag_maketsv, flag_bash = parser_setup()
    if flag_maketsv:
        tsv_access(gene_id, root_path)
    elif flag_remove:
        tsv_remove(gene_id, root_path)
    elif flag_bash:
        bash_command_maker(root_path, python_path, FAS_handler_path, fas_path)
        
    
    
    
    

if __name__ == "__main__":
    main()