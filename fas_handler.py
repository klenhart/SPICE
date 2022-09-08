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

from library_class import Library

def concat_FAS_output(fas_lib):
    distance_master_path = fas_lib.get_config("distance_master_path")
    fas_buffer_path = fas_lib.get_config("fas_buffer_path")
    file_names = os.listdir(fas_buffer_path)
    print(len(file_names))
    print(fas_buffer_path)
    for filename in file_names:
        path = fas_buffer_path + filename
        with open(path, "r") as f_in:
            query = f_in.read().split("\n")
            query = query[1:]
            if len(query[-1]) == 0:
                query = query[0:-1]
            query = "\n".join(query)
            with open(distance_master_path, "a") as f_out:
                f_out.write(query + "\n")
        os.remove(path)

def tsv_collection_maker(header_dict, fas_lib):
    """
    Generates a pairings_tsv.json file in the root_path that contains all pairings.tsv files.
    Parameters
    ----------
    header_dict : dict.keys() == [str], dict.values() == [str]
        dictionary containing all headers indexed by their ENS gene ID
    fas_lib : FAS Library class object
        path to the root of specific species library (e.g. /home/FAS_library/homo_sapiens/release_107/)
    Returns
    -------
    None.
    """
    tsv_dict = dict()
    for gene_id in header_dict.keys():
        tsv_dict[gene_id] = ""
        pairs = itertools.product(header_dict[gene_id], header_dict[gene_id])
        pairs = [pair for pair in pairs if pair[0] <= pair[1]]
        for header_1, header_2 in pairs:
            tsv_dict[gene_id] += header_1 + "\t" + header_2 + "\n"
    with open(fas_lib.get_config("pairings_tsv_json_path"), 'w') as fp:
        json.dump(tsv_dict, fp,  indent=4)

def tsv_access(gene_id, fas_lib):
    """
    Extracts pairing tsv from pairings_tsv.json into root_path/tsv_buffer/gene_id.tsv
    Parameters
    ----------
    gene_id : str
        ENS gene ID
    root_path : TYPE
        path to the root of specific species library (e.g. /home/FAS_library/homo_sapiens/release_107/)
    Returns
    -------
    None.
    """
    with open(fas_lib.get_config("pairings_tsv_json_path"), "r") as fp: 
        tsv_data = json.load(fp)[gene_id]
    with open(fas_lib.get_config("tsv_buffer_path") + gene_id + ".tsv", "w") as fp:
        fp.write(tsv_data)

def tsv_remove(gene_id, fas_lib):
    """
    Removes pairing tsv from pairings_tsv.json into root_path/tsv_buffer/gene_id.tsv
    Parameters
    ----------
    gene_id : TYPE
        ENS gene ID
    root_path : TYPE
        path to the root of specific species library (e.g. /home/FAS_library/homo_sapiens/release_107/)
    Returns
    -------
    None.
    """
    os.remove(fas_lib.get_config("tsv_buffer_path") + gene_id + ".tsv")
        

def parser_setup():
    """
    Reads the parser input.
    Returns
    -------
    Get run options if either to create a tsv file, to delete it again or to concatenate FAS results into a main file.
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library. Is required for library creation.")

    parser.add_argument("-g", "--gene", type=str, default=None,
                        help="""Ensembl gene ID of the gene that the tsv file should be created or deleted for.""")
    
    parser.add_argument("-r", "--remove", action="store_true",
                        help="Remove the tsv file of the given gene.")

    parser.add_argument("-m", "--maketsv", action="store_true",
                        help="Make the tsv file for the given gene.")
    
    parser.add_argument("-j", "--join", action="store_true",
                        help="Join all the FAS_output into the distance_master.phyloprofile.")
    
    args = parser.parse_args()

    config_path = args.config
    gene_id = args.gene

    flag_remove = args.remove
    flag_maketsv = args.maketsv
    flag_join = args.join
    

    return config_path, gene_id, flag_remove, flag_maketsv, flag_join

def main():
    """
    Create tsv output for FAS or delete a tsv that was used by FAS already.
    """
    config_path, gene_id, flag_remove, flag_maketsv, flag_join = parser_setup()
    fas_lib = Library(config_path, False)
    if flag_maketsv:
        tsv_access(gene_id, fas_lib)
    elif flag_remove:
        tsv_remove(gene_id, fas_lib)
    elif flag_join:
        concat_FAS_output(fas_lib)

if __name__ == "__main__":
    main()