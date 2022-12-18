#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################################################
# Copyright (C) 2022 Christian, Blümel, Julian Dosch
#
# This file is part of grand-trumpet.
#
#  grand-trumpet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  grand-trumpet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with grand-trumpet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
"""
Created on Wed Jul 27 10:00:07 2022
@author: chrisbl
"""

import json
import argparse
import os

from valves.library_class import Library

def concat_FAS_output(fas_lib, flag_lcr, flag_tmhmm):
    if flag_lcr:
        distance_master_path = fas_lib.get_config("fas_lcr_path")
        forward_domain_master_path = fas_lib.get_config("forward_domain_lcr_path")
        reverse_domain_master_path = fas_lib.get_config("reverse_domain_lcr_path")
        fas_buffer_path = fas_lib.get_config("fas_buffer_lcr_path")
    elif flag_tmhmm:
        distance_master_path = fas_lib.get_config("fas_tmhmm_path")
        forward_domain_master_path = fas_lib.get_config("forward_domain_tmhmm_path")
        reverse_domain_master_path = fas_lib.get_config("reverse_domain_tmhmm_path")
        fas_buffer_path = fas_lib.get_config("fas_buffer_tmhmm_path")
    else:
        distance_master_path = fas_lib.get_config("fas_all_path")
        forward_domain_master_path = fas_lib.get_config("forward_domain_all_path")
        reverse_domain_master_path = fas_lib.get_config("reverse_domain_all_path")
        fas_buffer_path = fas_lib.get_config("fas_buffer_path")

    tsv_path = fas_lib.get("tsv_buffer_path")
    still_running_names = [ name[:-4] for name in os.listdir(tsv_path) ]
    file_names = [name for name in os.listdir(fas_buffer_path) if not (name in still_running_names)]
    file_names_forward_domain = [ name for name in file_names if name.endswith("forward.domains") ]
    file_names_reverse_domain = [ name for name in file_names if name.endswith("reverse.domains") ]
    file_names_fas = [ name for name in file_names if name.endswith(".phyloprofile") ]
    print("Forward Domain files: ", len(file_names_forward_domain))
    print("Reverse Domain files: ", len(file_names_reverse_domain))
    print("FAS files: ", len(file_names_fas))
    for filename in file_names_fas:
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
    for filename in file_names_forward_domain:
        path = fas_buffer_path + filename
        with open(path, "r") as f_in:
            query = f_in.read().split("\n")
            if len(query[-1]) == 0:
                query = query[0:-1]
            query = "\n".join(query)
            with open(forward_domain_master_path, "a") as f_out:
                f_out.write(query + "\n")
        os.remove(path)
    for filename in file_names_reverse_domain:
        path = fas_buffer_path + filename
        with open(path, "r") as f_in:
            query = f_in.read().split("\n")
            if len(query[-1]) == 0:
                query = query[0:-1]
            query = "\n".join(query)
            with open(reverse_domain_master_path, "a") as f_out:
                f_out.write(query + "\n")
        os.remove(path)    

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
    
    parser.add_argument("-t", "--tmhmm", action="store_true",
                        help="""Only do the FAS runs using TMHMM and SignalP domains.""")

    parser.add_argument("-l", "--lcr", action="store_true",
                        help="""Only do the FAS runs using flPS and SEG domains.""")
    
    args = parser.parse_args()

    config_path = args.config
    gene_id = args.gene

    flag_remove = args.remove
    flag_maketsv = args.maketsv
    flag_join = args.join
    flag_lcr = args.lcr
    flag_tmhmm = args.tmhmm

    return config_path, gene_id, flag_remove, flag_maketsv, flag_join, flag_lcr, flag_tmhmm

def main():
    """
    Create tsv output for FAS or delete a tsv that was used by FAS already.
    """
    config_path, gene_id, flag_remove, flag_maketsv, flag_join, flag_lcr, flag_tmhmm = parser_setup()
    fas_lib = Library(config_path, False)
    if flag_maketsv:
        tsv_access(gene_id, fas_lib)
    elif flag_remove:
        tsv_remove(gene_id, fas_lib)
    elif flag_join:
        concat_FAS_output(fas_lib, flag_lcr, flag_tmhmm)

if __name__ == "__main__":
    main()