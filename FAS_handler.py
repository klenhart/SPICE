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

# python /home/chrisbl/project/FAS_Pipe/Scripts/grand-trumpet/FAS_handler.py -b -c /share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/config.tsv -p /home/chrisbl/miniconda3/envs/FAS/bin/python -s /home/chrisbl/project/FAS_Pipe/Scripts/grand-trumpet/FAS_handler.py -f /home/chrisbl/miniconda3/envs/FAS/bin/fas.run

TEST = {"A" : ["A", "B", "C", "D"], "B" : ["1", "2", "3", "4"]}

RAW_SLURM = """!/bin/bash

SBATCH --partition=all,special,inteli7
SBATCH --cpus-per-task=1
SBATCH --mem-per-cpu=2G
SBATCH --job-name="fas_{4}{12}"
SBATCH --output=/dev/null 
SBATCH --error=/dev/null
SBATCH --array={10}-{11}

gene=$(awk FNR==$SLURM_ARRAY_TASK_ID "{13}gene_ids{12}.txt")
{0} {1} \
-m \
-g $gene \
-c {2} \
&& \
{3} \
--seed {5} \
--query {5} \
--annotation_dir {6} \
--out_dir {7} \
--bidirectional \
--pairwise {8}$gene.tsv \
--out_name $gene \
--tsv \
--phyloprofile {9} \
--domain \
--empty_as_1 \
; \
{0} {1} \
-r \
-g $gene \
-c {2}"""

def start_stop_range(length, n):
    """
    Splits length into equal chunks of size n.
    
    Parameters
    ----------
    length : int
        length to be split apart.
    n : int
        size of chunks.

    Yields
    ------
    generator of (int, int)
        Splits length into equal chunks of size n.
    """
    for i in range(1, length, n):
        yield (i, min(i+n-1, length))

def bash_command_maker(fas_lib, python_path, FAS_handler_path, fas_path):
    """
    Generates SLURM job files for every 1000 genes in the library instance.

    Parameters
    ----------
    root_path : str
        path to the root of specific species library (e.g. /home/FAS_library/homo_sapiens/release_107/)
    python_path : str
        path to the python binary.
    FAS_handler_path : str
        path to FAS_handler.py
    fas_path : str
        path to fas.run binary.

    Returns
    -------
    None.

    """
    with open(fas_lib.get_config("gene_ids_path"), "r") as f:
        gene_ids = f.read().split("\n")
        gene_count = len(gene_ids) - 1
    jobs_ranges = start_stop_range(gene_count, 1000)
    for i, entry in enumerate(jobs_ranges):
        start = 1
        if entry[1] % 1000 == 0:
            stop = 1000
        else:
            stop = (entry[1] % 1000)
        output_ids = ["gene"] + gene_ids[entry[0]:entry[1]+1]
        with open(fas_lib.get_config("slurm_path") + "gene_ids{0}.txt".format(str(i)), "w") as gene_chunk:
            gene_chunk.write("\n".join(output_ids))
        output = RAW_SLURM.format(python_path,                      #0
                                  FAS_handler_path,                 #1
                                  fas_lib.get_config("self_path"),  #2
                                  fas_path,                         #3
                                  fas_lib.get_config("species"),    #4
                                  fas_lib.get_config("isoforms_path"), #5
                                  fas_lib.get_config("annotation_path"), #6
                                  fas_lib.get_config("fas_buffer_path"), #7
                                  fas_lib.get_config("tsv_buffer_path"), #8
                                  fas_lib.get_config("phyloprofile_ids_path"), #9
                                  str(start), #10
                                  str(stop), #11
                                  str(i),   #12
                                  fas_lib.get_config("slurm_path")) #13
        with open(fas_lib.get_config("slurm_path") + "FAS_job{0}.job".format(str(i)), "w") as f:
            f.write(output)

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
        pairs = [pair for pair in pairs if pair[0] < pair[1]]
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
    
    parser.add_argument("-p", "--python", type=str, default=None,
                        help="""Specify the location of the pythonversion to use. This is only necessary when using the --bash option.""")

    parser.add_argument("-s", "--handler", type=str, default=None,
                        help="""Specify the location FAS_handler.py. This is only necessary when using the --bash option.""")
    
    parser.add_argument("-f", "--fas", type=str, default=None,
                        help="""Specify the location of fas.run. This is only necessary when using the --bash option.""")
    
    parser.add_argument("-g", "--gene", type=str, default=None,
                        help="""Ensembl gene ID of the gene that the tsv file should be created or deleted for.""")
    
    parser.add_argument("-r", "--remove", action="store_true",
                        help="Remove the tsv file of the given gene.")

    parser.add_argument("-m", "--maketsv", action="store_true",
                        help="Make the tsv file for the given gene.")
    
    parser.add_argument("-b", "--bash", action="store_true",
                        help="Make one FAS bash command for each gene in a text file.")
    

    args = parser.parse_args()

    config_path = args.config
    gene_id = args.gene
    python_path = args.python
    FAS_handler_path = args.handler
    fas_path = args.fas
    flag_remove = args.remove
    flag_maketsv = args.maketsv
    flag_bash = args.bash

    return config_path, python_path, FAS_handler_path, fas_path, gene_id, flag_remove, flag_maketsv, flag_bash

def main():
    """
    Create tsv output for FAS or delete a tsv that was used by FAS already.
    """
    config_path, python_path, FAS_handler_path, fas_path, gene_id, flag_remove, flag_maketsv, flag_bash = parser_setup()
    fas_lib = Library(config_path, False)
    if flag_maketsv:
        tsv_access(gene_id, fas_lib)
    elif flag_remove:
        tsv_remove(gene_id, fas_lib)
    elif flag_bash:
        bash_command_maker(fas_lib, python_path, FAS_handler_path, fas_path)
        

if __name__ == "__main__":
    main()