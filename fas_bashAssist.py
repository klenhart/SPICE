#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:31:05 2022

@author: chrisbl
"""

# standard modules
import argparse
import os
import sys

# self-made modules
import library_class
import fas_utility


RAW_SLURM_1 = """#!/bin/bash

#SBATCH --partition={14}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={15}G
#SBATCH --job-name="fas_{4}{12}"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
#SBATCH --array={10}-{11}

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
"""

RAW_SLURM_2 = """
; \
{0} {1} \
-r \
-g $gene \
-c {2}"""


def bash_command_maker(fas_lib, python_path, FAS_handler_path, fas_path, lcr_flag,
tmhmm_flag, partition_list=["all","special","inteli7"], mem_per_cpu="2"):
    """
    Generates SLURM job files for every 1000 genes in the library instance.

    Parameters
    ----------
    fas_lib : Library class object

    python_path : str
        path to the python binary.
    FAS_handler_path : str
        path to fas_handler.py
    fas_path : str
        path to fas.run binary.
    partition_list : list of str
        the partitions to be used for the calculation.
    mem_per_cpu : str
        memory in G to be used for calculation.

    Returns
    -------
    None.

    """
    partitions = ",".join(partition_list)
    with open(fas_lib.get_config("gene_ids_path"), "r") as f:
        gene_ids = f.read().split("\n")
        if gene_ids[0] == "gene":
            gene_ids = gene_ids[1:]
        gene_count = len(gene_ids) - 1
    jobs_ranges = fas_utility.start_stop_range(gene_count, 1000)
    for i, entry in enumerate(jobs_ranges):
        start = 1
        if (entry[1] + 1) % 1000 == 0:
            stop = 1000
        else:
            stop = ((entry[1] + 1) % 1000)
        output_ids = gene_ids[entry[0]:entry[1]+1]
        with open(fas_lib.get_config("slurm_path") + "gene_ids{0}.txt".format(str(i)), "w") as gene_chunk:
            gene_chunk.write("\n".join(output_ids))
        output = RAW_SLURM_1.format(python_path,                      #0
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
                                  fas_lib.get_config("slurm_path"), #13
                                  partitions, #14
                                  mem_per_cpu) #15
        output_2 = RAW_SLURM_2.format(python_path,                      #0
                                      FAS_handler_path,                 #1
                                      fas_lib.get_config("self_path"))  #2
        if lcr_flag or tmhmm_flag:
            if lcr_flag:
                mod_path = fas_lib.get_config("lcr_path")
                job_name = "lcr_FAS_job{0}.job"
            elif tmhmm_flag:
                mod_path = fas_lib.get_config("tmhmm_path")
                job_name = "tmhmm_FAS_job{0}.job"
            output = output + "-d " + mod_path + " \ "
        else:
            job_name = "FAS_job{0}.job"
        output + output_2
        with open(fas_lib.get_config("slurm_path") + job_name.format(str(i)), "w") as f:
            f.write(output)


def parser_setup():
    """
    Reads the parser input.

    Returns
    -------
    Get run options if either to create a tsv file, to delete it again,
    to concatenate FAS results into a main file or to automatically generate bash arrays.
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library.")
    
    parser.add_argument("-m", "--memory", type=str, default="2",
                        help="How much memory to distribute per CPU. Default=2")

    parser.add_argument("-p", "--partitions", nargs="*", action="append", default=["all","special","inteli7"],
                        help="What are the names of the partitions to be used. default=['all', 'special', 'inteli7']")
    
    parser.add_argument("-t", "--tmhmm", action="store_true",
                        help="""Only do the FAS runs using TMHMM and SignalP domains.""")

    parser.add_argument("-l", "--lcr", action="store_true",
                        help="""Only do the FAS runs using flPS and SEG domains.""")

    args = parser.parse_args()

    config_path = args.config
    mem_per_cpu = args.memory
    if len(args.partitions) == 4:
        partition_list = args.partitions[-1]
    else:
        partition_list = args.partitions
    lcr_flag = args.lcr
    tmhmm_flag = args.tmhmm

    return config_path, mem_per_cpu, partition_list, tmhmm_flag, lcr_flag


def main():
    
    # Read out the parser input.
    config_path, mem_per_cpu, partition_list, tmhmm_flag, lcr_flag = parser_setup()
    
    # Construct the fas_handler.py path
    fas_handler_path = os.path.abspath(__file__).split("/")[:-1]
    fas_handler_path.append("fas_handler.py")
    fas_handler_path = "/".join(fas_handler_path)
    if not os.path.isfile(fas_handler_path):
        raise FileNotFoundError("fas_handler.py not found in:", fas_handler_path, ". The file should be in the same directory as fas_bashAssist.py.")
        sys.exit(1)
        
    # Get python path
    python_path = sys.executable
    
    # Get the fas.run binary path from python binary path.
    fas_path= python_path.split("/")[:-1]
    fas_path.append("fas.run")
    fas_path = "/".join(fas_path)
    
    # Load fas_lib from config.tsv
    fas_lib = library_class.Library(config_path, False)
    
    # Assemble the bash arrays
    bash_command_maker(fas_lib,
                       python_path,
                       fas_handler_path,
                       fas_path,
                       lcr_flag,
                       tmhmm_flag,
                       partition_list,
                       mem_per_cpu)


if __name__ == "__main__":
    main()
