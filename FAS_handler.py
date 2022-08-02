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

# python /home/chrisbl/project/FAS_Pipe/Scripts/grand-trumpet/FAS_handler.py -b -o /share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107 -p /home/chrisbl/miniconda3/envs/FAS/bin/python -s /home/chrisbl/project/FAS_Pipe/Scripts/grand-trumpet/FAS_handler.py -f /home/chrisbl/miniconda3/envs/FAS/bin/fas.run

TEST = {"A" : ["A", "B", "C", "D"], "B" : ["1", "2", "3", "4"]}

RAW_SLURM = """!/bin/bash

SBATCH --partition=all,special,inteli7
SBATCH --cpus-per-task=1
SBATCH --mem-per-cpu=2G
SBATCH --job-name="fas_{4}"
SBATCH --output=/dev/null 
SBATCH --error=/dev/null
SBATCH --array={5}-{6}

gene=$(awk FNR==$SLURM_ARRAY_TASK_ID "{2}/gene_ids{7}.txt")
{0} {1} \
-m \
-g $gene \
-o {2} \
&& \
{3} \
--seed {2}/isoforms.fasta \
--query {2}/isoforms.fasta \
--annotation_dir {2} \
--out_dir {2}/FAS_buffer/ \
--bidirectional \
--pairwise {2}/tsv_buffer/$gene.tsv \
--out_name $gene \
--tsv \
--phyloprofile {2}/phyloprofile_ids.tsv \
--domain \
--empty_as_1 \
; \
{0} {1} \
-r \
-g $gene \
-o {2}"""

# {0} {1} -b -o {2} -p {0} -s {1} -f {3}

def start_stop_range(length, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(1, length, n):
        yield (i, min(i+n-1, length))

def bash_command_maker(root_path, python_path, FAS_handler_path, fas_path):
    slurm_path = root_path + "/SLURM/"
    if not os.path.exists(slurm_path):
        os.makedirs(slurm_path)
    with open(root_path + "/gene_ids.txt", "r") as f:
        gene_ids = f.read().split("\n")
        gene_count = len(gene_ids) - 1
    jobs_ranges = start_stop_range(gene_count, 1000)
    species = root_path.split("/")
    index = species.index("FAS_library") + 1
    species = species[index]
    for i, entry in enumerate(jobs_ranges):
        start = 1
        if entry[1] % 1000 == 0:
            stop = 1000
        else:
            stop = (entry[1] % 1000) - 1
        output_ids = ["gene"] + gene_ids[entry[0]:entry[1]+1]
        with open(slurm_path + "gene_ids{0}.txt".format(str(i)), "w") as gene_chunk:
            gene_chunk.write("\n".join(output_ids))
        output = RAW_SLURM.format(python_path, 
                                  FAS_handler_path, 
                                  root_path, 
                                  fas_path, 
                                  species, 
                                  str(start), 
                                  str(stop),
                                  str(i))
        with open(slurm_path + "FAS_job{0}.job".format(str(i)), "w") as f:
            f.write(output)

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
    access_command = " ". join(options_access)
    remove_command = " ".join(options_remove)
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
    fas_command = " ".join(options)
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

    root_path = args.output
    gene_id = args.gene
    python_path = args.python
    FAS_handler_path = args.handler
    fas_path = args.fas
    flag_remove = args.remove
    flag_maketsv = args.maketsv
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