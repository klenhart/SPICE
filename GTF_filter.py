#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:48:44 2022

@author: chrisbl
"""

import argparse
import os
import subprocess
import pyranges as pr

#/share/project/zarnack/chrisbl/FAS/utility/protein_lib/library/transcript_ids.txt

# python GTF_filter.py -t /home/chrisbl/project/FAS_Pipe/transcript_ids_test.txt -g /home/chrisbl/project/FAS_Pipe/gtf_paths_test.txt
def parser_setup():
    """
    

    Returns
    -------
    Path to transcript IDs
    
    and
    
    Path to paths of gtf files.

    """
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-t", "--transcripts", type=str, 
                        help="path of text file containng linebreak delimited transcript IDs.")
    
    parser.add_argument("-g", "--gtfs", type=str,
                        help="""path of text file containing linebreak delimited GTF file paths.
                        The filtered path will be put output in the same path as the original GTF files.""")
    
    args = parser.parse_args()
    
    gtf_paths_path = args.gtfs
    transcripts_path = args.transcripts
    
    return gtf_paths_path, transcripts_path

def extract_transcripts(path):
    """
    

    Parameters
    ----------
    path : path to a file that contains newline delimited transcript IDs

    Returns
    -------
    transcript_list : list of all transcripts found in the given file.

    """
    transcript_list = list()
    with open(path, "r") as f:
        for line in f:
            transcript_list.append(line.strip())
    return transcript_list


def extract_gtf_paths(path):
    """
    

    Parameters
    ----------
    path : path to a text file that contains newline delimited paths to gtf files

    Returns
    -------
    gtf_filenames : filenames of the gtf files
    gtf_directories : path to the directories of each gtf file.

    """
    gtf_paths = list()
    with open(path, "r") as f: 
        for line in f:
            gtf_paths.append(line.strip())
    gtf_paths = [path.split("/") for path in gtf_paths]
    gtf_filenames = [path[-1] for path in gtf_paths]
    gtf_directories = ["/".join(path[:-1]) + "/" for path in gtf_paths]
    return gtf_filenames, gtf_directories


def filter_gtf_files(gtf_directory_paths, gtf_filenames, transcript_list):
    stable_filter_id_vector = [transcript_id.split(".")[0] for transcript_id in transcript_list]
    print("Check1")
    print(stable_filter_id_vector[0:10])
    print("Check2")
    for i, path in enumerate(gtf_directory_paths):
        gtf = pr.read_gtf(path + gtf_filenames[i])
        gtf_df = gtf.df
        id_column = [transcript_id.split(".")[0] for transcript_id in gtf_df["transcript_id"]]
        inclusion_vector = [transcript_id in stable_filter_id_vector for transcript_id in id_column]
        print(any(inclusion_vector))
        break
        

def main():
    gtf_paths_path, transcripts_path = parser_setup()
    
    transcript_list = extract_transcripts(transcripts_path)

    gtf_filenames, gtf_directory_paths = extract_gtf_paths(gtf_paths_path)

    filter_gtf_files(gtf_directory_paths, gtf_filenames, transcript_list)
    


if __name__ == "__main__":
    main()