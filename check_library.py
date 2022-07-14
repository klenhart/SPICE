#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:23:30 2022

@author: chrisbl
"""

def extract_transcripts_from_fasta_header(path):
    with open(path, "r") as fasta:
        lines = fasta.readlines()[::2]
    lines = [line.split(" ")[1] for line in lines]
    return lines

def check_for_transcript(gene_id, transcript_id, library_path):
    path = library_path + "/" + gene_id + "/isoform.fasta"
    flag = transcript_id in extract_transcripts_from_fasta_header(path)
    return flag

def main():
    #Should be False
    print(check_for_transcript("ENSG00000000419", "ENST00000371488", "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/library_test/"))
    #Should be True
    print(check_for_transcript("ENSG00000000419", "ENST00000371588", "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/library_test/"))

if __name__ == "__main__":
    main()