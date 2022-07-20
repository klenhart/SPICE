#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:23:30 2022

@author: chrisbl
"""

def extract_transcripts_from_fasta_header(path):
    with open(path, "r") as fasta:
        lines = fasta.readlines()[::2]
    lines = [line.split("|")[0:2] for line in lines]
    int_ids = [(int(gene_id[5:]), int(transcript_id[4:])) for gene_id, transcript_id in lines]
    return int_ids

def check_for_transcript(gene_id, transcript_id, isoforms_path):
    gene_num = int(gene_id[4:])
    transcript_num = int(transcript_id[4:])
    int_ids = extract_transcripts_from_fasta_header(isoforms_path)
    flag = (gene_num, transcript_num) in int_ids
    return flag

def main():
    test_path = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/isoform.fasta"
    gene_id = "ENSG00000092068"
    transcript_id = "ENST00000529705"
    #Should be False
    print(check_for_transcript(gene_id, transcript_id, test_path))
    gene_id = "ENSG00000192068"
    transcript_id = "ENST00000525062"
    #Should be True
    print(check_for_transcript(gene_id, transcript_id, test_path))
    gene_id = "ENSG00000392068"
    transcript_id = "ENST00000525062"
    #Should be True
    print(check_for_transcript(gene_id, transcript_id, test_path))

if __name__ == "__main__":
    main()