#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 12:43:02 2022

@author: chrisbl
"""

import pyranges as pr

from ensembl_access import tsv_to_tuple_list

#generate_expression_summary("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/protein_coding_genes.tsv", "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf")

# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"
# "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/protein_coding_genes.tsv"

def load_expression_gtf(expression_path, flag_filter_unknown=True, flag_remove_transcript_prefix=True):
    gtf = pr.read_gtf(expression_path, as_df=True)
    gtf = [ [gene_id, transcript_id, float(fpkm)] for gene_id,
                                      transcript_id,
                                      fpkm,
                                      feature in zip(gtf['gene_id'],
                                                    gtf['transcript_id'],
                                                    gtf['FPKM'],
                                                    gtf['Feature']) if feature == "transcript"]
    if flag_remove_transcript_prefix:
        gtf = [ [entry[0], entry[1].split(":")[-1], entry[2]] for entry in gtf ]
    if flag_filter_unknown:
        gtf = [ entry for entry in gtf if entry[1].startswith("ENS") ]
    return gtf


def load_protein_coding_ids(path):
    #fas_lib.get_config("protein_coding_ids_path")
    protein_coding_ids = tsv_to_tuple_list(path)
    return protein_coding_ids


def filter_protein_coding(expression_data, protein_coding_ids, fpkm_threshold=0):
    print(len(expression_data))
    # FIlter by protein coding
    protein_coding_transcript_ids = [ transcript_id for gene_id, protein_id, transcript_id in protein_coding_ids ]
    expression_data = [ entry for entry in expression_data if entry[1] in protein_coding_transcript_ids and entry[2] > fpkm_threshold]
    print(len(expression_data))
    return expression_data


def fix_expression_ids(protein_coding_ids, expression_data):
    translate_trans_to_gene_dict = dict()
    for gene_id, protein_id, transcript_id in protein_coding_ids:
        translate_trans_to_gene_dict[transcript_id] = (gene_id, protein_id)
    
    new_expression_data = []
    for gene_id, transcript_id, fpkm in expression_data:
        new_gene_id, new_prot_id = translate_trans_to_gene_dict[transcript_id]
        # If the gene_id is missing.
        if not gene_id.startswith("ENS"):
            new_expression_data.append([new_gene_id, new_prot_id, transcript_id, fpkm])
        # If the gene_id is not missing.
        else:
            # If the gene_id does not belong to the transcript_id
            if new_gene_id != gene_id:
                print(gene_id, "incorrectly associated with", transcript_id)
                print("Dropping from Expression data...")
            # If all is in order
            else:
                new_expression_data.append([gene_id, new_prot_id, transcript_id, fpkm])
    return new_expression_data

              
def extract_total_expression_of_genes(expression_data):
    gene_expression_dict = dict()
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        gene_expression_dict[gene_id] = 0
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        gene_expression_dict[gene_id] += fpkm
    return gene_expression_dict
            

def extract_relative_expression_of_isoform(gene_expression_dict, expression_data):
    isoform_relative_expression_dict = dict()
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        gene_expression_dict[protein_id] = fpkm / gene_expression_dict[gene_id]
    return isoform_relative_expression_dict

def generate_expression_summary(protein_coding_path, expression_path):
    expression_data = load_expression_gtf(expression_path)
    protein_coding_ids = load_protein_coding_ids(protein_coding_path)
    
    expression_data = filter_protein_coding(expression_data, protein_coding_ids)
    
    expression_data = fix_expression_ids(protein_coding_ids, expression_data)
    
    gene_expression_dict = extract_total_expression_of_genes(expression_data)
    isoform_relative_expression_dict = extract_relative_expression_of_isoform(gene_expression_dict,
                                                                              expression_data)
    
    return expression_data, protein_coding_ids, gene_expression_dict, isoform_relative_expression_dict
    
    


