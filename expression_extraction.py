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


def filter_protein_coding(expression_data, protein_coding_ids):
    # FIlter by protein coding
    protein_coding_transcript_ids = [ transcript_id for gene_id, protein_id, transcript_id in protein_coding_ids ]
    expression_data = [ [gene_id, transcript_id, fpkm] for gene_id, transcript_id, fpkm in expression_data if transcript_id in protein_coding_transcript_ids]
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
        isoform_relative_expression_dict[protein_id] = fpkm / gene_expression_dict[gene_id]
    return isoform_relative_expression_dict

def filter_non_expressed(expression_data, fpkm_threshold=0):
    return [ [gene_id, transcript_id, fpkm] for gene_id, transcript_id, fpkm in expression_data if fpkm > fpkm_threshold ]

def make_isoform_protein_id_dict(expression_data):
    isoform_protein_id_dict = dict()
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        isoform_protein_id_dict[gene_id] = []
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        isoform_protein_id_dict[gene_id].append(protein_id)
    return isoform_protein_id_dict


def isoform_count_distribution_check(isoform_protein_id_dict):
    isoform_count_list= []
    three_or_more_count = 0
    for isoform_list in isoform_protein_id_dict.values():
        isoform_count_list.append(len(isoform_list))
        if len(isoform_list) >= 3:
            three_or_more_count += 1
    isoform_count_list.sort()
    three_or_more_ratio = three_or_more_count / len(isoform_protein_id_dict.values())
    median = isoform_count_list[int(len(isoform_count_list)/2)]
    mean = sum(isoform_count_list) / len(isoform_count_list)
    return median, mean, three_or_more_ratio

def generate_expression_summary(protein_coding_path, expression_path):
    expression_data = load_expression_gtf(expression_path)
    protein_coding_ids = load_protein_coding_ids(protein_coding_path)
    
    expression_data = filter_non_expressed(expression_data)
    
    expression_data = filter_protein_coding(expression_data, protein_coding_ids)
    
    expression_data = fix_expression_ids(protein_coding_ids, expression_data)
    
    gene_expression_dict = extract_total_expression_of_genes(expression_data)
    isoform_relative_expression_dict = extract_relative_expression_of_isoform(gene_expression_dict,
                                                                              expression_data)
    
    isoform_protein_id_dict = make_isoform_protein_id_dict(expression_data)
    
    # Check how many isoforms are expressed per gene.
    return  isoform_protein_id_dict #expression_data, gene_expression_dict, isoform_relative_expression_dict

def main():
    expression_paths = ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856511.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856512.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856513.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856514.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856515.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856516.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856517.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856518.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856519.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856520.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578910.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
"/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578911.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"]
    protein_coding_path = '/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/protein_coding_genes.tsv'
    stats = [generate_expression_summary(protein_coding_path, expression_path) for expression_path in expression_paths]
    
    median_list = []
    mean_list = []
    three_or_more_list = []
    for median, mean, three_or_more_ratio in stats:
        median_list.append(median)
        mean_list.append(mean)
        three_or_more_list.append(three_or_more_ratio)
    print("Median isoform count:", sum(median_list) / len(median_list))
    print("Mean isoform count:", sum(mean_list) / len(mean_list))
    print("three or more isoform count:", sum(three_or_more_list) / len(three_or_more_list))
        
        

if __name__ == "__main__":
    main()


