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


def fix_expression_ids( expression_data, protein_coding_ids):
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

def check_for_same_order(m1, m2):
    flag_check = True
    for i, entry in enumerate(m1):
        gene_id_1, transcript_id_1, fpkm_1 = entry
        gene_id_2, transcript_id_2, fpkm_2 = m2[i]
        if gene_id_1 != gene_id_2 or transcript_id_1 != transcript_id_2:
            flag_check = False
            print(i)
            print(gene_id_1, " x ", gene_id_2)
            print(transcript_id_1, " x ", transcript_id_2)
            break
    return flag_check

def join_expression(expression_path_list, protein_coding_path):
    protein_coding_ids = load_protein_coding_ids(protein_coding_path)
    expression_dict = dict()
    prot_to_gene_dict = dict()
    for expression_path in expression_path_list:
        expression_data = load_expression_gtf(expression_path)
        expression_data = filter_non_expressed(expression_data)
        expression_data = filter_protein_coding(expression_data, protein_coding_ids)
        expression_data = fix_expression_ids(expression_data, protein_coding_ids)
        for gene_id, prot_id, transcript_id, fpkm in expression_data:
            if prot_id not in expression_dict.keys():
                expression_dict[prot_id] = fpkm
            else:
                expression_dict[prot_id] += fpkm
            prot_to_gene_dict[prot_id] = gene_id
    isoforms_dict = dict()
    for prot_id in expression_dict.keys():
        gene_id = prot_to_gene_dict[prot_id]
        if gene_id in isoforms_dict.keys():
            isoforms_dict[gene_id].append(prot_id)
        else:
            isoforms_dict[gene_id] = [prot_id]
    return expression_dict, isoforms_dict


def load_dist_matrix(fas_lib):
    distance_master_path = fas_lib.get_config["distance_master_path"]
    with open(distance_master_path, "r") as f:
        distance_master = f.read()
    distance_master = distance_master.split("\n")
    distance_master = [ entry.split("\t") for entry in distance_master ]
    dist_matrix_dict = dict()
    for seed, query, fas_1, fas_2 in distance_master:
        fas_score = (float(fas_1) + float(fas_2)) / 2
        gene_id, prot_id_1, tax_id = seed.split("|")
        gene_id, prot_id_2, tax_id = query.split("|")
        if gene_id not in dist_matrix_dict.keys():
            dist_matrix_dict[gene_id] = {prot_id_1 : { prot_id_2 : fas_score} , prot_id_2 : { prot_id_1 : fas_score} }
        else:
            if prot_id_1 not in dist_matrix_dict[gene_id].keys():
                dist_matrix_dict[gene_id][prot_id_1] = { prot_id_2 : fas_score}
            else:
                dist_matrix_dict[gene_id][prot_id_1][prot_id_2] = fas_score
            if prot_id_2 not in dist_matrix_dict[gene_id].keys():
                dist_matrix_dict[gene_id][prot_id_2] = { prot_id_1 : fas_score}
            else:
                dist_matrix_dict[gene_id][prot_id_2][prot_id_1] = fas_score
    return dist_matrix_dict


def make_fas_graph_row(expression_dict, isoforms_dict, dist_matrix_dict, gene_id):
    # Start the row
    fas_graph_row = gene_id + "\t"
    
    # Extract the IDs of all isoforms of the gene in a list
    isoform_prot_ids = isoforms_dict[gene_id]
    
    # Extract the distance matrix containing th FAS scores
    dist_matrix = dist_matrix_dict[gene_id]
    
    # Initialize the total fpkm of the gene
    total_fpkm = 0
    
    # initialize dictionary to store only the fpkm of the protein_ids
    relative_expression_dict = dict()
    fas_sigma_dict = dict()
    
    # Sum up the total fpkm and store all single fpkm values for each protein id
    for prot_id in isoform_prot_ids:
        total_fpkm += expression_dict[prot_id]
        relative_expression_dict[prot_id] = expression_dict[prot_id]
        fas_sigma_dict[prot_id] = 0
    
    # Turn the fpkm values into relative values.
    for prot_id in isoform_prot_ids:
        relative_expression_dict[prot_id] = relative_expression_dict[prot_id] / total_fpkm

    # Calculate the Sigma values.
    for seed_id in isoform_prot_ids:
        for query_id in isoform_prot_ids:
            fas_sigma_dict[seed_id] += (dist_matrix[seed_id][query_id] * relative_expression_dict[query_id])
    
    for prot_id in isoform_prot_ids:
        fas_graph_row += prot_id + ":" + str(fas_sigma_dict[prot_id]) + ";"
    fas_graph_row = fas_graph_row[:-1]
    
    fas_graph_row = fas_graph_row + "\t" + str(total_fpkm)

    return fas_graph_row
    

def generate_FAS_polygon(fas_lib, expression_paths_path, name_path):
    dist_matrix_dict = load_dist_matrix(fas_lib)
    with open(expression_paths_path, "r") as f:
        expression_paths_list = f.read()
    expression_paths_list = expression_paths_list.split("\n")
    expression_paths_list = [ entry.split(";") for entry in expression_paths_list ]
    
    with open(name_path, "r") as f:
        name_list = f.read()
    name_list = name_list.split("\n")

    for i, expression_paths in enumerate(expression_paths_list):
        name = name_list[i]
        expression_dict, isoforms_dict = join_expression(expression_paths,
                                                         fas_lib.get_config("protein_coding_path"))
        polygon_output = "geneID\tisoformVectors\tgeneFPKM"
        for gene_id in isoforms_dict.keys():
            polygon_output += "\n" + make_fas_graph_row(expression_dict,
                                                        isoforms_dict,
                                                        dist_matrix_dict,
                                                        gene_id)
        with open(fas_lib.get_config("root_path") + "polygonFAS_{0}.tsv".format(name), "w") as f:
            f.write(polygon_output)
    


def main():
    expression_path_list = ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
        "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856511.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
        "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856512.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"]
    protein_coding_path = '/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/protein_coding_genes.tsv'
    expression_dict, isoforms_dict = join_expression(expression_path_list, protein_coding_path)
#     expression_paths = ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856511.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856512.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856513.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856514.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856515.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856516.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856517.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856518.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856519.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856520.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578910.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
# "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578911.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"]
#     protein_coding_path = '/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/protein_coding_genes.tsv'
#     stats = [generate_expression_summary(expression_path, protein_coding_path) for expression_path in expression_paths]
    
#     for stat_dict in stats:
#         for gene in stat_dict.keys():
            
    
#     median_list = []
#     mean_list = []
#     three_or_more_list = []
#     for median, mean, three_or_more_ratio in stats:
#         median_list.append(median)
#         mean_list.append(mean)
#         three_or_more_list.append(three_or_more_ratio)
#     print("Median isoform count:", sum(median_list) / len(median_list))
#     print("Mean isoform count:", sum(mean_list) / len(mean_list))
#     print("three or more isoform count:", sum(three_or_more_list) / len(three_or_more_list))
        
        

if __name__ == "__main__":
    main()


