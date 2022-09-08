#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 12:43:02 2022

@author: chrisbl
"""

import pyranges as pr
import os
import json

import ensembl_access


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
    protein_coding_ids = ensembl_access.tsv_to_tuple_list(path)
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


def filter_exempt_genes(expression_data, exempt_genes):
    if len(exempt_genes) == 0:
        pass
    else:
        expression_data = [ entry for entry in expression_data if entry[0] not in exempt_genes ]
    return expression_data

def join_expression(expression_path_list, protein_coding_path, exempt_genes):
    protein_coding_ids = load_protein_coding_ids(protein_coding_path)
    protein_coding_ids = [ entry for entry in protein_coding_ids if entry[0] not in exempt_genes ]
    expression_dict = dict()
    prot_to_gene_dict = dict()
    for gene_id, prot_id, transcript_id in protein_coding_ids:
        expression_dict[prot_id] = 0
        prot_to_gene_dict[prot_id] = gene_id
    # Extract all the expression from the gtf files and filter out everything unnecessary.
    for expression_path in expression_path_list:
        expression_data = load_expression_gtf(expression_path)
        expression_data = filter_non_expressed(expression_data)
        expression_data = filter_protein_coding(expression_data, protein_coding_ids)
        expression_data = fix_expression_ids(expression_data, protein_coding_ids)
        expression_data = filter_exempt_genes(expression_data, exempt_genes)
        for gene_id, prot_id, transcript_id, fpkm in expression_data:
            expression_dict[prot_id] += fpkm
    isoforms_dict = dict()
    for prot_id in expression_dict.keys():
        gene_id = prot_to_gene_dict[prot_id]
        if gene_id in isoforms_dict.keys():
            isoforms_dict[gene_id].append(prot_id)
        else:
            isoforms_dict[gene_id] = [prot_id]
    return expression_dict, isoforms_dict

def load_dist_matrix(fas_lib):
    if "distance_master.json" in os.listdir(fas_lib.get_config("root_path")):
        with open(fas_lib.get_config("root_path") + "distance_master.json", "r") as f: 
            dist_matrix_dict = json.load(f)
    else:
        distance_master_path = fas_lib.get_config("distance_master_path")
        with open(distance_master_path, "r") as f:
            distance_master = f.read()
        distance_master = distance_master.split("\n")   
        distance_master = distance_master[1:]
        if distance_master[-1] == "":
            distance_master = distance_master[:-1]
        distance_master = [ entry.split("\t") for entry in distance_master ]
        dist_matrix_dict = dict()
        for seed, tax_id, query, fas_1, fas_2 in distance_master:
            #fas_score = (float(fas_1) + float(fas_2)) / 2
            gene_id, prot_id_1, tax_id = seed.split("|")
            gene_id, prot_id_2, tax_id = query.split("|")
            if gene_id not in dist_matrix_dict.keys():
                dist_matrix_dict[gene_id] = {prot_id_1 : { prot_id_2 : fas_1} , prot_id_2 : { prot_id_1 : fas_2} }
            else:
                if prot_id_1 not in dist_matrix_dict[gene_id].keys():
                    dist_matrix_dict[gene_id][prot_id_1] = { prot_id_2 : fas_1}
                else:
                    dist_matrix_dict[gene_id][prot_id_1][prot_id_2] = fas_1
                if prot_id_2 not in dist_matrix_dict[gene_id].keys():
                    dist_matrix_dict[gene_id][prot_id_2] = { prot_id_1 : fas_2}
                else:
                    dist_matrix_dict[gene_id][prot_id_2][prot_id_1] = fas_2
        with open(fas_lib.get_config("root_path") + "distance_master.json", 'w') as f:
            json.dump(dist_matrix_dict, f,  indent=4)
    return dist_matrix_dict


def make_fas_graph_row(expression_dict, isoforms_dict, dist_matrix_dict, gene_id):
    # Start the row
    fas_graph_row = gene_id + "\t"
    
    # Extract the IDs of all isoforms of the gene in a list
    isoform_prot_ids = isoforms_dict[gene_id]
    
    # Extract the distance matrix containing th FAS scores
    dist_matrix = dist_matrix_dict[gene_id]
    
    # Initialize the total fpkm of the gene
    total_fpkm = 0.0
    
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
        if total_fpkm == 0.0:
            relative_expression_dict[prot_id] = 0.0
        else:
            relative_expression_dict[prot_id] = relative_expression_dict[prot_id] / total_fpkm

    # Calculate the Sigma values.
    for seed_id in isoform_prot_ids:
        for query_id in isoform_prot_ids:
            fas_sigma_dict[seed_id] += (dist_matrix[seed_id][query_id] * relative_expression_dict[query_id])
    
    for prot_id in isoform_prot_ids:
        fas_graph_row += prot_id + ":" + str(fas_sigma_dict[prot_id]) + ":" + str(relative_expression_dict[prot_id]) + ";"
    fas_graph_row = fas_graph_row[:-1]
    
    fas_graph_row = fas_graph_row + "\t" + str(total_fpkm)
    return fas_graph_row
    

def generate_FAS_polygon(fas_lib, expression_paths_path, name_path):
    dist_matrix_dict = load_dist_matrix(fas_lib)
    with open(expression_paths_path, "r") as f:
        expression_paths_list = f.read()
    expression_paths_list = expression_paths_list.split("\n")
    expression_paths_list = [ expression_paths for expression_paths in expression_paths_list if len(expression_paths) > 0 ]
    expression_paths_list = [ entry.split(";") for entry in expression_paths_list ]
    
    with open(name_path, "r") as f:
        name_list = f.read()
    name_list = name_list.split("\n")
    for i, expression_paths in enumerate(expression_paths_list):
        name = name_list[i]
        expression_dict, isoforms_dict = join_expression(expression_paths,
                                                         fas_lib.get_config("protein_coding_ids_path"),
                                                         ["ENSG00000155657"])
        polygon_output = "geneID\tisoformVectors\tgeneFPKM"
        for gene_id in isoforms_dict.keys():
            polygon_output += "\n" + make_fas_graph_row(expression_dict,
                                                        isoforms_dict,
                                                        dist_matrix_dict,
                                                        gene_id)
        if not os.path.exists(fas_lib.get_config("root_path") + "FAS_polygon"):
            os.makedirs(fas_lib.get_config("root_path") + "FAS_polygon")
        with open(fas_lib.get_config("root_path") +  "FAS_polygon/polygonFAS_{0}.tsv".format(name), "w") as f:
            f.write(polygon_output)

def main():
    pass

if __name__ == "__main__":
    main()
