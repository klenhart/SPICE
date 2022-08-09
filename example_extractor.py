#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:30:12 2022

@author: chrisbl
"""

from expression_extraction import join_expression
from expression_extraction import load_dist_matrix
from expression_extraction import make_fas_graph_row

from library_class import Library

import json


# Long read data Link: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7334/samples/

# P633L A2 : ERR4578910
# P633L A2B : ERR4578911

# R634Q A2 : ERR2856510 ERR2856511
# R634Q A2B : ERR2856512 ERR2856513

# WTNC A2 : ERR2856514
# WTNC A2B : ERR2856515

# WT A2 : ERR2856516 ERR2856517
# WT    ERR2856518 ERR2856519 ERR2856520

# bis ENSG00000135000 habe ich die FAS Scores, also kann ich RECK, ARNTL, DPF2 vergleichen.

flag = True

test_ids_dict = { "RECK" : "ENSG00000122707",
                 "ARNTL" : "ENSG00000133794",
                 "DPF2" : "ENSG00000133884",
                 "PITX2" : "ENSG00000164093",
                 "TPM4" : "ENSG00000167460",
                 "SETD5" : "ENSG00000168137",
                 "UPP1" : "ENSG00000183696",
                 "MROH7" : "ENSG00000184313",
                 "BCO2" : "ENSG00000197580",     
                 "TPM2" : "ENSG00000198467" }

fas_lib = Library("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/config.tsv", False)

dist_matrix = load_dist_matrix(fas_lib)

# dist_matrix_DPF2 = dist_matrix["ENSG00000133884"]
# dist_matrix_ARNTL = dist_matrix["ENSG00000133794"]

#make_fas_graph_row(expression_dict, isoforms_dict, dist_matrix_dict, gene_id)

expression_paths = {"R634Q_A2" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856510.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
              "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856511.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"R634Q_A2B" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856512.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
               "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856513.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"WTNC_A2" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856514.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"WTNC_A2B" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856515.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"WT_A2" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856516.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
           "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856517.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"WT_A2B" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856518.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
            "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856519.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf",
            "/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR2856520.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"P633L_A2" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578910.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"],
"P633L_A2B" : ["/share/gluster/Projects/FeatureArchitectureUniverse/gtf/ERR4578911.fastq.gz_desalt.sort.bam.out_stringtie_recount.gtf"] }

if flag:
    expression_dict_R634Q_A2, isoforms_dict_R634Q_A2 = join_expression(expression_paths["R634Q_A2"], fas_lib.get_config("protein_coding_ids_path"))
    expression_dict_R634Q_A2B, isoforms_dict_R634Q_A2B = join_expression(expression_paths["R634Q_A2B"], fas_lib.get_config("protein_coding_ids_path"))
    
    expression_dict_WT_A2, isoforms_dict_WT_A2 = join_expression(expression_paths["WT_A2"], fas_lib.get_config("protein_coding_ids_path"))
    expression_dict_WT_A2B, isoforms_dict_WT_A2B = join_expression(expression_paths["WT_A2B"], fas_lib.get_config("protein_coding_ids_path"))
    
    expression_dict_WTNC_A2B, isoforms_dict_WTNC_A2B = join_expression(expression_paths["WTNC_A2B"], fas_lib.get_config("protein_coding_ids_path"))
    expression_dict_WTNC_A2, isoforms_dict_WTNC_A2 = join_expression(expression_paths["WTNC_A2"], fas_lib.get_config("protein_coding_ids_path"))
    
    expression_dict_P633L_A2, isoforms_dict_P633L_A2 = join_expression(expression_paths["P633L_A2"], fas_lib.get_config("protein_coding_ids_path"))
    expression_dict_P633L_A2B, isoforms_dict_P633L_A2B = join_expression(expression_paths["P633L_A2B"], fas_lib.get_config("protein_coding_ids_path"))
    
    
    DPF2_dict = {"R634Q_A2" : make_fas_graph_row(expression_dict_R634Q_A2, isoforms_dict_R634Q_A2, dist_matrix, test_ids_dict["DPF2"]),
                 "R634Q_A2B" : make_fas_graph_row(expression_dict_R634Q_A2B, isoforms_dict_R634Q_A2B, dist_matrix, test_ids_dict["DPF2"]),
                 "P633L_A2" : make_fas_graph_row(expression_dict_P633L_A2, isoforms_dict_P633L_A2, dist_matrix, test_ids_dict["DPF2"]),
                 "P633L_A2B" : make_fas_graph_row(expression_dict_P633L_A2B, isoforms_dict_P633L_A2B, dist_matrix, test_ids_dict["DPF2"]),
                 "WTNC_A2" : make_fas_graph_row(expression_dict_WTNC_A2, isoforms_dict_WTNC_A2, dist_matrix, test_ids_dict["DPF2"]),
                 "WTNC_A2B" : make_fas_graph_row(expression_dict_WTNC_A2B, isoforms_dict_WTNC_A2B, dist_matrix, test_ids_dict["DPF2"]),
                 "WT_A2" : make_fas_graph_row(expression_dict_WT_A2, isoforms_dict_WT_A2, dist_matrix, test_ids_dict["DPF2"]),
                 "WT_A2B" : make_fas_graph_row(expression_dict_WT_A2B, isoforms_dict_WT_A2B, dist_matrix, test_ids_dict["DPF2"])}
    
    with open(fas_lib.get_config("root_path") + "DPF2_FAS_graph_dict.json", 'w') as f:
        json.dump(DPF2_dict, f,  indent=4)
    
    ARNTL_dict = {"R634Q_A2" : make_fas_graph_row(expression_dict_R634Q_A2, isoforms_dict_R634Q_A2, dist_matrix, test_ids_dict["ARNTL"]),
                 "R634Q_A2B" : make_fas_graph_row(expression_dict_R634Q_A2B, isoforms_dict_R634Q_A2B, dist_matrix, test_ids_dict["ARNTL"]),
                 "P633L_A2" : make_fas_graph_row(expression_dict_P633L_A2, isoforms_dict_P633L_A2, dist_matrix, test_ids_dict["ARNTL"]),
                 "P633L_A2B" : make_fas_graph_row(expression_dict_P633L_A2B, isoforms_dict_P633L_A2B, dist_matrix, test_ids_dict["ARNTL"]),
                 "WTNC_A2" : make_fas_graph_row(expression_dict_WTNC_A2, isoforms_dict_WTNC_A2, dist_matrix, test_ids_dict["ARNTL"]),
                 "WTNC_A2B" : make_fas_graph_row(expression_dict_WTNC_A2B, isoforms_dict_WTNC_A2B, dist_matrix, test_ids_dict["ARNTL"]),
                 "WT_A2" : make_fas_graph_row(expression_dict_WT_A2, isoforms_dict_WT_A2, dist_matrix, test_ids_dict["ARNTL"]),
                 "WT_A2B" : make_fas_graph_row(expression_dict_WT_A2B, isoforms_dict_WT_A2B, dist_matrix, test_ids_dict["ARNTL"])}
    
    
    with open(fas_lib.get_config("root_path") + "ARNTL_FAS_graph_dict.json", 'w') as f:
        json.dump(ARNTL_dict, f,  indent=4)


