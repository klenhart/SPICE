#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:13:06 2022

@author: chrisbl
"""

import plotly.graph_objects as go
import json

from library_class import Library

########
def main():
    fas_lib = Library("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/config.tsv", False)
    
    with open(fas_lib.get_config("root_path") + "ARNTL_FAS_graph_dict.json", "r") as f: 
        ARNTL_dict = json.load(f)
    
    with open(fas_lib.get_config("root_path") + "DPF2_FAS_graph_dict.json", "r") as f: 
        DPF2_dict = json.load(f)
    
    DPF2_transcript_ids = set()
    
    DPF2_vector_dict = dict()
    DPF2_fpkm_dict = dict()
    
    for sample in DPF2_dict.keys():
        row = DPF2_dict[sample]
        gene_id, transcripts, fpkm = row.split("\t")
    
        DPF2_fpkm_dict[sample] = float(fpkm)
        DPF2_vector_dict[sample] = dict()
        
        protein_id_expression_tuple_list = [entry.split(":") for entry in transcripts.split(";")]
    
        for protein_id, rel_expr in protein_id_expression_tuple_list:
            DPF2_transcript_ids.add(protein_id)
            DPF2_vector_dict[sample][protein_id] = float(rel_expr)
    
    DPF2_transcript_ids = list(DPF2_transcript_ids)
    
    #######
    
    ARNTL_transcript_ids = set()
    
    ARNTL_vector_dict = dict()
    ARNTL_fpkm_dict = dict()
    
    for sample in ARNTL_dict.keys():
        row = ARNTL_dict[sample]
        gene_id, transcripts, fpkm = row.split("\t")
    
        ARNTL_fpkm_dict[sample] = float(fpkm)
        ARNTL_vector_dict[sample] = dict()
        
        protein_id_expression_tuple_list = [entry.split(":") for entry in transcripts.split(";")]
    
        for protein_id, rel_expr in protein_id_expression_tuple_list:
            ARNTL_transcript_ids.add(protein_id)
            ARNTL_vector_dict[sample][protein_id] = float(rel_expr)
    
    ARNTL_transcript_ids = list(ARNTL_transcript_ids)
    
    #########
    
    #ARNTL
    #Example for 
    make_graph_unscaled("ARNTL", ["WT_A2B","R634Q_A2B"], ARNTL_transcript_ids, ARNTL_vector_dict, fas_lib)
    # Example for unscaled while several samples dont work well without scaling
    make_graph_unscaled("ARNTL", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, fas_lib)
    # Same graph but scaled
    make_graph_scaled("ARNTL", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WTNC_A2", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WTNC_A2B", "R634Q_A2B"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WT_A2", "WTNC_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["R634Q_A2", "P633L_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    
    #DPF2
    #Example for 
    make_graph_unscaled("DPF2", ["WT_A2B","R634Q_A2B"], DPF2_transcript_ids, DPF2_vector_dict, fas_lib)
    # Example for unscaled while several samples dont work well without scaling
    make_graph_unscaled("DPF2", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, fas_lib)
    # Same graph but scaled
    make_graph_scaled("DPF2", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["P633L_A2", "P633L_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["WTNC_A2", "R634Q_A2B"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["WTNC_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)

def make_graph_unscaled(gene_name, figure_samples, categories, vector_dict, fas_lib):
    title = gene_name + "_unscaled_" + "x".join([ name for name in figure_samples])
    fig = go.Figure()        
    vecs = dict()

    for sample in figure_samples:
        vecs[sample] = []
        for prot_id in categories:
            if prot_id in vector_dict[sample].keys():
                vecs[sample].append(vector_dict[sample][prot_id])
            else:
                vecs[sample].append(0.0)
        fig.add_trace(go.Scatterpolar(
            r=vecs[sample],
            theta=categories,
            fill="toself",
            name=sample
            ))
        
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True
        )
    fig.show()
    fig.write_image(file=fas_lib.get_config("root_path") + title + ".svg")
    
def make_graph_scaled(gene_name, figure_samples, categories, vector_dict,fpkm_dict , fas_lib):
    title = gene_name + "_scaled_" + "x".join([ name for name in figure_samples])
    fig = go.Figure()        
    vecs = dict()
    
    # Find largest fpkm
    max_fpkm = 0
    for sample in figure_samples:
        if fpkm_dict[sample] > max_fpkm:
            max_fpkm = fpkm_dict[sample]
    
    # Calculate Scales in comparison to largest fpkm.
    scales = dict()    
    for sample in figure_samples:
        scales[sample] = fpkm_dict[sample] / max_fpkm
    
    for sample in figure_samples:
        vecs[sample] = []
        for prot_id in categories:
            if prot_id in vector_dict[sample].keys():
                vecs[sample].append(vector_dict[sample][prot_id])
            else:
                vecs[sample].append(0.0)
        vecs[sample] = scale_list(scales[sample], vecs[sample])
        fig.add_trace(go.Scatterpolar(
            r=vecs[sample],
            theta=categories,
            fill="toself",
            name=sample
            ))
        
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True
        )
    fig.show()
    fig.write_image(file=fas_lib.get_config("root_path") + title + ".svg")

def scale_list(scale_factor, float_list):
    output_list = []
    for value in float_list:
        output_list.append(value * scale_factor)
    return output_list

if __name__ == "__main__":
    main()

