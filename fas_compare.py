#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:31:47 2022

@author: chrisbl
"""

import argparse
import os

import fas_polygon as poly
import ensembl_access
import library_class
import fas_utility

"""
Move the generation of a sample comparison file part of the
"--visualize"
option from fas_handler.py to this
"""

def extract_all_graph(fas_lib, path, exempt=["ENSG00000155657"]):
    fas_graphs_dict = dict()
    gene_ids = ensembl_access.load_gene_ids_txt(fas_lib.get_config("gene_ids_path"))
    gene_ids = [ gene_id for gene_id in gene_ids if gene_id not in exempt]
    for gene_id in gene_ids:
        fas_graphs_dict[gene_id] = poly.extract_graph(path, gene_id)
    return fas_graphs_dict


def generate_comparison(polygon_dict_list, path_list, fas_lib, name_list):
    title = "x".join(name_list) + "/"
    pictures_path = fas_lib.get_config("root_path") + "pictures/"
    comparison_path = pictures_path + title
    file_path = comparison_path + "polygonFAS_" + "x".join(name_list) + ".tsv"
    
    if not os.path.exists(pictures_path):
        os.makedirs(pictures_path)
    if not os.path.exists(comparison_path):
        os.makedirs(comparison_path)
    
    polygon_dict1, polygon_dict2 = polygon_dict_list
    output = "gene_id\tsample_names\tprot_id\tunscaled_expression\tscaled_expression\tunscaled_rmsd\tscaled_rmsd\n"
    for gene_name in polygon_dict1.keys():
        polygon_list = [ polygon_dict1[gene_name], polygon_dict2[gene_name]]
        polygon_dict = poly.prepare_polygon_pair_visual(polygon_list, name_list, gene_name)

        output_row_list = [polygon_dict["gene_id"],
                           ";".join(polygon_dict["name_list"])]
        unscaled_expr_list = []
        scaled_expr_list = []
        output_row_list.append(";".join(polygon_dict["categories"]))
        for i, name in enumerate(name_list):
            print(polygon_dict["unscaled_rel_expr"])
            print(2.0*"")
            unscaled_expr_list.append(":".join([ str(value) for value in polygon_dict["unscaled_rel_expr"]]))
            scaled_expr_list.append(":".join([ str(value) for value in polygon_dict["scaled_rel_expr"]]))
        output_row_list.append(";".join(unscaled_expr_list))
        output_row_list.append(";".join(scaled_expr_list))
        output_row_list.append(str(polygon_dict["unscaled_rmsd"]))
        output_row_list.append(str(polygon_dict["scaled_rmsd"]))

        output_row = "\t".join(output_row_list)
        
        output += output_row + "\n"
    with open(file_path, "w") as f:
        f.write(output)
    return file_path

def sort_by_rmsd(fas_lib, path, flag_more_than_2=True):
    output = "gene_id\tsample_names\tprot_id\tunscaled_expression\tscaled_expression\tunscaled_rmsd\tscaled_rmsd\n"
    new_path = path[:-4] + "_sorted.tsv"
    with open(path, "r") as f:
        file = f.read()
        comparisons = file.split("\n")
    comparisons = comparisons[1:]
    comparisons = [ comparison.split("\t") for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if len(comparison) == 7:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = [ comparison[:-2] + [float(comparison[-2]), float(comparison[-1])] for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if comparison[-2] < 1 and comparison[-2] > 0 and comparison[-2] != comparison[-1]:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = sorted(comparisons, key = lambda x: (-x[-2], -x[-1]))
    comparisons = [ comparison[:-2] + [str(comparison[-2]), str(comparison[-1])] for comparison in comparisons ]
    comparisons = [ "\t".join(comparison) for comparison in comparisons ]
    file = "\n".join(comparisons)
    file = output + file
    with open(new_path, "w") as f:
        f.write(file)


def parser_setup():
    """
    Reads the parser input.

    Returns
    -------
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library.")

    parser.add_argument("-i", "--input", nargs="+", action="append",
                        help="Paths to two fas_polygon files for which all comparisons shall be calculated.")
    
    args = parser.parse_args()
    
    config_path = args.config
    path_list = args.input[0]

    return config_path, path_list
    
def main():
    config_path, path_list = parser_setup()
    fas_lib = library_class.Library(config_path, False)
    path_list = sorted(path_list)
    fas_graphs_dict_list = [ extract_all_graph(fas_lib, path) for path in path_list ]
    name_list = [ fas_utility.get_name(path) for path in path_list ]
    file_path = generate_comparison(fas_graphs_dict_list, path_list, fas_lib, name_list)
    sort_by_rmsd(fas_lib, file_path, flag_more_than_2=True)

if __name__ == "__main__":
    main()
   