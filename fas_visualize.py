#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:29:49 2022

@author: chrisbl
"""

"""
Move the actual visualization part of option
"--visualize",
from fas_handler.py to this.
"""

# Standard modules
import argparse
import sys

# self made modules
import library_class
import fas_utility
import fas_polygon

def parser_setup():
    """
    Reads the parser input.

    Returns
    -------
    Get run options if either to create a tsv file, to delete it again or to concatenate FAS results into a main file.
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library.")
    
    parser.add_argument("-g", "--gene", type=str,
                        help="""Ensembl gene ID of the gene that shall be visualized.""")
    
    parser.add_argument("-f", "--outFormat", type=str, default="svg",
                        help="""Define the output format. Default: svg""")
    
    parser.add_argument("-q", "--quick", action="store_true",
                        help="Visualize from two FAS polygon files instead of one precalculated FAS file.")
    
    parser.add_argument("-p", "--path", nargs="+", action="append",
                        help="""Path either to the polygon files that shall be compared or to the already calculated comparison file.
                        If the former is the case you must also choose the option --quick""")
    
    


    args = parser.parse_args()

    config_path = args.config
    gene_id = args.gene
    path_list = args.path[0]
    outFormat = args.outFormat
    flag_quick = args.quick

    return config_path, gene_id, path_list, outFormat, flag_quick


def main():
    # Get parser input
    config_path, gene_id, path_list, outFormat, flag_quick = parser_setup()
    
    #Import library settings.
    fas_lib = library_class.Library(config_path, False)
    
    # Do not use pre calculated file.
    if flag_quick:
        fas_polygon.visualize_fas_polygon(path_list, fas_lib, gene_id, outFormat)
    
    # Use precalculated file.
    else:
        # Generate filenames.
        filename = fas_utility.get_name(path_list[0]).split("/")
        filename_scaled = filename[:-1] + [gene_id + "_scaled." + outFormat]
        filename_unscaled = filename[:-1] + [gene_id + "_unscaled." + outFormat]
        
        filename_scaled = "/".join(filename_scaled)
        filename_unscaled = "/".join(filename_unscaled)
        
        # Import pregenerated file.
        with open(path_list[0], "r") as f:
                all_polygons = f.read().split("\n")
                all_polygons = [ polygon.split("\t") for polygon in all_polygons ]
                all_polygons = [ entry for entry in all_polygons if entry[0] == gene_id]
        # Check if the file actually contains the gene. This should technically always work if the library is intact.
        if len(all_polygons) > 0:
            gene_id, sample_names, categories, unscaled_expression, scaled_expression, unscaled_rmsd, scaled_rmsd = all_polygons[0]
            sample_names = sample_names.split(";")
            categories = categories.split(";")
            # Check if there is no expression for any isoforms in both graphs.
            if categories == "":
                raise Exception(path_list[0], "has no entry of gene", gene_id)
                sys.exit(1)
            else:
                unscaled1, unscaled2 = unscaled_expression.split(";")
                scaled1, scaled2 = scaled_expression.split(";")
                unscaled1 = [ float(entry) for entry in unscaled1.split(":") ]
                unscaled2 = [ float(entry) for entry in  unscaled2.split(":") ]
                scaled1 = [ float(entry) for entry in  scaled1.split(":") ]
                scaled2 = [ float(entry) for entry in  scaled2.split(":") ]
                
                filepath_unscaled = fas_lib.get_config("root_path") + "/pictures/" + filename[-1] + "/" + filename_unscaled
                filepath_scaled = fas_lib.get_config("root_path") + "/pictures/" + filename[-1] + "/" + filename_scaled
                # Actually draw the graphs.
                fas_polygon.make_graph(fas_lib,
                                       gene_id,
                                       sample_names,
                                       categories,
                                       [unscaled1, unscaled2],
                                       unscaled_rmsd,
                                       filepath_unscaled)
                fas_polygon.make_graph(fas_lib,
                                       gene_id,
                                       sample_names,
                                       categories,
                                       [scaled1, scaled2],
                                       scaled_rmsd,
                                       filepath_scaled)


if __name__ == "__main__":
    main()