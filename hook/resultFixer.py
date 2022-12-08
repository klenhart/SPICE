# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:58:21 2022

@author: chris
"""

import argparse
import json

def parser_setup():
    """
    

    Returns
    -------
    Get paths and flags required to run.

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()


    parser.add_argument("-d", "--dir", type=str,
                        help="""Directory of the result folder to be fixed.""")

                        
    args = parser.parse_args()
    dir_path = args.dir

    return dir_path


def main():
    dir_path = parser_setup()
    with open(dir_path + "/result_config.json", "r") as f:
        info_dict = json.load(f)
    list_of_samples = sorted(list(info_dict["conditions"].keys()))
    
    non_expressed_genes_dict = dict()
    
    # Collected all non expressed genes
    for name in list_of_samples:
        non_expressed_genes_dict[name] = set()
        with open(dir_path + "/expression/expression_" + name + "_ENSv107.json", "r") as f:
            expression_dict = json.load(f)
        gene_list = list(expression_dict["expression"].keys())
        for gene in gene_list:
            not_expressed_flag = True
            replicate_list = list(expression_dict["expression"][gene].keys())
            for replicate in replicate_list:
                if expression_dict["expression"][gene][replicate] != 0:
                    not_expressed_flag = False
                    break
            if not_expressed_flag:
                non_expressed_genes_dict[name].add(gene)
      
    comparison_dict = dict()
    for i, sample1 in enumerate(list_of_samples):
        for j, sample2 in enumerate(list_of_samples):
            if i >= j:
                continue
            else:
              comparison_dict[sample1 + "@" + sample2 + "_all"] = list(non_expressed_genes_dict[sample1].union(non_expressed_genes_dict[sample2]))
    
    for comparison in list(comparison_dict.keys()):
        with open(dir_path + "/main_comparison/" + comparison + "/result_" + comparison + ".tsv", "r") as f:
            file = f.read()
            
            
        rows = file.split("\n")
        headerInfo = "\n".join(rows[0:3])
        header = rows[3]
        colNames = header.split("\t")
        colNames[7] = "conditionRMSD_lowerthan_interReplicateRMSD(max)(1)"
        colNames[8] = "conditionRMSD_lowerthan_interReplicateRMSD(avg+std)(1)"
        colNames[9] = "conditionRMSD_lowerthan_interReplicateRMSD(max)(2)"
        colNames[10] = "conditionRMSD_lowerthan_interReplicateRMSD(avg+std)(2)"
        colNames.append("notExpressedFlag")
        header = "\t".join(colNames)
        
        output = [ headerInfo, header ]
        
        for row in rows[4:]:
            entries = row.split("\t")
            if entries[0] in comparison_dict[comparison]:
                entries.append("1")
            else:
                entries.append("0")
            new_row = "\t".join(entries)
            output.append(new_row)
            
        finished_new_result = "\n".join(output)
        
        with open(dir_path + "/main_comparison/" + comparison + "/result_" + comparison + ".tsv", "w") as f:
            f.write(finished_new_result)
            
if __name__ == "__main__":
    main()
    