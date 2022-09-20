#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:34:52 2022

@author: chrisbl
"""

"""
Move all minor tools from all modules to this.
"""

def tsv_to_tuple_list(path):
    tuple_list = []
    with open(path, "r") as f:
        tsv = f.read()
        tsv = tsv.split("\n")
        tsv = [ entry for entry in tsv if len(entry) > 0 ]
    tuple_list = [ tuple(entry.split("\t")) for entry in tsv ]
    return tuple_list

def transcript_id_to_protein_id(transcript_ids, ids_path):
    ids_triple = tsv_to_tuple_list(ids_path)
    return [protein_id for _, protein_id, transcript_id in ids_triple if transcript_id in transcript_ids]

def build_fasta_from_ids_and_fasta(id_list, input_fasta_path, output_fasta_path):
    output_list = []
    with open(input_fasta_path, "r") as f:
        input_fasta = f.read()
        input_fasta = input_fasta.split(">")
    for protein_id in id_list:
        for fasta_entry in input_fasta:
            if protein_id in fasta_entry:
                output_list.append(fasta_entry)
                break
    output_fasta = ">".join(output_list)
    with open(output_fasta_path, "w") as f:
        f.write(output_fasta)
    

def start_stop_range(length, n):
    """
    Splits length into equal chunks of size n.
    
    Parameters
    ----------
    length : int
        length to be split apart.
    n : int
        size of chunks.

    Yields
    ------
    generator of (int, int)
        Splits length into equal chunks of size n.
    """
    for i in range(0, length, n):
        yield (i, min(i+n-1, length))

def get_name(path):
    start = path.index("polygonFAS_") + 11
    end = path.index(".tsv")
    return path[start:end]

def make_request_data(id_list, id_index=1):
    id_list = [ids[id_index] for ids in id_list]
    request_data = '{ "ids" : ['
    for entry in id_list:
        request_data += '"' + entry + '", '
    request_data = request_data[:-2]
    request_data += ' ] }'
    return request_data

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def load_gene_ids_txt(gene_ids_path):
    with open(gene_ids_path, "r") as f:
        gene_ids = f.read()
        gene_ids = gene_ids.split("\n")
        gene_ids = [gene_id for gene_id in gene_ids if len(gene_id) > 0]
        if gene_ids[0] == "gene":
            gene_ids = gene_ids[1:]
    return gene_ids