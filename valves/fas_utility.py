#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################################################
# Copyright (C) 2022 Christian, Blümel, Julian Dosch
#
# This file is part of grand-trumpet.
#
#  grand-trumpet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  grand-trumpet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with grand-trumpet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
"""
Created on Thu Aug 18 11:34:52 2022

@author: chrisbl
"""

"""
Move all minor tools from all modules to this.
"""

import itertools
import json

def tsv_collection_maker(header_dict, fas_lib):
    """
    Generates a pairings_tsv.json file in the root_path that contains all pairings.tsv files.
    Parameters
    ----------
    header_dict : dict.keys() == [str], dict.values() == [str]
        dictionary containing all headers indexed by their ENS gene ID
    fas_lib : FAS Library class object
        path to the root of specific species library (e.g. /home/FAS_library/homo_sapiens/release_107/)
    Returns
    -------
    None.
    """
    tsv_dict = dict()
    for gene_id in header_dict.keys():
        tsv_dict[gene_id] = ""
        pairs = itertools.product(header_dict[gene_id], header_dict[gene_id])
        pairs = [pair for pair in pairs if pair[0] <= pair[1]]
        for header_1, header_2 in pairs:
            tsv_dict[gene_id] += header_1 + "\t" + header_2 + "\n"
    with open(fas_lib.get_config("pairings_tsv_json_path"), 'w') as fp:
        json.dump(tsv_dict, fp,  indent=4)

def find_max_tsl(fas_lib, protein_id_list):
    if len(protein_id_list) == 0:
        return "7"
    tuple_list = tsv_to_tuple_list(fas_lib.get_config("protein_coding_ids_path"))
    return str(max([int(entry[-2]) for entry in tuple_list if entry[1] in protein_id_list]))

def tuple_list_to_tsv(tuple_list):
    tsv = ""
    for str_x, str_y in tuple_list:
        tsv += str_x + "\t" + str_y + "\n"
    return tsv

def triple_list_to_tsv(triple_list):
    tsv = ""
    for str_x, str_y, str_z in triple_list:
        tsv += str_x + "\t" + str_y + "\t" + str_z + "\n"
    return tsv

def quadruple_list_to_tsv(quadruple_list):
    tsv = ""
    for str_w, str_x, str_y, str_z in quadruple_list:
        tsv += str_w + "\t" + str_x + "\t" + str_y + "\t" + str_z + "\n"
    return tsv

def quintuple_list_to_tsv(quintuple_list):
    tsv = ""
    for str_v, str_w, str_x, str_y, str_z in quintuple_list:
        tsv += str_v + "\t" + str_w + "\t" + str_x + "\t" + str_y + "\t" + str_z + "\n"
    return tsv

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
    return [protein_id for _, protein_id, transcript_id, _ in ids_triple if transcript_id in transcript_ids]



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

def longest_common_prefix(str_list):
    prefix = ""
    i = 0
    flag_break = False
    while i < len(str_list[0]):
        for string in str_list[1:]:
            if str_list[0][i] != string[i]:
                flag_break = True
                break
        if flag_break:
            break
        else:
            prefix += str_list[0][i]
            i += 1
    return prefix
                

def load_gene_ids_txt(gene_ids_path):
    with open(gene_ids_path, "r") as f:
        gene_ids = f.read()
        gene_ids = gene_ids.split("\n")
        gene_ids = [gene_id for gene_id in gene_ids if len(gene_id) > 0]
        if gene_ids[0] == "gene":
            gene_ids = gene_ids[1:]
    return gene_ids