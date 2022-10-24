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
Created on Thursday October 20 10:41:03 2022

@author: chrisbl
"""

import json
import pyranges as pr

tsl_dict = {"1": "1", "2": "2", "3": "3",
    "4": "4", "5": "5", "n": "6", "N": "6"}


def get_all_domains(path):
    with open(path, "r") as fp:
        domain_dict = json.load(fp)
    return list(domain_dict.keys())


def make_ensembl_isoform_dataframe(ensembl_path):
    gtf = pr.read_gtf(ensembl_path, as_df=True)
    filtered = [(feature, str(p_id), str(t_id), g_id,
                 g_biotype, str(t_biotype), tsl_dict[str(tsl)[0]],
                 str(tag)) for g_id, p_id, t_id,
                          g_biotype, t_biotype, tsl,
                          tag, feature in zip(gtf['gene_id'], gtf['protein_id'], gtf['transcript_id'],
                                              gtf['gene_biotype'], gtf['transcript_biotype'], gtf['transcript_support_level'],
                                              gtf['tag'], gtf['Feature']) if g_biotype in ["protein_coding",
                                                                                           "protein_coding_LoF"] and t_biotype in ["protein_coding",
                                                                                                                                   "protein_coding_LoF"] and feature in ["gene", "transcript"]]
    filtered = [ entry for entry in filtered if entry[1:3] == ("nan", "nan") or entry[1] != "nan" and entry[2] != "nan"]                                                                                                                           
    tuple_list_to_tsv(filtered, ensembl_path[:-3] + "tsv", ['Feature', 'protein_id', 'transcript_id', 
                                                            'gene_id', 'gene_biotype', 'transcript_biotype', 
                                                            'transcript_support_level', 'tag'])


def tuple_list_to_tsv(tuple_list, path, col_names):
    row_list = ["\t".join(col_names)]
    row_list += ["\t".join(entry) for entry in tuple_list]
    tsv= "\n".join(row_list)
    with open(path, "w") as f:
        f.write(tsv)

def main():
    make_ensembl_isoform_dataframe("/home/chrisbl/backup/Bachelor/Ensembl_stats/test.gtf")

if __name__ == "__main__":
    main()