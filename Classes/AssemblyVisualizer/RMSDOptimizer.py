#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  RMSDOptimizer is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RMSDOptimizer is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import numpy as np
from typing import Dict, List, Any, Set
from scipy.optimize import minimize
import argparse
import copy
import json
import os
import random


class RMSDOptimizer:

    def __init__(self,
                 distance_dict: Dict[str, Dict[str, float]],
                 info_dict: Dict[str, Any],
                 protein_coding_flag: bool = True,
                 complete_flag: bool = True,
                 pre_ewfd_flag: bool = False):

        self.matrix = RMSDOptimizer.distance_dict_to_matrix(distance_dict,
                                                            protein_coding_flag,
                                                            complete_flag,
                                                            pre_ewfd_flag,
                                                            info_dict)

        if info_dict["_id"] == "ENSG00000163081":
            print(self.matrix)
            print(distance_dict)

        self.max_rmsd: float = RMSDOptimizer.extract_max_rmsd(self.matrix)

        if info_dict["_id"] == "ENSG00000163081":
            print(self.max_rmsd)
            print("---")

    @staticmethod
    def extract_max_rmsd(matrix):
        vector: List[float] = list()
        for i, _ in enumerate(matrix):
            for j, _ in enumerate(matrix[i]):
                if i != j:
                    vector.append((matrix[i][j] + matrix[j][i]) / 2)
        complement: List[float] = [1-entry for entry in vector]
        if len(complement) == 0:
            return 0.0
        else:
            return max(complement)

    @staticmethod
    def distance_dict_to_matrix(distance_dict: Dict[str, Dict[str, float]],
                                protein_coding_flag: bool,
                                complete_flag: bool,
                                pre_ewfd_flag: bool,
                                info_dict: Dict[str, Any]) -> np.array:
        raw_matrix: List[List[float]] = list()
        complete_flag = complete_flag and pre_ewfd_flag
        protein_coding_flag = protein_coding_flag and pre_ewfd_flag
        for protein_id_outer in distance_dict.keys():
            if complete_flag and "incomplete" in info_dict["transcripts"][protein_id_outer]["tags"]:
                continue
            elif protein_coding_flag and info_dict["transcripts"][protein_id_outer]["biotype"] != "protein_coding":
                continue
            raw_matrix.append(list())
            for protein_id_inner in distance_dict[protein_id_outer].keys():
                if complete_flag and "incomplete" in info_dict["transcripts"][protein_id_inner]["tags"]:
                    continue
                elif protein_coding_flag and info_dict["transcripts"][protein_id_inner]["biotype"] != "protein_coding":
                    continue
                raw_matrix[-1].append(distance_dict[protein_id_outer][protein_id_inner])
        return np.array(raw_matrix)


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-l",
                        "--library",
                        type=str,
                        action="store",
                        help="Path to library.")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        action="store",
                        help="Name of the output file.")
    parser.add_argument("-a",
                        "--already",
                        type=str,
                        action="store",
                        help="Path to tsv file containing already calced genes.")
    parser.add_argument("-m",
                        "--mode",
                        action="store_true",
                        help="""Alternative version that truly reduces the library by non_coding and incomplete
                        transcript.""")

    argument_dict: Dict[str, Any] = vars(parser.parse_args())

    if argument_dict["mode"] is None:
        argument_dict["mode"] = False

    already_calced: Set[str] = set()

    if argument_dict["already"] is not None:
        with open(argument_dict["already"], "r") as f:
            output_file: List[str] = f.read().split("\n")
            output_list: List[List[str]] = [entry.split("\t") for entry in output_file]
            for entry in output_list[1:]:
                already_calced.add(entry[0])
    else:
        output_list: List[List[str]] = [["geneid", "all", "no_incomplete", "no_non_coding", "no_both"]]

    with open(os.path.join(argument_dict["library"], "transcript_data", "transcript_info.json"), "r") as f:
        info_dict = json.load(f)

    with open(os.path.join(argument_dict["library"], "fas_data", "fas_index.json"), "r") as f:
        fas_file_list = list(set(json.load(f).values()))

    total: int = len(fas_file_list)

    for i, fas_file in enumerate(fas_file_list):
        print(str(i) + "/" + str(total), fas_file)
        with open(os.path.join(argument_dict["library"], "fas_data", "fas_scores", fas_file), "r") as f:
            distance_dicts: Dict[str, Any] = json.load(f)
        for gene_id in distance_dicts.keys():
            if gene_id not in already_calced:
                stages: List[bool] = [False, False, False, False]

                optimizer = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                          info_dict[gene_id],
                                          False,
                                          False,
                                          argument_dict["mode"])
                stages[0] = True
                optimizer_no_incomplete = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                        info_dict[gene_id], False, True,
                                                        argument_dict["mode"])
                stages[1] = True
                optimizer_no_non_coding = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                        info_dict[gene_id], True, False,
                                                        argument_dict["mode"])
                stages[2] = True
                optimizer_both = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                               info_dict[gene_id],
                                               True,
                                               True,
                                               argument_dict["mode"])
                stages[3] = True
                output_list.append([gene_id,
                                    str(round(optimizer.max_rmsd, 2)),
                                    str(round(optimizer_no_incomplete.max_rmsd, 2)),
                                    str(round(optimizer_no_non_coding.max_rmsd, 2)),
                                    str(round(optimizer_both.max_rmsd, 2))])

        with open(argument_dict["outfile"], "w") as f:
            f.write("\n".join(["\t".join(entry) for entry in output_list]))
    #  print("\n".join(["\t".join(entry) for entry in output_list]))
    with open(argument_dict["outfile"], "w") as f:
        f.write("\n".join(["\t".join(entry) for entry in output_list]))


if __name__ == "__main__":
    main()
