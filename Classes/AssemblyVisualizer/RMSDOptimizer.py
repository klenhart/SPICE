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


class RMSDOptimizerAlt:

    def __init__(self,
                 distance_dict: Dict[str, Dict[str, float]],
                 info_dict: Dict[str, Any],
                 protein_coding_flag: bool = True,
                 complete_flag: bool = True):
        self.matrix = RMSDOptimizerAlt.distance_dict_to_matrix(distance_dict,
                                                               protein_coding_flag,
                                                               complete_flag,
                                                               info_dict)
        self.length: int = len(self.matrix)
        if self.length == 0:
            self.result = None
        else:
            self.bounds = [(0, 1)] * (2 * self.length)
            while True:
                new_result = self.__optimize__()
                if self.test_new_result_for_corruption(new_result):
                    self.result = new_result
                    break
            count: int = 0
            while count != 5:
                new_result = self.__optimize__()
                if self.test_new_result_for_corruption(new_result):
                    if -self.objective_function(new_result.x) > -self.objective_function(self.result.x):
                        self.result = new_result
                        count = 0
                    else:
                        count += 1

    def objective_function(self, combined_vector):
        v = combined_vector[:self.length]
        w = combined_vector[self.length:]
        expressed_indices = list()
        for i, x in enumerate(v):
            if x > 0 or w[i] > 0:
                expressed_indices.append(i)
        if len(expressed_indices) == 0:
            return 0.0
        diff = np.dot(self.matrix, v) - np.dot(self.matrix, w)
        expressed_diff_list = list()
        for i in expressed_indices:
            expressed_diff_list.append(diff[i])
        expressed_diff = np.array(expressed_diff_list)
        return -np.sqrt(np.sum(expressed_diff ** 2) / len(expressed_diff))

    def constraint_function_v(self, combined_vector):
        v = combined_vector[:self.length]
        return np.sum(v) - 1.0

    def constraint_function_w(self, combined_vector):
        w = combined_vector[self.length:]
        return np.sum(w) - 1.0

    def __optimize__(self):
        v_init = RMSDOptimizerAlt.make_random_vector(self.length)
        w_init = RMSDOptimizerAlt.make_random_vector(self.length)
        x0 = np.concatenate((v_init, w_init))
        problem = {
            'fun': self.objective_function,
            'x0': x0,
            'constraints': [{'type': 'eq', 'fun': self.constraint_function_v},
                            {'type': 'eq', 'fun': self.constraint_function_w}],
            'method': 'SLSQP',
            'bounds': self.bounds
        }
        return minimize(**problem)

    def get_results(self):
        return self.result

    @staticmethod
    def make_random_vector(n: int) -> np.array:
        vector = np.zeros(n)
        remain = 100
        number_pool = list()

        for i in range(n):
            if remain == 0:
                number_pool.append(0.0)
            elif i == n-1:
                number_pool.append(round(remain*0.01, 2))
            else:
                x = random.randint(1, remain)
                number_pool.append(round(x*0.01, 2))
                remain = remain - x

        for i, _ in enumerate(vector):
            number = random.choice(number_pool)
            number_pool.remove(number)
            vector[i] = number
        return vector

    @staticmethod
    def distance_dict_to_matrix(distance_dict: Dict[str, Dict[str, float]],
                                protein_coding_flag: bool,
                                complete_flag: bool,
                                info_dict: Dict[str, Any]) -> np.array:
        raw_matrix: List[List[float]] = list()
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
                raw_matrix[-1].append(1 - distance_dict[protein_id_outer][protein_id_inner])
        return np.array(raw_matrix)

    @staticmethod
    def is_non_coding(transcript_id: str) -> bool:
        if "T" in [transcript_id[3], transcript_id[6]]:
            return True
        elif transcript_id[3:8] == "noORF":
            return True
        elif transcript_id[3:8] == "noDia":
            return True
        else:
            return False

    @staticmethod
    def is_incomplete(transcript_id: str, info_dict: Dict[str, Any]):
        return "incomplete" in info_dict["transcripts"][transcript_id]["tags"]

    def calc_rmsd(self):
        if self.result is None:
            return 0
        else:
            return -self.objective_function(self.result.x)

    def test_new_result_for_corruption(self, new_result):
        if abs(self.objective_function(new_result.x)) <= 1:
            return True
        elif 0.98 <= sum(new_result.x[:self.length]) <= 1 and 0.98 <= sum(new_result.x[self.length:]) <= 1:
            return True
        else:
            return False


class RMSDOptimizer:

    def __init__(self,
                 distance_dict: Dict[str, Dict[str, float]],
                 info_dict: Dict[str, Any],
                 protein_coding_flag: bool = True,
                 complete_flag: bool = True):
        self.compliant_indices_set = RMSDOptimizer.check_distance_dict(distance_dict,
                                                                       info_dict,
                                                                       protein_coding_flag,
                                                                       complete_flag)
        self.matrix = RMSDOptimizer.distance_dict_to_matrix(distance_dict)
        self.length: int = len(self.matrix)
        if self.length == 0:
            self.result = None
        else:
            self.bounds = [(0, 1)] * (2 * self.length)
            self.result = self.__optimize__()
            count: int = 0
            while count != 5:
                new_result = self.__optimize__()
                if self.test_new_result_for_corruption(new_result):
                    if -self.objective_function(new_result.x) > -self.objective_function(self.result.x):
                        self.result = new_result
                        count = 0
                    else:
                        count += 1

    def test_new_result_for_corruption(self, new_result):
        if abs(self.objective_function(new_result.x)) <= 1:
            return True
        elif 0.98 <= sum(new_result.x[:self.length]) <= 1 and 0.98 <= sum(new_result.x[self.length:]) <= 1:
            return True
        else:
            return False

    def objective_function(self, combined_vector):
        v = combined_vector[:self.length]
        w = combined_vector[self.length:]
        expressed_indices = list()
        for i, x in enumerate(v):
            if (x > 0 or w[i] > 0) and i in self.compliant_indices_set:
                expressed_indices.append(i)
        if len(expressed_indices) == 0:
            return 0.0
        diff = np.dot(self.matrix, v) - np.dot(self.matrix, w)
        expressed_diff_list = list()
        for i in expressed_indices:
            expressed_diff_list.append(diff[i])
        expressed_diff = np.array(expressed_diff_list)
        return -np.sqrt(np.sum(expressed_diff ** 2) / len(expressed_diff))

    def constraint_function_v(self, combined_vector):
        v = combined_vector[:self.length]
        return np.sum(v) - 1.0

    def constraint_function_w(self, combined_vector):
        w = combined_vector[self.length:]
        return np.sum(w) - 1.0

    def __optimize__(self):
        v_init = RMSDOptimizer.make_random_vector(self.length)
        w_init = RMSDOptimizer.make_random_vector(self.length)
        x0 = np.concatenate((v_init, w_init))
        problem = {
            'fun': self.objective_function,
            'x0': x0,
            'constraints': [{'type': 'eq', 'fun': self.constraint_function_v},
                            {'type': 'eq', 'fun': self.constraint_function_w}],
            'method': 'SLSQP',
            'bounds': self.bounds
        }
        return minimize(**problem)

    def get_results(self):
        return self.result

    @staticmethod
    def make_random_vector(n: int) -> np.array:
        vector = np.zeros(n)
        remain = 100
        number_pool = list()

        for i in range(n):
            if remain == 0:
                number_pool.append(0.0)
            elif i == n-1:
                number_pool.append(round(remain*0.01, 2))
            else:
                x = random.randint(1, remain)
                number_pool.append(round(x*0.01, 2))
                remain = remain - x

        for i, _ in enumerate(vector):
            number = random.choice(number_pool)
            number_pool.remove(number)
            vector[i] = number
        return vector

    @staticmethod
    def check_distance_dict(distance_dict: Dict[str, Dict[str, float]],
                            info_dict: Dict[str, Any],
                            protein_coding_flag,
                            complete_flag):
        compliant_indices_set = set()
        if protein_coding_flag or complete_flag:
            if protein_coding_flag and complete_flag:
                for i, key in enumerate(distance_dict.keys()):
                    if not (RMSDOptimizer.is_non_coding(key) or RMSDOptimizer.is_incomplete(key, info_dict)):
                        compliant_indices_set.add(i)
            elif protein_coding_flag:
                for i, key in enumerate(distance_dict.keys()):
                    if not RMSDOptimizer.is_non_coding(key):
                        compliant_indices_set.add(i)
            elif complete_flag:
                for i, key in enumerate(distance_dict.keys()):
                    if not RMSDOptimizer.is_incomplete(key, info_dict):
                        compliant_indices_set.add(i)
        else:
            for i in range(len(distance_dict.keys())):
                compliant_indices_set.add(i)
        return compliant_indices_set

    @staticmethod
    def distance_dict_to_matrix(distance_dict: Dict[str, Dict[str, float]]) -> np.array:
        raw_matrix: List[List[float]] = list()
        for protein_id_outer in distance_dict.keys():
            raw_matrix.append(list())
            for protein_id_inner in distance_dict[protein_id_outer].keys():
                raw_matrix[-1].append(1 - distance_dict[protein_id_outer][protein_id_inner])
        return np.array(raw_matrix)

    @staticmethod
    def is_non_coding(transcript_id: str) -> bool:
        if "T" in [transcript_id[3], transcript_id[6]]:
            return True
        elif transcript_id[3:8] == "noORF":
            return True
        elif transcript_id[3:8] == "noDia":
            return True
        else:
            return False

    @staticmethod
    def is_incomplete(transcript_id: str, info_dict: Dict[str, Any]):
        return "incomplete" in info_dict["transcripts"][transcript_id]["tags"]

    def calc_rmsd(self):
        if self.result is None:
            return 0
        else:
            return -self.objective_function(self.result.x)


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

    argument_dict: Dict[str, str] = vars(parser.parse_args())

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

    if argument_dict["mode"]:
        print("Running alt mode")
        for i, fas_file in enumerate(fas_file_list):
            print(str(i) + "/" + str(total), fas_file)
            with open(os.path.join(argument_dict["library"], "fas_data", "fas_scores", fas_file), "r") as f:
                distance_dicts: Dict[str, Any] = json.load(f)
            for gene_id in distance_dicts.keys():
                if gene_id not in already_calced:
                    stages: List[bool] = [False, False, False, False]
                    try:
                        optimizer = RMSDOptimizerAlt(copy.deepcopy(distance_dicts[gene_id]),
                                                     info_dict[gene_id],
                                                     False,
                                                     False)
                        stages[0] = True
                        optimizer_no_incomplete = RMSDOptimizerAlt(copy.deepcopy(distance_dicts[gene_id]),
                                                                   info_dict[gene_id], False, True)
                        stages[1] = True
                        optimizer_no_non_coding = RMSDOptimizerAlt(copy.deepcopy(distance_dicts[gene_id]),
                                                                   info_dict[gene_id], True, False)
                        stages[2] = True
                        optimizer_both = RMSDOptimizerAlt(copy.deepcopy(distance_dicts[gene_id]),
                                                          info_dict[gene_id],
                                                          True,
                                                          True)
                        stages[3] = True
                        output_list.append([gene_id,
                                            str(round(optimizer.calc_rmsd(), 2)),
                                            str(round(optimizer_no_incomplete.calc_rmsd(), 2)),
                                            str(round(optimizer_no_non_coding.calc_rmsd(), 2)),
                                            str(round(optimizer_both.calc_rmsd(), 2))])
                    except ValueError:
                        print(gene_id, stages)
            with open(argument_dict["outfile"], "w") as f:
                f.write("\n".join(["\t".join(entry) for entry in output_list]))
        print("\n".join(["\t".join(entry) for entry in output_list]))
        with open(argument_dict["outfile"], "w") as f:
            f.write("\n".join(["\t".join(entry) for entry in output_list]))
    else:
        for i, fas_file in enumerate(fas_file_list):
            print(str(i) + "/" + str(total), fas_file)
            with open(os.path.join(argument_dict["library"], "fas_data", "fas_scores", fas_file), "r") as f:
                distance_dicts: Dict[str, Any] = json.load(f)
            for gene_id in distance_dicts.keys():
                if gene_id not in already_calced:
                    stages: List[bool] = [False, False, False, False]
                    try:
                        optimizer = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                  info_dict[gene_id],
                                                  False,
                                                  False)
                        stages[0] = True
                        optimizer_no_incomplete = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                                info_dict[gene_id], False, True)
                        stages[1] = True
                        optimizer_no_non_coding = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                                info_dict[gene_id], True, False)
                        stages[2] = True
                        optimizer_both = RMSDOptimizer(copy.deepcopy(distance_dicts[gene_id]),
                                                       info_dict[gene_id],
                                                       True,
                                                       True)
                        stages[3] = True
                        output_list.append([gene_id,
                                            str(round(optimizer.calc_rmsd(), 2)),
                                            str(round(optimizer_no_incomplete.calc_rmsd(), 2)),
                                            str(round(optimizer_no_non_coding.calc_rmsd(), 2)),
                                            str(round(optimizer_both.calc_rmsd(), 2))])
                    except ValueError:
                        print(gene_id, stages)
            with open(argument_dict["outfile"], "w") as f:
                f.write("\n".join(["\t".join(entry) for entry in output_list]))
        print("\n".join(["\t".join(entry) for entry in output_list]))
        with open(argument_dict["outfile"], "w") as f:
            f.write("\n".join(["\t".join(entry) for entry in output_list]))


if __name__ == "__main__":
    main()
