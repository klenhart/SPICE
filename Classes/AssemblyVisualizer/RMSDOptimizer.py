#!/bin/env python
import argparse
import json
import os

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ResultBuddy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ResultBuddy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import numpy as np
from typing import Dict, List, Any, Tuple
from scipy.optimize import minimize
import math
import random


class RMSDOptimizer:

    def __init__(self,
                 distance_dict: Dict[str, Dict[str, float]],
                 info_dict: Dict[str, Any],
                 protein_coding_flag: bool = True,
                 complete_flag: bool = True):
        distance_dict = RMSDOptimizer.cull_distance_dict(distance_dict, info_dict, protein_coding_flag, complete_flag)
        self.matrix = RMSDOptimizer.distance_dict_to_matrix(distance_dict)
        self.length: int = len(self.matrix)
        self.bounds = [(0, 1)] * (2 * self.length)
        self.result = self.__optimize__()
        count: int = 0
        while count != 30:
            new_result = self.__optimize__()
            if -self.objective_function(new_result.x) > -self.objective_function(self.result.x):
                self.result = new_result
            else:
                count = count+1

    def objective_function(self, combined_vector):
        v = combined_vector[:self.length]
        w = combined_vector[self.length:]
        diff = np.dot(self.matrix, v) - np.dot(self.matrix, w)
        return -np.sqrt(np.sum(diff ** 2) / self.length)

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
        remain = 10

        for i in range(n):
            if remain == 0:
                break
            elif i == n-1:
                vector[i] = round(remain*0.1, 1)
            else:
                x = random.randint(1, remain)
                vector[i] = round(x*0.1, 1)
                remain = remain - x
        return vector

    @staticmethod
    def cull_distance_dict(distance_dict: Dict[str, Dict[str, float]],
                           info_dict: Dict[str, Any],
                           protein_coding_flag,
                           complete_flag):
        if protein_coding_flag or complete_flag:
            delete_set = set()
            if protein_coding_flag:
                for key in distance_dict.keys():
                    if RMSDOptimizer.is_non_coding(key):
                        delete_set.add(key)
            if complete_flag:
                for key in distance_dict.keys():
                    if RMSDOptimizer.is_incomplete(key, info_dict):
                        delete_set.add(key)
            for entry in delete_set:
                del distance_dict[entry]
            for key in distance_dict.keys():
                for entry in delete_set:
                    del distance_dict[key][entry]
        return distance_dict

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
    argument_dict: Dict[str, str] = vars(parser.parse_args())

    output_list: List[List[str]] = [["geneid", "all", "no_incomplete", "no_non_coding", "no_both"]]

    with open(os.path.join(argument_dict["library"], "transcript_data", "transcript_info.json"), "r") as f:
        info_dict = json.load(f)

    with open(os.path.join(argument_dict["library"], "fas_data", "fas_index.json"), "r") as f:
        fas_file_list = list(set(json.load(f).values()))

    for fas_file in fas_file_list:
        with open(os.path.join(argument_dict["library"], "fas_data", "fas_scores", fas_file), "r") as f:
            distance_dicts: Dict[str, Any] = json.load(f)




if __name__ == "__main__":
    main()
