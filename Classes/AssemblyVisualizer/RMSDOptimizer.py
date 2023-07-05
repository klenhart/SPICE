#!/bin/env python
import math
import random

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
from typing import Dict, List
from scipy.optimize import minimize


class RMSDOptimizer:

    def __init__(self, matrix: np.array):
        self.matrix: np.array = matrix
        self.length: int = len(matrix)
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
    def distance_dict_to_matrix(distance_dict: Dict[str, Dict[str, float]]) -> np.array:
        raw_matrix: List[List[float]] = list()
        for protein_id_outer in distance_dict.keys():
            if not RMSDOptimizer.is_non_coding(protein_id_outer):
                raw_matrix.append(list())
                for protein_id_inner in distance_dict[protein_id_outer].keys():
                    if not RMSDOptimizer.is_non_coding(protein_id_inner):
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

    def calc_rmsd(self):
        return -self.objective_function(self.result.x)


def main():

    distance_dict: Dict[str, Dict[str, float]] = {
        "ENSP00000411638": {
            "ENSP00000411638": 1.0,
            "ENSP00000442360": 0.9916,
            "ENSP00000320343": 0.7362
        },
        "ENSP00000442360": {
            "ENSP00000411638": 0.9916,
            "ENSP00000442360": 1.0,
            "ENSP00000320343": 0.7372
        },
        "ENSP00000320343": {
            "ENSP00000411638": 0.9816,
            "ENSP00000442360": 0.983,
            "ENSP00000320343": 1.0
        }
    }

    optimizer = RMSDOptimizer(RMSDOptimizer.distance_dict_to_matrix(distance_dict))
    print(optimizer.matrix)
    print(optimizer.result.x)
    print(optimizer.calc_rmsd())


if __name__ == "__main__":
    main()
