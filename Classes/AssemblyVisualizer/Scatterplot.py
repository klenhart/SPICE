#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  Scatterplot is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Scatterplot is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import argparse
import os
from typing import Dict, List, Set
import numpy as np

from matplotlib import pyplot as plt


class Scatterplot:

    def __init__(self, path1: str, path2: str, outpath: str, label1: str, label2: str):
        self.median_rank_1: Dict[str, int] = Scatterplot.extract_median_ranks(path1)
        self.median_rank_2: Dict[str, int] = Scatterplot.extract_median_ranks(path2)

        self.data_1_median_ranks, self.data_2_median_ranks = self.__fuse_median_ranks__()

        diagonal_line = np.array([min(self.data_1_median_ranks), max(self.data_1_median_ranks)])

        plt.scatter(self.data_1_median_ranks, self.data_2_median_ranks, s=10)

        plt.plot(diagonal_line, diagonal_line, color='red', linestyle='--')

        plt.ylim(-5, 1000)
        plt.xlim(-5, 1000)

        plt.xlabel(label1)
        plt.ylabel(label2)
        plt.title(r"$RMSD_{EWFD}$ Ranking comparison:" + "\n post EWFD vs pre EWFD filtering")

        plt.savefig(outpath, format='svg')

    def __fuse_median_ranks__(self):
        gene_set_1: Set[str] = set(self.median_rank_1.keys())
        gene_set_2: Set[str] = set(self.median_rank_2.keys())
        gene_set_1.update(gene_set_2)
        data_1_median_ranks = list()
        data_2_median_ranks = list()
        for gene_id in gene_set_1:
            if gene_id in self.median_rank_1.keys() and gene_id in self.median_rank_2.keys():
                data_1_median_ranks.append(self.median_rank_1[gene_id])
                data_2_median_ranks.append(self.median_rank_2[gene_id])
        return data_1_median_ranks, data_2_median_ranks

    @staticmethod
    def extract_median_ranks(path):
        filenames = os.listdir(path)
        median_rank_dict: Dict[str, List[int]] = dict()
        for filename in filenames:
            with open(os.path.join(path, filename), "r") as f:
                lines = f.read().split("\n")
                for i, line in enumerate(lines[1:]):
                    gene_id = line.split(",")[0]
                    if gene_id not in median_rank_dict.keys():
                        median_rank_dict[gene_id] = list()
                    median_rank_dict[gene_id].append(i)
        median_rank: Dict[str, int] = dict()
        for key in median_rank_dict.keys():
            median_rank_dict[key].sort()
            median_rank[key] = median_rank_dict[key][int(len(median_rank_dict[key]) / 2)]
        return median_rank


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        action="store",
                        nargs=2,
                        help="Paths to result folders.")
    parser.add_argument("-l",
                        "--label",
                        type=str,
                        action="store",
                        nargs=2,
                        help="Both Axis label for scatterplot.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        action="store",
                        help="Path to output file.")

    argument_dict: Dict[str, str] = vars(parser.parse_args())

    scatterplot: Scatterplot = Scatterplot(argument_dict["input"][0],
                                           argument_dict["input"][1],
                                           argument_dict["output"],
                                           argument_dict["label"][0],
                                           argument_dict["label"][1])


if __name__ == "__main__":
    main()
