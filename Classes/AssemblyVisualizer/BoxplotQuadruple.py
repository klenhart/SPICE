#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  BoxplotTriple is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BoxplotTriple is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import List, Dict, Any
import argparse
import json
import math
import numpy as np
import os
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline


class BoxPlotQuadruple:

    def __init__(self, path_list: List[str],
                 descriptors: List[str],
                 values: List[str]):
        self.descriptors: List[str] = descriptors
        self.values: List[str] = values
        self.path_list: List[str] = path_list

    def plot_max_distribution(self):
        data = self.extract_data()

        fig, ax = plt.subplots()

        concat_data = list()

        for path in self.path_list:
            for key in self.values:
                concat_data.append(data[path][key])

        bp = ax.boxplot(concat_data, showfliers=False)

        ax.set_xticklabels(self.values * len(self.path_list))

        for i, box in enumerate(bp['boxes']):
            if i < 2:
                box.set(color='green')
            elif i < 4:
                box.set(color='red')
            elif i < 6:
                box.set(color='blue')
            elif i < 8:
                box.set(color='yellow')
        for i, whisker in enumerate(bp['whiskers']):
            if i < 2:
                whisker.set(color='green')
            elif i < 4:
                whisker.set(color='red')
            elif i < 6:
                whisker.set(color='blue')
            elif i < 8:
                whisker.set(color='yellow')

        for i, cap in enumerate(bp['caps']):
            if i < 2:
                cap.set(color='green')
            elif i < 4:
                cap.set(color='red')
            elif i < 6:
                cap.set(color='blue')
            elif i < 8:
                cap.set(color='yellow')

        ax.set_title('Boxplot of Four Lists')

        plt.show()

    def extract_data(self):
        data: Dict[str, Dict[str, List[float]]] = dict()
        for path in self.path_list:
            data[path] = dict()
            for key in self.values:
                data[path][key] = list()
        for path in self.path_list:
            with open(path, "r") as f:
                lines = f.read().split("\n")
                col_names = lines[0].split("\t")
            for line in lines[1:]:
                entry_list = line.split("\t")
                for key in self.values:
                    i = col_names.index(key)
                    if float(entry_list[i]) <= 1.0:
                        data[path][key].append(float(entry_list[i]))
        return data




def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-m",
                        "--maximum",
                        type=str,
                        action="store",
                        nargs="+",
                        help="Filepath to the tsvs containing the maximum RMSD for each gene.")
    parser.add_argument("-d",
                        "--descriptors",
                        type=str,
                        action="store",
                        nargs="+",
                        help="Descriptors for the given Max RMSD tsvs.")
    parser.add_argument("-v",
                        "--values",
                        type=str,
                        action="store",
                        nargs="+",
                        help="Which values shall be contained. [all], [no_incomplete], [no_non_coding], [no_both]")
    argument_dict: Dict[str, Any] = vars(parser.parse_args())

    max_visualizer = BoxPlotQuadruple(argument_dict["maximum"],
                                      argument_dict["descriptors"],
                                      argument_dict["values"])
    max_visualizer.plot_max_distribution()


if __name__ == "__main__":
    main()
