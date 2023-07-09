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


class BoxPlotTriple:

    def __init__(self, path_list: List[str],
                 descriptors: List[str],
                 values: List[str]):
        self.descriptors: List[str] = descriptors
        self.values: List[str] = values
        self.path_list: List[str] = path_list

    def plot_max_distribution(self):
        pass



def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-m",
                        "--maximum",
                        type=str,
                        action="store",
                        help="Filepath to the tsvs containing the maximum RMSD for each gene.")
    parser.add_argument("-d",
                        "--descriptors",
                        type=str,
                        action="store",
                        help="Descriptors for the given Max RMSD tsvs.")
    parser.add_argument("-v",
                        "--values",
                        type=str,
                        action="store",
                        help="Which values shall be contained. [all], [no_incomplete], [no_non_coding], [no_both]")
    argument_dict: Dict[str, Any] = vars(parser.parse_args())

    max_visualizer = BoxPlotTriple(argument_dict["maximum"],
                                   argument_dict["descriptors"],
                                   argument_dict["values"])
    max_visualizer.plot_max_distribution()


if __name__ == "__main__":
    main()
