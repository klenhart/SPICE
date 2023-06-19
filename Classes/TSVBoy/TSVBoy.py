#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  TSVBoy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TSVBoy is distributed in the hope that it will be useful,
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


class DiamondTSVBoy:

    field_names: List[str] = ["query", "target", "seq_identity", "length", "mismatches", "gap_openings",
                              "query_start", "query_end", "target_start", "target_end", "e-value", "bit_score"]

    def __init__(self, tsv_path: str, group_key: str = "query", key_split: str = ".",
                 outer_split_index: int = 0, inner_split_index: int = 1, best_hit_flag: bool = True):
        self.tsv_path = tsv_path
        self.tsv_dict: Dict[str, Dict[str, List[Dict[str, Any]]]] = dict()
        self.group_key: str = group_key
        self.key_split: str = key_split
        self.out_split_index: int = outer_split_index
        self.in_split_index: int = inner_split_index
        self.best_hit_flag: bool = best_hit_flag

    def __iter__(self):
        with open(self.tsv_path, "r") as f:
            for line in f:
                yield line

    def get_dict(self) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
        return self.tsv_dict

    def parse_tsv(self):
        for line in self:
            if line == "":
                continue
            line = line.replace("\n", "")
            entry_list = line.split("\t")
            line_dict: Dict[str, Any] = dict(zip(DiamondTSVBoy.field_names, entry_list))
            outer_group_value: str = line_dict[self.group_key].split(self.key_split)[self.out_split_index]
            inner_group_value: str = line_dict[self.group_key].split(self.key_split)[self.in_split_index]


            line_dict["seq_identity"] = float(line_dict["seq_identity"])
            line_dict["length"] = int(line_dict["length"])
            line_dict["mismatches"] = int(line_dict["mismatches"])
            line_dict["gap_openings"] = int(line_dict["gap_openings"])
            line_dict["query_start"] = int(line_dict["query_start"])
            line_dict["query_end"] = int(line_dict["query_end"])
            line_dict["target_start"] = int(line_dict["target_start"])
            line_dict["target_end"] = int(line_dict["target_end"])
            line_dict["e-value"] = float(line_dict["e-value"])
            line_dict["bit_score"] = float(line_dict["bit_score"])

            if self.best_hit_flag:
                if outer_group_value in self.tsv_dict.keys():
                    if inner_group_value in self.tsv_dict[outer_group_value].keys():
                        bit_score: float = self.tsv_dict[outer_group_value][inner_group_value][0]["bit_score"]
                        if line_dict["bit_score"] > bit_score:
                            self.tsv_dict[outer_group_value][inner_group_value][0] = line_dict
                        continue
                    else:
                        self.tsv_dict[outer_group_value][inner_group_value] = list()
                else:
                    self.tsv_dict[outer_group_value] = dict()
                    self.tsv_dict[outer_group_value][inner_group_value] = list()

            else:
                if outer_group_value not in self.tsv_dict.keys():
                    self.tsv_dict[outer_group_value] = dict()
                if inner_group_value not in self.tsv_dict[outer_group_value].keys():
                    self.tsv_dict[outer_group_value][inner_group_value] = list()

            self.tsv_dict[outer_group_value][inner_group_value].append(line_dict)

    def save(self, out_path: str):
        with open(out_path, "w") as f:
            json.dump(self.tsv_dict, f, indent=4)


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        action="store",
                        help="Path to a tsv file.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        action="store",
                        help="Path to json output file.")
    parser.add_argument("-m",
                        "--mode",
                        type=str,
                        action="store",
                        help="Only 'diamond' for now.")

    argument_dict: Dict[str, str] = vars(parser.parse_args())

    if argument_dict["mode"] == "diamond":
        tsv_iterator: DiamondTSVBoy = DiamondTSVBoy(argument_dict["input"])
        tsv_iterator.parse_tsv()
        tsv_iterator.save(argument_dict["output"])


if __name__ == "__main__":
    main()
