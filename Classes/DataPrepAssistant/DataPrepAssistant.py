#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  DataPrepAssistant is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  DataPrepAssistant is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import argparse
import json
import os
from typing import Dict, List, Any


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        action="store",
                        help="Path to directory that shall be parsed.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        action="store",
                        help="Path to store the lists.")
    parser.add_argument("-p",
                        "--prefix",
                        type=str,
                        action="store",
                        help="Prefix of the folders that shall be searched.")

    argument_dict: Dict[str, str] = vars(parser.parse_args())

    input_path: str = argument_dict["input"]
    output_path: str = argument_dict["output"]
    prefix: str = argument_dict["prefix"]

    condition_dict: Dict[str, List[str]] = dict()
    expression_list: List[str] = list()
    annotation_list: List[str] = list()
    name_list: List[str] = list()

    candidates: List[str] = list()

    for entry in os.listdir(argument_dict["input"]):
        if os.path.isdir(os.path.join(input_path, entry)) and entry.startswith(prefix):
            candidates.append(os.path.join(input_path, entry))

    for entry in candidates:
        with open(os.path.join(entry, "info.json"), "r") as f:
            info_dict: Dict[str, Any] = json.load(f)
        print("-------->")
        print("Entry:", info_dict["description"])
        print(info_dict["biosample_summary"])
        print()
        decision_flag: bool = False
        while not decision_flag:
            key_input = input("Add? [y/n]")
            if key_input == "y":
                replicate_list: List[int] = list(info_dict["replicate_coverage_relation"].keys())
                for value in replicate_list:
                    expression_list.append(info_dict["replicate_coverage_relation"][value])
                for value in replicate_list:
                    annotation_list.append(os.path.join(entry, info_dict["replicate_anno_relation"][value] + ".gtf"))
                decision_flag = True
                name_flag = False
                replicate_names: List[str] = list()
                while not name_flag:
                    print("Get current names with ##")
                    name: str = input("Choose name:")
                    if name == "##":
                        print(name_list)
                    elif name == "":
                        continue
                    else:
                        if name + "_rep1" in replicate_names:
                            print("Name was already chosen. Retry.")
                            continue
                        name_flag = True
                        for value in replicate_list:
                            replicate_names.append(name + "_rep" + str(value))
                        name_list += replicate_names

                condition_flag = False
                while not condition_flag:
                    print("Get current conditions with ##")
                    condition: str = input("Choose condition:")
                    if condition == "##":
                        print(list(condition_dict.keys()))
                    elif condition == "":
                        continue
                    else:
                        if condition not in condition_dict.keys():
                            condition_dict[condition] = list()
                        condition_dict[condition] += replicate_names
                        done_flag: bool = False
                        while not done_flag:
                            key_done: str = input("Done? [y/n]")
                            if key_done == "y":
                                done_flag = True
                                condition_flag = True
                            elif key_done == "n":
                                done_flag = True
                            else:
                                continue

            elif key_input == "n":
                decision_flag = True
            else:
                continue
        print("<--------")

    names_file = "\n".join(name_list)
    coverage_file = "\n".join(expression_list)
    annotations_file = "\n".join(annotation_list)

    conditions_temp_list = list()
    replicates_temp_list = list()
    for condition in condition_dict.keys():
        conditions_temp_list.append(condition)
        replicates_temp_list.append(" ".join(condition_dict[condition]))

    conditions_file = "\n".join(conditions_temp_list)
    replicates_file = "\n".join(replicates_temp_list)

    names_path = os.path.join(output_path, "names.txt")
    coverage_path = os.path.join(output_path, "coverage_list.txt")
    anno_path = os.path.join(output_path, "annotation_list.txt")
    condition_path = os.path.join(output_path, "conditions.txt")
    replicates_path = os.path.join(output_path, "replicates.txt")

    with open(names_path, "w") as f:
        f.write(names_file)

    with open(coverage_path, "w") as f:
        f.write(coverage_file)

    with open(anno_path, "w") as f:
        f.write(annotations_file)

    with open(condition_path, "w") as f:
        f.write(conditions_file)

    with open(replicates_path, "w") as f:
        f.write(replicates_file)


if __name__ == "__main__":
    main()