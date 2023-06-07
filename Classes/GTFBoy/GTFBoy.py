#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  GTFBoy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GTFBoy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Iterator, Dict, List
import argparse
import os
from pathlib import Path

class GTFBoy:

    start_incomplete_tags: List[str] = ['mRNA_start_NF', 'cds_start_NF']
    end_incomplete_tags: List[str] = ['mRNA_end_NF', 'cds_end_NF']

    GTF_MASK: List[str] = ["seqname", "source", "feature",
                           "start", "end", "score",
                           "strand", "frame", "attribute"]

    def __init__(self, gtf_path: str):
        self.gtf_path: str = gtf_path
        with open(gtf_path) as f:
            self.total_lines = len(f.readlines())

    def __iter__(self) -> Iterator[str]:
        with open(self.gtf_path, "r") as f:
            for line in f:
                yield line

    def filter_gtf(self):
        pass

    def process_gtf(self):
        pass

    @staticmethod
    def build_attribute_dict(attribute_entry: str) -> Dict[str, str]:
        attribute_dict: Dict[str, str] = dict()
        attribute_list: List[str] = attribute_entry.split(";")
        tag_list: List[str] = list()
        for entry in attribute_list:
            pair: List[str] = entry.split(" ")
            if len(pair) == 1:
                continue
            elif len(pair) == 2:
                pass
            elif len(pair) == 3:
                pair = pair[1:]
            elif len(pair) > 3:
                pair = pair[1:3]

            pair[1] = pair[1].replace("\"", "")

            if pair[0] == "tag":
                tag_list.append(pair[1])
                continue
            elif pair[0] == "transcript_support_level":
                if pair[1][0] != "N":
                    attribute_dict["transcript_support_level"] = pair[1][0]
            else:
                attribute_dict[pair[0]] = pair[1]
        # Search for tags indicating a 3'/5' incomplete transcript.
        if any([tag in GTFBoy.start_incomplete_tags for tag in tag_list]):
            tag_list.append("start_incomplete")
            tag_list.append("incomplete")
        elif any([tag in GTFBoy.end_incomplete_tags for tag in tag_list]):
            tag_list.append("end_incomplete")
            tag_list.append("incomplete")
        else:
            tag_list.append("complete")
        tag_string = ";".join(tag_list)
        attribute_dict["tag"] = tag_string
        if "transcript_support_level" not in attribute_dict.keys():
            attribute_dict["transcript_support_level"] = "6"
        return attribute_dict

    @staticmethod
    def build_dict(gtf_split_line: List[str]) -> Dict[str, str]:
        full_dict: Dict[str, str] = dict()
        attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(gtf_split_line[8])
        for i, field_name in enumerate(GTFBoy.GTF_MASK[:-1]):
            full_dict[field_name] = gtf_split_line[i]
        full_dict.update(attribute_dict)
        return full_dict

    @staticmethod
    def has_attribute_value(attribute_key: str, attribute_value: str, attribute_entry: str):
        attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(attribute_entry)
        #if attribute_key not in attribute_dict.keys():
        #    print("TEST MARKER")
        #    print(attribute_dict)
        #    return False
        return attribute_dict[attribute_key] == attribute_value

    @staticmethod
    def has_values(inclusion_filter_dict: Dict[str, List[str]], gtf_split_line: List[str]) -> bool:
        flag: bool = True
        full_dict: Dict[str, str] = GTFBoy.build_dict(gtf_split_line)
        for key, values in inclusion_filter_dict.items():
            if key in full_dict.keys():
                if full_dict[key] not in values:
                    flag = False
                    break
        return flag


def discriminate_line_dict(line_dict: Dict[str, str]) -> bool:
    seqname_flag: bool = line_dict["seqname"].startswith("chr")
    gene_id_flag: bool = line_dict["gene_id"].startswith("ENS")
    return seqname_flag and gene_id_flag


def make_memory_dict_entry(line_dict: Dict[str, str]) -> str:
    return "#".join([line_dict["gene_id"].split(".")[0],
                     line_dict["seqname"],
                     line_dict["feature"],
                     line_dict["start"],
                     line_dict["end"]])


def main() -> None:
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="store", help="""File path to the txt containing all file
     paths to the GTFs that shall get merged.""")
    parser.add_argument("-o", "--out_path", action="store", help="""Path to the directory that shall
     contain the merged annotation file.""")

    args_dict: Dict[str, str] = vars(parser.parse_args())

    with open(args_dict["input"], "r") as f:
        annotation_files: List[str] = f.read().split("\n")

    Path(os.path.join(args_dict["out_path"], "temp")).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(args_dict["out_path"], "temp", "memory.temp")).touch()
    Path(os.path.join(args_dict["out_path"], "merged_annotation.gtf")).touch()

    memory_set: set = set()  # gene_id#seqname#feature#start#end

    for gtf_path in annotation_files:
        gtf_boy: GTFBoy = GTFBoy(gtf_path)
        for line in gtf_boy:
            if line.startswith("#"):
                pass
            else:
                line_dict: Dict[str, str] = GTFBoy.build_dict(line.split("\t"))
                if discriminate_line_dict(line_dict, memory_set):

        break

if __name__ == "__main__":
    main()
