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
    def line_dict_to_gtf_line(line_dict: Dict[str, str]) -> str:
        gtf_line = ""
        gtf_line += line_dict["seqname"] + "\t"
        gtf_line += "SPICE\t"
        gtf_line += line_dict["feature"] + "\t"
        gtf_line += line_dict["start"] + "\t"
        gtf_line += line_dict["end"] + "\t"
        gtf_line += line_dict["score"] + "\t"
        gtf_line += line_dict["strand"] + "\t"
        gtf_line += line_dict["frame"] + "\t"
        gtf_line += "gene_id \"" + line_dict["gene_id"] + "\"; "
        gtf_line += "transcript_id \"" + line_dict["transcript_id"] + "\"; "
        gtf_line += "gene_biotype \"protein_coding\"; "
        gtf_line += "transcript_biotype \"protein_coding\"; "
        gtf_line += "gene_status \"NOVEL\"; "
        if line_dict["feature"] == "exon":
            gtf_line += "exon_id \"" + line_dict["exon_id"] + "\"; "
        return gtf_line

    @staticmethod
    def has_attribute_value(attribute_key: str, attribute_value: str, attribute_entry: str):
        attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(attribute_entry)
        # if attribute_key not in attribute_dict.keys():
        #     print("TEST MARKER")
        #     print(attribute_dict)
        #     return False
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


def main() -> None:
    pass


if __name__ == "__main__":
    main()
