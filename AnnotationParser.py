#!/bin/env python
import sys

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  AnnotationParser is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AnnotationParser is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.GTFBoy.GTFBoy import GTFBoy
from typing import List, Dict, Any

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse


class AnnotationParser:

    def __init__(self, annotation_path_list: List[str]):
        self.header: str = "# " + "Merged annotation of novel transcripts of protein coding genes."
        self.annotation_path_list: List[str] = annotation_path_list
        self.transcript_dict: Dict[str, Dict[str, Any]] = dict()
        self.novel_transcript_count: int = 0
        self.gene_count: int = 0

    def parse_annotations(self):
        total: int = len(self.annotation_path_list)
        for i, annotation_path in enumerate(self.annotation_path_list):
            print(str(i+1) + "/" + str(total), "Parsing ", annotation_path)
            gtf_iterator: GTFBoy = GTFBoy(annotation_path)
            for line in gtf_iterator:
                if not line.startswith("#"):
                    line_dict: Dict[str, str] = GTFBoy.build_dict(line.split("\t"))
                    line_dict["gene_id"] = line_dict["gene_id"].split(".")[0]
                    if line_dict["feature"] == "exon":
                        line_dict["exon_id"] = line_dict["exon_id"].split(".")[0]
                    if AnnotationParser.check_if_candidate(line_dict):
                        gene_id: str = line_dict["gene_id"]
                        transcript_id: str = line_dict["transcript_id"]
                        feature: str = line_dict["feature"]
                        if gene_id not in self.transcript_dict.keys():
                            self.transcript_dict[gene_id] = dict()
                        if transcript_id not in self.transcript_dict[gene_id].keys():
                            self.transcript_dict[gene_id][transcript_id] = dict()

                        gtf_line = GTFBoy.line_dict_to_gtf_line(line_dict)
                        if feature == "exon":
                            exon_id: str = line_dict["exon_id"]
                            self.transcript_dict[gene_id][transcript_id][exon_id] = gtf_line
                        else:
                            self.transcript_dict[gene_id][transcript_id]["GTF_line"] = gtf_line
        print("Done with parsing.")
        print("Updating stats of the merged annoation.")
        self.__update__()
        print("Done with update.")

    def __iter__(self):
        yield "# Spice Annotation Parser Collection of Novel transcripts\n"
        yield "# " + str(self.novel_transcript_count) + " new transcripts\n"
        yield "# across " + str(self.gene_count) + " genes.\n"
        total = len(self.transcript_dict.keys())
        for i, gene_id in enumerate(self.transcript_dict.keys()):
            if i % 1000 == 0:
                print("Gene:", str(i) + "/" + str(total))
            for transcript_id in self.transcript_dict[gene_id].keys():
                yield self.transcript_dict[gene_id][transcript_id]["GTF_line"]  + "\n"
                for exon_id in self.transcript_dict[gene_id][transcript_id].keys():
                    if exon_id != "GTF_line":
                        yield self.transcript_dict[gene_id][transcript_id][exon_id] + "\n"

    def __update__(self):
        self.gene_count = len(self.transcript_dict.keys())
        self.novel_transcript_count = 0
        for gene_id in self.transcript_dict.keys():
            self.novel_transcript_count += len(self.transcript_dict[gene_id].keys())

    def save(self, out_path: str):
        with open(out_path, "w") as f:
            f.write("")
        with open(out_path, "a") as f:
            for line in self:
                f.write(line)

    @staticmethod
    def check_if_candidate(line_dict: Dict[str, str]) -> bool:
        if line_dict["feature"] not in ["exon", "transcript"]:
            return False
        elif line_dict["gene_status"] == "NOVEL":
            return False
        elif "gene_type" not in line_dict.keys():
            print(line_dict)
            sys.exit()
        elif line_dict["gene_type"] != "protein_coding":
            return False
        elif line_dict["transcript_status"] == "NOVEL":
            return True
        else:
            return False


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--input", "--out_path"],
                                                   [str, str],
                                                   ["store", "store"],
                                                   [1, 1],
                                                   ["""Path to a text file containing the paths to all
                                                   gtfs that shall be merged. One path per line.""",
                                                    "Path to the output file.gtf."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, str] = argument_parser.get_args()
    argument_dict["input"] = argument_dict["input"][0]
    argument_dict["out_path"] = argument_dict["out_path"][0]

    with open(argument_dict["input"], "r") as f:
        anno_list: List[str] = f.read().split("\n")

    annotation_parser: AnnotationParser = AnnotationParser(anno_list)
    annotation_parser.parse_annotations()
    annotation_parser.save(argument_dict["out_path"])


if __name__ == "__main__":
    main()
