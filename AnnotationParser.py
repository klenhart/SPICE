#!/bin/env python

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


class AnnotationParser:

    def __init__(self, annotation_path_list: List[str]):
        self.header: str = "# " + "Merged annotation of novel transcripts of protein coding genes."
        self.annotation_path_list: List[str] = annotation_path_list
        self.transcript_dict: Dict[str, Dict[str, Any]] = dict()

    def parse_annotations(self):
        for annotation_path in self.annotation_path_list:
            gtf_iterator: GTFBoy = GTFBoy(annotation_path)
            count = 0
            for line in gtf_iterator:
                if not line.startswith("#"):
                    line_dict: Dict[str, str] = GTFBoy.build_dict(line.split("\t"))
                    line_dict["gene_id"] = line_dict["gene_id"].split(".")[0]
                    if line_dict["feature"] == "exon":
                        line_dict["exon_id"] = line_dict["exon_id"].split(".")[0]
                    if AnnotationParser.check_if_candidate(line_dict):
                        count += 1
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
                if count >= 20:
                    count = 0
                    break
        print(self)

    def __str__(self) -> str:
        output_string: str = ""
        for gene_id in self.transcript_dict.keys():
            output_string += gene_id + "\n"
            for transcript_id in self.transcript_dict[gene_id].keys():
                output_string += "\t" + transcript_id + "\n"
                for exon_id in self.transcript_dict[gene_id][transcript_id].keys():
                    output_string += "\t\t" + exon_id + "\n"
                    output_string += "\t\t\t" + self.transcript_dict[gene_id][transcript_id][exon_id] + "\n"
        return output_string

    def save(self, out_path: str):
        pass

    @staticmethod
    def check_if_candidate(line_dict: Dict[str, str]) -> bool:
        if line_dict["feature"] not in ["exon", "transcript"]:
            return False
        elif line_dict["gene_type"] != "protein_coding":
            return False
        elif line_dict["gene_status"] == "NOVEL":
            return False
        elif line_dict["transcript_status"] == "NOVEL":
            return True
        else:
            return False




def main():
    anno_list: List[str] = ["C:/Users/chris/Desktop/ENCSR902GAF/ENCFF298VFJ.gtf",
                            "C:/Users/chris/Desktop/ENCSR902GAF/ENCFF509CYO.gtf",
                            "C:/Users/chris/Desktop/ENCSR902GAF/ENCFF596NHC.gtf"]
    annotation_parser: AnnotationParser = AnnotationParser(anno_list)
    annotation_parser.parse_annotations()

if __name__ == "__main__":
    main()
