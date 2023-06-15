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
import json
import hashlib
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse


class AnnotationParser:

    def __init__(self, annotation_path_list: List[str], expression_path_list: List[str], threshold: float):
        self.annotation_path_list: List[str] = annotation_path_list
        self.expression_path_list: List[str] = expression_path_list
        self.threshold = threshold
        if self.threshold is None:
            self.threshold: float = 1.0
        self.transcript_dict: Dict[str, Dict[str, Any]] = dict()
        self.novel_transcript_count: int = 0
        self.gene_count: int = 0

    def parse_annotations(self):
        total: int = len(self.annotation_path_list)
        for i, annotation_path in enumerate(self.annotation_path_list):
            if annotation_path == "":
                print("Encountered empty entry. This is the whole list of gtf paths:\n",
                      self.annotation_path_list,
                      "\nThe current index of annotation parsing is", i)
                continue
            print(str(i+1) + "/" + str(total), "Parsing ", annotation_path)
            gtf_iterator: GTFBoy = GTFBoy(annotation_path)
            current_anno_dict: Dict[str, Dict[str, Any]] = dict()
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
                        if gene_id not in current_anno_dict.keys():
                            current_anno_dict[gene_id] = dict()
                        if transcript_id not in current_anno_dict[gene_id].keys():
                            current_anno_dict[gene_id][transcript_id] = dict()

                        if feature == "exon":
                            exon_id: str = line_dict["exon_id"]
                            current_anno_dict[gene_id][transcript_id][exon_id] = line_dict
                        else:
                            current_anno_dict[gene_id][transcript_id]["transcript"] = line_dict

            print("\tGenerating unique genome coordinate strings.")
            for gene_id in current_anno_dict.keys():
                if gene_id not in self.transcript_dict.keys():
                    self.transcript_dict[gene_id] = dict()
                    self.transcript_dict[gene_id]["id_map"] = dict()
                gene_dict = current_anno_dict[gene_id]

                for transcript_id in gene_dict.keys():
                    current_transcript_dict = gene_dict[transcript_id]

                    transcript_coord_id: str = AnnotationParser.make_coord_string(current_transcript_dict)
                    hash_coord_id: str = md5_hash(transcript_coord_id)

                    synonym_field: List[str] = [current_transcript_dict["transcript"]["transcript_id"]]
                    current_transcript_dict["transcript"]["synonyms"] = synonym_field
                    current_transcript_dict["transcript"]["coord_id"] = transcript_coord_id
                    current_transcript_dict["transcript"]["FPKM"] = 0.0

                    for key in current_transcript_dict.keys():
                        current_transcript_dict[key]["transcript_id"] = hash_coord_id

                    if hash_coord_id not in self.transcript_dict[gene_id].keys():
                        self.transcript_dict[gene_id][hash_coord_id] = current_transcript_dict
                    else:
                        if transcript_id not in self.transcript_dict[gene_id][hash_coord_id]["transcript"]["synonyms"]:
                            self.transcript_dict[gene_id][hash_coord_id]["transcript"]["synonyms"].append(transcript_id)
                    self.transcript_dict[gene_id]["id_map"][transcript_id] = hash_coord_id

        if len(self.expression_path_list) > 0:
            print("Filtering transcripts. Expression (FPKM) must be greater than", str(self.threshold) + ".")
            print("Assigning maximum FPKMs.")
            for i, path in enumerate(self.expression_path_list):
                if path == "":
                    print("Encountered empty expression path entry. This is the whole list of gtf paths:\n",
                          self.annotation_path_list,
                          "\nThe current index of annotation parsing is", i)
                    continue
                gtf_iterator: GTFBoy = GTFBoy(path)
                for line in gtf_iterator:
                    if line.startswith("#") or line == "":
                        continue
                    line_dict: Dict[str, str] = GTFBoy.build_dict(line.split("\t"))
                    if line_dict["feature"] == "transcript":
                        gene_id = line_dict["gene_id"].split(".")[0]
                        if gene_id in self.transcript_dict.keys():
                            if line_dict["transcript_id"] in self.transcript_dict[gene_id]["id_map"].keys():
                                hash_id = self.transcript_dict[gene_id]["id_map"][line_dict["transcript_id"]]
                                new_fpkm: float = float(line_dict["FPKM"])
                                if self.transcript_dict[gene_id][hash_id]["transcript"]["FPKM"] < new_fpkm:
                                    self.transcript_dict[gene_id][hash_id]["transcript"]["FPKM"] = new_fpkm
            print("Deleting transcripts below FPKM threshold.")
            gene_delete_list: List[str] = list()
            for gene_id in self.transcript_dict.keys():
                delete_list: List[str] = list()
                for transcript_id in self.transcript_dict[gene_id].keys():
                    if transcript_id != "id_map":
                        if self.transcript_dict[gene_id][transcript_id]["transcript"]["FPKM"] < self.threshold:
                            delete_list.append(transcript_id)
                for transcript_id in delete_list:
                    self.remove_transcript(gene_id, transcript_id)
                if len(self.transcript_dict[gene_id].keys()) == 1:
                    gene_delete_list.append(gene_id)
            for gene_id in gene_delete_list:
                del self.transcript_dict[gene_id]
            print("Done.")
        print("Updating stats of the merged annoation.")
        self.__update__()
        print("Done with update.")

    def remove_transcript(self, gene_id, transcript_id):
        synonym_list: List[str, str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["synonyms"]
        for synonym in synonym_list:
            del self.transcript_dict[gene_id]["id_map"][synonym]
        del self.transcript_dict[gene_id][transcript_id]

    def __iter__(self):
        yield "# Spice Annotation Parser Collection of Novel transcripts\n"
        yield "# " + str(self.novel_transcript_count) + " new transcripts\n"
        yield "# across " + str(self.gene_count) + " genes.\n"
        total = len(self.transcript_dict.keys())
        for i, gene_id in enumerate(self.transcript_dict.keys()):
            if i % 1000 == 0:
                print("Gene:", str(i) + "/" + str(total))
            for transcript_id in self.transcript_dict[gene_id].keys():
                if transcript_id == "id_map":
                    continue
                yield GTFBoy.line_dict_to_gtf_line(self.transcript_dict[gene_id][transcript_id]["transcript"]) + "\n"
                for exon_id in self.transcript_dict[gene_id][transcript_id].keys():
                    if exon_id != "transcript":
                        yield GTFBoy.line_dict_to_gtf_line(self.transcript_dict[gene_id][transcript_id][exon_id]) + "\n"

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

    def save_json(self, out_path: str):
        with open(out_path, "w") as f:
            json.dump(self.transcript_dict, f, indent=4)

    @staticmethod
    def check_if_candidate(line_dict: Dict[str, str]) -> bool:
        if line_dict["feature"] not in ["exon", "transcript"]:
            return False
        elif not line_dict["seqname"].startswith("chr"):
            return False
        elif line_dict["seqname"].startswith("chrM"):
            return False
        elif line_dict["seqname"].endswith("random"):
            return False
        elif line_dict["seqname"].startswith("chrU"):
            return False
        elif line_dict["gene_status"] == "NOVEL":
            return False
        elif line_dict["gene_type"] != "protein_coding":
            return False
        elif line_dict["transcript_status"] == "NOVEL":
            return True
        else:
            return False

    @staticmethod
    def make_coord_string(transcript_dict: Dict[str, Any]):
        chromosome: str = ""
        exon_list = list()
        for key in transcript_dict.keys():
            if key == "transcript":
                chromosome = transcript_dict[key]["seqname"]
            else:
                exon = transcript_dict[key]["start"] + "-" + transcript_dict[key]["end"]
                exon_list.append(exon)
        return "_".join(sorted(exon_list)) + "_" + chromosome


def md5_hash(data):
    hasher = hashlib.md5()
    hasher.update(data.encode("utf-8"))
    hash_string = hasher.hexdigest()
    return hash_string


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--input", "--out_path", "--json", "--expression", "--threshold"],
                                                   [str, str, str, str, float],
                                                   ["store", "store", "store", "store", "store"],
                                                   [1, 1, None, "?", "?"],
                                                   ["""Path to a text file containing the paths to all
                                                   gtfs that shall be merged. One path per line.""",
                                                    "Path to the output file.gtf.",
                                                    """If a json version of the annotation is required
                                                     a path needs to be added to this argument.""",
                                                    """Path to a text file containing the paths to gtf files
                                                     including transcript abundances. Only transcripts will be kept
                                                     which exceed the float given in the --threshold parameter.""",
                                                    "Expression threshold, which will be used for transcript curation."
                                                    ])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict["input"] = argument_dict["input"][0]
    argument_dict["out_path"] = argument_dict["out_path"][0]

    with open(argument_dict["input"], "r") as f:
        anno_list: List[str] = f.read().split("\n")

    if argument_dict["expression"] is None:
        expression_list: List[str] = list()
    else:
        with open(argument_dict["expression"], "r") as f:
            expression_list: List[str] = f.read().split("\n")

    annotation_parser: AnnotationParser = AnnotationParser(anno_list,
                                                           expression_list,
                                                           argument_dict["threshold"])
    annotation_parser.parse_annotations()
    annotation_parser.save(argument_dict["out_path"])

    if argument_dict["json"] is not None:
        annotation_parser.save_json(argument_dict["json"])


if __name__ == "__main__":
    main()
