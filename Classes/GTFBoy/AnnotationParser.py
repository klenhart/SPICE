#!/bin/env python
import os

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


class AnnotationParser:

    def __init__(self, annotation_path_list: List[str], expression_path_list: List[str],
                 threshold: float, name: str):
        self.name: str = name
        self.annotation_path_list: List[str] = annotation_path_list
        self.expression_path_list: List[str] = expression_path_list
        self.threshold = threshold
        if self.threshold is None:
            self.threshold: float = 1.0
        self.transcript_dict: Dict[str, Dict[str, Any]] = dict()

        self.new_id_map: Dict[str, Any] = dict() # Maps from new IDs to their genes and sequences.
        self.old_id_map: Dict[str, Any] = dict() # Maps from old ids to their new IDs.
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
                if gene_id not in self.old_id_map.keys():
                    self.old_id_map[gene_id] = dict()
                gene_dict = current_anno_dict[gene_id]

                for transcript_id in gene_dict.keys():
                    current_transcript_dict = gene_dict[transcript_id]

                    transcript_coord_id: str = AnnotationParser.make_coord_string(current_transcript_dict)
                    hash_coord_id: str = md5_hash(transcript_coord_id, 15)

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
                    self.old_id_map[gene_id][transcript_id] = hash_coord_id

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
                            if line_dict["transcript_id"] in self.old_id_map[gene_id].keys():
                                hash_id = self.old_id_map[gene_id][line_dict["transcript_id"]]
                                new_fpkm: float = float(line_dict["FPKM"])
                                if self.transcript_dict[gene_id][hash_id]["transcript"]["FPKM"] < new_fpkm:
                                    self.transcript_dict[gene_id][hash_id]["transcript"]["FPKM"] = new_fpkm
            print("Deleting transcripts below FPKM threshold.")
            gene_delete_list: List[str] = list()
            for gene_id in self.transcript_dict.keys():
                delete_list: List[str] = list()
                for transcript_id in self.transcript_dict[gene_id].keys():
                    if self.transcript_dict[gene_id][transcript_id]["transcript"]["FPKM"] < self.threshold:
                        delete_list.append(transcript_id)
                for transcript_id in delete_list:
                    self.remove_transcript(gene_id, transcript_id)
                if len(self.transcript_dict[gene_id].keys()) == 0:
                    gene_delete_list.append(gene_id)
            for gene_id in gene_delete_list:
                del self.transcript_dict[gene_id]
                del self.old_id_map[gene_id]
            print("Done.")
        print("Generating transcript id to gene id map.")
        self.__generate_id_map__()
        print("Updating stats of the merged annoation.")
        self.__update__()
        print("Done with update.")

    def remove_transcript(self, gene_id, transcript_id):
        synonym_list: List[str, str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["synonyms"]
        for synonym in synonym_list:
            del self.old_id_map[gene_id][synonym]
        del self.transcript_dict[gene_id][transcript_id]

    def __generate_id_map__(self):
        for gene_id in self.transcript_dict.keys():
            for transcript_id in self.transcript_dict[gene_id].keys():
                strand: str = self.transcript_dict[gene_id][transcript_id]["transcript"]["strand"]
                self.new_id_map[transcript_id] = dict()
                self.new_id_map[transcript_id]["gene_id"] = gene_id
                self.new_id_map[transcript_id]["strand"] = strand
                self.new_id_map[transcript_id]["peptides"] = list()

    def __iter__(self):
        yield "# Spice Annotation Parser Collection of Novel transcripts\n"
        yield "# " + str(self.novel_transcript_count) + " new transcripts\n"
        yield "# across " + str(self.gene_count) + " genes.\n"
        total = len(self.transcript_dict.keys())
        for i, gene_id in enumerate(self.transcript_dict.keys()):
            if i % 1000 == 0:
                print("Gene:", str(i) + "/" + str(total))
            for transcript_id in self.transcript_dict[gene_id].keys():
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
        with open(os.path.join(out_path, self.name + ".gtf"), "w") as f:
            f.write("")
        with open(os.path.join(out_path, self.name + ".gtf"), "a") as f:
            for line in self:
                f.write(line)

    def save_json(self, out_path: str):
        with open(os.path.join(out_path, self.name + "_idMap.json"), "w") as f:
            json.dump(self.old_id_map, f, indent=4)
        with open(os.path.join(out_path, self.name + ".json"), "w") as f:
            json.dump(self.new_id_map, f, indent=4)

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


def md5_hash(data, length: int = 32):
    hasher = hashlib.md5()
    hasher.update(data.encode("utf-8"))
    hash_string = hasher.hexdigest()
    return hash_string[:length]


def main():
    print("Hello World!")


if __name__ == "__main__":
    main()
