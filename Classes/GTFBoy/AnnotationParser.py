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
from typing import List, Dict, Any, Set
import json
import hashlib
import os


class AnnotationParser:

    def __init__(self, annotation_path_list: List[str], protein_coding_gene_set: Set[str],
                 threshold: float, name: str, species_name: str, check_threshold_flag: bool):
        self.name: str = name
        self.annotation_path_list: List[str] = annotation_path_list
        self.protein_coding_gene_set: Set[str] = protein_coding_gene_set
        self.check_threshold_flag: bool = check_threshold_flag
        self.threshold = threshold
        self.transcript_dict: Dict[str, Dict[str, Any]] = dict()
        if len(species_name) > 3:
            self.species_prefix: str = species_name[:3].upper()
        elif len(species_name) == 0:
            self.species_prefix: str = ""
        else:
            self.species_prefix: str = species_name.upper()

        self.new_id_map: Dict[str, Any] = dict()  # Maps from new IDs to their genes and sequences.
        self.old_id_map: Dict[str, Any] = dict()  # Maps from old ids to their new IDs.
        self.novel_transcript_count: int = 0

        self.gene_count: int = 0

    def parse_annotations(self):
        transcripts_threshold_dropped: Set[str] = set()
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
                    line_dict["transcript_id"] = line_dict["transcript_id"].split(".")[0]
                    if line_dict["gene_id"] not in self.protein_coding_gene_set:
                        continue
                    is_transcript: bool = line_dict["feature"] == "transcript"
                    if is_transcript and self.check_threshold_flag and float(line_dict["FPKM"]) < self.threshold:
                        transcripts_threshold_dropped.add(line_dict["transcript_id"])
                        continue
                    if line_dict["feature"] == "exon":
                        if line_dict["transcript_id"] in transcripts_threshold_dropped:
                            continue
                        if "exon_number" in line_dict.keys() and "exon_id" not in line_dict.keys():
                            line_dict["exon_id"] = line_dict["exon_number"]
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
                    hash_coord_id: str = self.species_prefix + md5_hash(transcript_coord_id, 15)

                    synonym_field: List[str] = [current_transcript_dict["transcript"]["transcript_id"]]
                    current_transcript_dict["transcript"]["synonyms"] = synonym_field
                    current_transcript_dict["transcript"]["coord_id"] = transcript_coord_id

                    for key in current_transcript_dict.keys():
                        current_transcript_dict[key]["transcript_id"] = hash_coord_id

                    if hash_coord_id not in self.transcript_dict[gene_id].keys():
                        self.transcript_dict[gene_id][hash_coord_id] = current_transcript_dict
                    else:
                        if transcript_id not in self.transcript_dict[gene_id][hash_coord_id]["transcript"]["synonyms"]:
                            self.transcript_dict[gene_id][hash_coord_id]["transcript"]["synonyms"].append(transcript_id)
                    self.old_id_map[gene_id][transcript_id] = hash_coord_id

        print("Generating transcript id to gene id map.")
        self.__generate_id_map__()
        print("Updating stats of the merged annotation.")
        self.__update__()
        print("Done with update.")

    def remove_transcript(self, gene_id, transcript_id):
        del self.transcript_dict[gene_id][transcript_id]

    def __generate_id_map__(self):
        for gene_id in self.transcript_dict.keys():
            for transcript_id in self.transcript_dict[gene_id].keys():
                strand: str = self.transcript_dict[gene_id][transcript_id]["transcript"]["strand"]
                synonyms: List[str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["synonyms"]
                start: List[str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["start"]
                end: List[str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["end"]
                chromosome: List[str] = self.transcript_dict[gene_id][transcript_id]["transcript"]["seqname"]
                if transcript_id in self.new_id_map.keys():
                    print("ERROR:", transcript_id, "already given out to", self.new_id_map[transcript_id])
                    print("Wanted to write", self.transcript_dict[gene_id][transcript_id]["transcript"], "to key.")
                    raise KeyError
                else:
                    self.new_id_map[transcript_id] = dict()
                    self.new_id_map[transcript_id]["gene_id"] = gene_id
                    self.new_id_map[transcript_id]["strand"] = strand
                    self.new_id_map[transcript_id]["start"] = start
                    self.new_id_map[transcript_id]["end"] = end
                    self.new_id_map[transcript_id]["start_orf"] = None
                    self.new_id_map[transcript_id]["end_orf"] = None
                    self.new_id_map[transcript_id]["chromosome"] = chromosome
                    self.new_id_map[transcript_id]["synonyms"] = synonyms
                    self.new_id_map[transcript_id]["biotype"] = ""
                    self.new_id_map[transcript_id]["peptide"] = ""

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
        elif line_dict["transcript_id"].startswith("ENS"):
            return False
        else:
            return True

    @staticmethod
    def make_coord_string(transcript_dict: Dict[str, Any]):
        chromosome: str = ""
        strand: str = ""
        name: str = ""
        exon_list = list()
        for key in transcript_dict.keys():
            if key == "transcript":
                chromosome = transcript_dict[key]["seqname"]
                strand = transcript_dict[key]["strand"]
                name = transcript_dict[key]["gene_id"]
            else:
                exon = transcript_dict[key]["start"] + "-" + transcript_dict[key]["end"]
                exon_list.append(exon)
        return "_".join(sorted(exon_list)) + "_" + chromosome + "_" + strand + "_" + name


def md5_hash(data, length: int = 32):
    hasher = hashlib.md5()
    hasher.update(data.encode("utf-8"))
    hash_string = hasher.hexdigest()
    return hash_string[:length]


def main():
    print("Hello World!")


if __name__ == "__main__":
    main()
