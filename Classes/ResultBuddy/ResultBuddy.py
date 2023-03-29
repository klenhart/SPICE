#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  ResultBuddy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ResultBuddy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
import json

from typing import Dict, Any, List

from tqdm import tqdm

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript
from Classes.TreeGrow.TreeGrow import TreeGrow
from Classes.WriteGuard.WriteGuard import WriteGuard


class ResultBuddy:

    def __init__(self, library_path: str, result_path: str, initial_flag: bool = False):
        self.library_path: str = library_path
        self.result_path: str = result_path
        if initial_flag:
            library_info: LibraryInfo = LibraryInfo(os.path.join(self.library_path, "info.yaml"))
            self.result_info: Dict[str, Any] = dict()

            self.result_info["species"] = library_info["info"]["species"]
            self.result_info["taxon_id"] = library_info["info"]["taxon_id"]
            self.result_info["release"] = library_info["info"]["release"]
            self.result_info["library_integrity_flag"] = all(library_info["status"].values())

            self.result_info["expression_imports"]: Dict[str, Dict[str, Dict[str, str]]] = dict()
            self.result_info["expression_imports"]["conditions"]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["replicates"]: Dict[str, Dict[str, Dict[str, str]]] = dict()

            self.result_paths: Dict[str, Any] = dict()
            self.result_paths["library_path"] = self.library_path
            self.result_paths["root"] = self.result_path
            self.result_paths["result_info"] = os.path.join(self.result_path, "info.json")
            self.result_paths["expression"] = os.path.join(self.result_path, "expression")
            self.result_paths["replicates"] = os.path.join(self.result_path, "expression", "replicates")
            self.result_paths["conditions"] = os.path.join(self.result_path, "expression", "conditions")
            self.result_paths["movement"] = os.path.join(self.result_path, "movement")
            self.result_paths["comparison"] = os.path.join(self.result_path, "comparison")

            tree_grow: TreeGrow = TreeGrow(self.result_paths)
            tree_grow.create_folders()
            tree_grow.put_path_json()

            self.__save_info()

        else:
            self.result_info: Dict[str, Any] = self.__load_info()
            self.result_paths: Dict[str, Any] = self.__load_paths()

    def __load_info(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "info.json"), "r") as f:
            result_info: Dict[str, Any] = json.load(f)
        return result_info

    def __save_info(self) -> None:
        with open(os.path.join(self.result_path, "info.json"), "w") as f:
            json.dump(self.result_info, f, indent=4)

    def __load_paths(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "paths.json"), "r") as f:
            result_paths: Dict[str, Any] = json.load(f)
        return result_paths

    def build_condition(self, condition_name: str, replicate_names: List[str]):
        self.result_info = self.__load_info()
        condition: ConditionAssembler = ConditionAssembler(os.path.join(self.library_path,
                                                                        "transcript_data",
                                                                        "transcript_set.json"),
                                                           condition_name,
                                                           True)
        for name in replicate_names:
            expression: ExpressionAssembler = ExpressionAssembler(os.path.join(self.library_path,
                                                                               "transcript_data",
                                                                               "transcript_set.json"),
                                                                  name)
            expression.load(self.result_info["expression_imports"]["replicates"][name]["path"])
            condition.insert_expression(expression)

        condition.cleanse_assembly()

        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path):
            self.result_info = self.__load_info()
            new_condition_dict: Dict[str, Any] = {"replicates": replicate_names,
                                                  "path": os.path.join(self.result_paths["conditions"],
                                                                       condition_name + ".json")}
            self.result_info["expression_imports"]["conditions"][condition_name]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["conditions"][condition_name] = new_condition_dict
            self.__save_info()
            condition.save(self.result_info["expression_imports"]["conditions"][condition_name]["path"])

    def import_expression_gtf(self, expression_path: str,
                              expression_name: str,
                              normalization: str,
                              expression_threshold: float = 1.0) -> None:
        transcript_to_protein_dict: Dict[str, str] = self.transcript_to_protein_map()
        expression_gtf: GTFBoy = GTFBoy(expression_path)
        expression_assembler: ExpressionAssembler = ExpressionAssembler(os.path.join(self.library_path,
                                                                                     "transcript_data",
                                                                                     "transcript_set.json"),
                                                                        expression_name,
                                                                        expression_path,
                                                                        normalization,
                                                                        True,
                                                                        expression_threshold)
        for line in tqdm(expression_gtf,
                         ncols=100,
                         total=expression_gtf.total_lines,
                         desc=expression_name + " GTF extraction progress"):
            if line.startswith("#"):
                continue
            else:
                split_line: List[str] = line.split("\t")
                line_dict: Dict[str, str] = GTFBoy.build_dict(split_line)
                transcript_flag: bool = line_dict["feature"] == "transcript"
                if all([key in line_dict.keys() for key in ["reference_id", "ref_gene_id"]]) and transcript_flag:
                    line_dict["gene_id"] = line_dict["ref_gene_id"].split(".")[0]
                    line_dict["transcript_id"] = line_dict["reference_id"].split(".")[0]
                    # This implicitly checks if the transcript is PROTEIN CODING or NMD bio-typed.
                    if line_dict["transcript_id"] in transcript_to_protein_dict.keys():
                        line_dict["protein_id"] = transcript_to_protein_dict[line_dict["transcript_id"]]
                        expression_assembler.insert_expression_dict(line_dict)
        expression_assembler.cleanse_assembly()
        expression_assembler.calc_relative_expression()
        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path):
            self.result_info = self.__load_info()
            new_expression_dict: Dict[str, str] = {"origin": expression_path,
                                                   "path": os.path.join(self.result_paths["replicates"],
                                                                        expression_name + ".json")}
            self.result_info["expression_imports"]["replicates"][expression_name]: Dict[str, Dict[str, str]] = dict()
            self.result_info["expression_imports"]["replicates"][expression_name] = new_expression_dict
            self.__save_info()
        expression_assembler.save(self.result_info["expression_imports"]["replicates"][expression_name]["path"])

    def transcript_to_protein_map(self) -> Dict[str, str]:
        transcript_to_protein_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(os.path.join(self.library_path, "transcript_data", "transcript_set.json"))
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_protein_map[transcript.get_id_transcript()] = transcript.get_id()
            elif isinstance(transcript, Transcript):
                transcript_to_protein_map[transcript.get_id()] = ""
        return transcript_to_protein_map


def main():

    library_path: str = "C:/Users/chris/Desktop/git/fade_lib_homo_sapiens_107"
    result_buddy: ResultBuddy = ResultBuddy(library_path, "C:/Users/chris/Desktop/git/result", True)
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF023EXJ.gtf", "EXJ", "FPKM")
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF082OHO.gtf", "OHO", "FPKM")
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF180OMN.gtf", "OMN", "FPKM")
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF263YFG.gtf", "YFG", "FPKM")
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF277PKW.gtf", "PKW", "FPKM")
    result_buddy.import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF304JRO.gtf", "JRO", "FPKM")

    # result_buddy: ResultBuddy = ResultBuddy(library_path, "C:/Users/chris/Desktop/git/result")
    result_buddy.build_condition("COND_0", ["EXJ", "OHO"])
    result_buddy.build_condition("COND_13", ["OMN", "JRO"])
    result_buddy.build_condition("COND_2", ["YFG", "PKW"])
    result_buddy.build_condition("COND_ALL", ["YFG", "PKW", "EXJ", "OHO", "OMN", "JRO"])


if __name__ == "__main__":
    main()
