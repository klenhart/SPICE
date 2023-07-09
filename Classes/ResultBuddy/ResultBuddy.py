#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
import json

from typing import Dict, Any, List

from tqdm import tqdm

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ComparisonHandling.ComparisonAssembler import ComparisonAssembler
from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.ResultBuddy.EWFDHandling.EWFDAssembler import EWFDAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript
from Classes.TreeGrow.TreeGrow import TreeGrow
from Classes.WriteGuard.WriteGuard import WriteGuard


class ResultBuddy:

    def __init__(self, library_path: str, output_path: str, initial_flag: bool = False, suffix: str = ""):
        self.library_path: str = library_path
        with open(os.path.join(library_path, "paths.json"), "r") as f:
            self.library_pass_path: PassPath = PassPath(json.load(f))
        library_info: LibraryInfo = LibraryInfo(self.library_pass_path["info"])
        species: str = library_info["info"]["species"]
        release: str = library_info["info"]["release"]
        fas_hex: str = library_info["info"]["fas_mode"]
        if len(suffix) > 0:
            filename: str = "spice_result_" + species + "_" + release + "_" + fas_hex + "_" + suffix
        else:
            filename: str = "spice_result_" + species + "_" + release + "_" + fas_hex
        self.result_path: str = output_path + "/" + filename

        if initial_flag:
            self.result_info: Dict[str, Any] = dict()

            self.result_info["species"] = species
            self.result_info["taxon_id"] = library_info["info"]["taxon_id"]
            self.result_info["release"] = release
            self.result_info["library_integrity_flag"] = all(library_info["status"].values())

            self.result_info["expression_imports"]: Dict[str, Dict[str, Dict[str, str]]] = dict()
            self.result_info["expression_imports"]["conditions"]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["replicates"]: Dict[str, Dict[str, Dict[str, str]]] = dict()

            result_paths: Dict[str, Any] = dict()
            result_paths["library_path"] = self.library_path
            result_paths["root"] = self.result_path
            result_paths["result_info"] = "info.json"
            result_paths["expression"] = "expression"
            result_paths["expression_replicates"] = "expression/replicates"
            result_paths["expression_conditions"] = "expression/conditions"
            result_paths["ewfd"] = "ewfd"
            result_paths["ewfd_replicates"] = "ewfd/replicates"
            result_paths["ewfd_conditions"] = "ewfd/conditions"
            result_paths["comparison"] = "comparison"

            self.result_pass_path: PassPath = PassPath(result_paths)
            tree_grow: TreeGrow = TreeGrow(result_paths)
            tree_grow.create_folders()
            tree_grow.put_path_json()

            self.__save_info()

        else:
            self.result_info: Dict[str, Any] = self.__load_info()
            self.result_pass_path: PassPath = PassPath(self.__load_paths())

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

    def generate_ewfd_file(self, name: str, condition_flag: bool = False):
        if condition_flag:
            ewfd_type: str = "conditions"
        else:
            ewfd_type: str = "replicates"
        self.result_info = self.__load_info()
        species: str = self.result_info["species"]
        taxon_id: int = self.result_info["taxon_id"]
        expr_path: str = self.result_info["expression_imports"][ewfd_type][name]["expression_path"]
        ewfd: EWFDAssembler = EWFDAssembler(species, taxon_id, self.library_pass_path, expr_path, True, condition_flag)

        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, name):
            self.result_info = self.__load_info()
            filename: str = "/ewfd_" + name + ".json"
            ewfd_path: str = self.result_pass_path["ewfd_" + ewfd_type] + filename
            self.result_info["expression_imports"][ewfd_type][name]["ewfd_path"] = ewfd_path
            self.__save_info()
            ewfd.save(self.result_info["expression_imports"][ewfd_type][name]["ewfd_path"])

    def build_condition(self, condition_name: str, replicate_names: List[str]):
        self.result_info = self.__load_info()

        condition: ConditionAssembler = ConditionAssembler(self.library_pass_path,
                                                           condition_name,
                                                           True)
        for name in replicate_names:
            expression: ExpressionAssembler = ExpressionAssembler(self.library_pass_path,
                                                                  name)
            expression.load(self.result_info["expression_imports"]["replicates"][name]["expression_path"])
            condition.insert_expression(expression)

        # Keep all genes and transcripts, doesn't matter if they have expression or not.
        # condition.cleanse_assembly()

        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, condition_name):
            self.result_info = self.__load_info()
            filename: str = "/expression_" + condition_name + ".json"
            condition_path: str = self.result_pass_path["expression_conditions"] + filename
            new_condition_dict: Dict[str, Any] = {"replicates": replicate_names,
                                                  "expression_path": condition_path,
                                                  "ewfd_path": ""}
            self.result_info["expression_imports"]["conditions"][condition_name]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["conditions"][condition_name] = new_condition_dict
            self.__save_info()
            condition.save(self.result_info["expression_imports"]["conditions"][condition_name]["expression_path"])

    def compare(self, condition_pair: List[str]):
        condition_pair.sort()
        comparison: ComparisonAssembler = ComparisonAssembler(condition_pair[0],
                                                              condition_pair[1],
                                                              self.result_pass_path,
                                                              self.result_info)
        return comparison

    def import_expression_gtf(self,
                              expression_path: str,
                              expression_name: str,
                              normalization: str,
                              expression_threshold: float = 1.0) -> None:
        transcript_to_protein_dict: Dict[str, str] = self.transcript_to_protein_map()
        transcript_to_gene_dict: Dict[str, str] = self.transcript_to_gene_map()
        synonym_to_transcript_dict: Dict[str, str] = self.synonym_to_transcript_map()
        expression_gtf: GTFBoy = GTFBoy(expression_path)
        expression_assembler: ExpressionAssembler = ExpressionAssembler(self.library_pass_path,
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

                transcript_flag: bool = line_dict["feature"] in ["transcript", "novel_transcript"]
                if "transcript_id" in line_dict.keys() and transcript_flag:
                    line_dict["transcript_id"] = line_dict["transcript_id"].split(".")[0].split(":")[-1]
                    transcript_in_lib_flag: bool = line_dict["transcript_id"] in synonym_to_transcript_dict.keys()
                    if transcript_in_lib_flag:
                        line_dict["transcript_id"] = synonym_to_transcript_dict[line_dict["transcript_id"]]
                        gene_in_lib_flag: bool = line_dict["transcript_id"] in transcript_to_gene_dict.keys()
                    else:
                        gene_in_lib_flag: bool = False
                    # This implicitly checks if the transcript is PROTEIN CODING or NMD bio-typed.
                    if gene_in_lib_flag:
                        line_dict["transcript_id"] = synonym_to_transcript_dict[line_dict["transcript_id"]]
                        line_dict["gene_id"] = transcript_to_gene_dict[line_dict["transcript_id"]]
                        line_dict["protein_id"] = transcript_to_protein_dict[line_dict["transcript_id"]]
                        expression_assembler.insert_expression_dict(line_dict)
        # Keep all genes and transcripts, doesn't matter if they have expression or not.
        # expression_assembler.cleanse_assembly()
        expression_assembler.calc_relative_expression()
        filename: str = "/expression_" + expression_name + ".json"
        expression_json_path: str = self.result_pass_path["expression_replicates"] + filename
        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, expression_name):
            self.result_info = self.__load_info()
            new_expression_dict: Dict[str, str] = {"origin": expression_path,
                                                   "expression_path": expression_json_path,
                                                   "ewfd_path": ""}
            self.result_info["expression_imports"]["replicates"][expression_name]: Dict[str, Dict[str, str]] = dict()
            self.result_info["expression_imports"]["replicates"][expression_name] = new_expression_dict
            self.__save_info()

        expression_assembler.save(expression_json_path)

    def transcript_to_biotype_map(self) -> Dict[str, str]:
        transcript_to_biotype_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_biotype_map[transcript.id_transcript] = "protein_coding"
            elif isinstance(transcript, Transcript):
                transcript_to_biotype_map[transcript.id_transcript] = "nonsense_mediated_decay"
        return transcript_to_biotype_map

    def synonym_to_transcript_map(self) -> Dict[str, str]:
        synonym_to_transcript_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                synonym_to_transcript_map[transcript.get_id_transcript()] = transcript.get_id_transcript()
                for synonym in transcript.get_synonyms():
                    synonym_to_transcript_map[synonym] = transcript.get_id_transcript()
            elif isinstance(transcript, Transcript):
                synonym_to_transcript_map[transcript.get_id()] = transcript.get_id()
                for synonym in transcript.get_synonyms():
                    synonym_to_transcript_map[synonym] = transcript.get_id()
        return synonym_to_transcript_map

    def transcript_to_protein_map(self) -> Dict[str, str]:
        transcript_to_protein_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_protein_map[transcript.get_id_transcript()] = transcript.get_id()
            elif isinstance(transcript, Transcript):
                transcript_to_protein_map[transcript.get_id()] = ""
        return transcript_to_protein_map

    def transcript_to_gene_map(self) -> Dict[str, str]:
        transcript_to_gene_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_gene_map[transcript.get_id_transcript()] = transcript.get_id_gene()
            elif isinstance(transcript, Transcript):
                transcript_to_gene_map[transcript.get_id()] = transcript.get_id_gene()
        return transcript_to_gene_map


def main():
    pass


if __name__ == "__main__":
    main()
