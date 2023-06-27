#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ConditionAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ConditionAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import os

from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ConditionAssembler:

    def __init__(self,
                 library_pass_path: PassPath,
                 condition_name: str = "",
                 initial_flag: bool = False):
        if initial_flag:
            self.library_pass_path: PassPath = library_pass_path
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = library_pass_path["root"]
            self.condition_assembly["normalization"]: str = ""
            self.condition_assembly["replicates"]: List[str] = list()
            self.condition_assembly["replicate_count"]: int = 0
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                self.condition_assembly["data"][gene_id]: Dict[str, List[Any]] = dict()
                self.condition_assembly["data"][gene_id]["ids"]: List[str] = list()
                self.condition_assembly["data"][gene_id]["synonyms"]: List[List[str]] = list()
                self.condition_assembly["data"][gene_id]["biotypes"]: List[str] = list()
                self.condition_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = list()
                self.condition_assembly["data"][gene_id]["tags"]: List[List[str]] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = list()
                self.condition_assembly["data"][gene_id]["expression_all"]: List[List[float]] = list()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.condition_assembly["data"][gene_id]["ids"].append(transcript.get_id())
                    self.condition_assembly["data"][gene_id]["synonyms"].append(transcript.get_synonyms())
                    biotype: str = transcript.get_biotype()
                    self.condition_assembly["data"][gene_id]["biotypes"].append(biotype)
                    tsl: int = transcript.get_transcript_support_level()
                    self.condition_assembly["data"][gene_id]["transcript_support_levels"].append(tsl)
                    self.condition_assembly["data"][gene_id]["tags"].append(transcript.get_tags())
                    self.condition_assembly["data"][gene_id]["expression_rel_avg"].append(1.0)
                    self.condition_assembly["data"][gene_id]["expression_rel_all"].append([])
                    self.condition_assembly["data"][gene_id]["expression_all"].append([])
        else:
            self.library_pass_path: PassPath = library_pass_path
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = library_pass_path["root"]
            self.condition_assembly["normalization"]: str = ""
            self.condition_assembly["replicates"]: List[str] = list()
            self.condition_assembly["replicate_count"]: int = 0
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

    def save(self, output_path: str):
        with open(output_path, "w") as f:
            json.dump(self.condition_assembly, f, indent=4)

    def load(self, input_path: str):
        with open(input_path, "r") as f:
            self.condition_assembly = json.load(f)

    def insert_expression(self, expression_assembler: ExpressionAssembler):
        expr_assembly: Dict[str, Any] = expression_assembler.expression_assembly
        self.condition_assembly["replicates"].append(expr_assembly["name"])
        self.condition_assembly["replicate_count"] += 1
        for gene_id in self.condition_assembly["data"].keys():
            if gene_id in expr_assembly["data"].keys():
                gene_found_flag: bool = True
            else:
                gene_found_flag: bool = False
            for i, transcript_id in enumerate(self.condition_assembly["data"][gene_id]["ids"]):
                if not gene_found_flag:
                    expr_value = 0.0
                    rel_expr_value = 0.0
                elif transcript_id in expr_assembly["data"][gene_id]["ids"]:
                    index: int = expr_assembly["data"][gene_id]["ids"].index(transcript_id)
                    rel_expr_value = expr_assembly["data"][gene_id]["expression_rel"][index]
                    expr_value = expr_assembly["data"][gene_id]["expression"][index]
                else:
                    expr_value = 0.0
                    rel_expr_value = 0.0

                # Append the relative and absolute expression to the list of all relative expressions.
                self.condition_assembly["data"][gene_id]["expression_rel_all"][i].append(rel_expr_value)
                self.condition_assembly["data"][gene_id]["expression_all"][i].append(expr_value)

                rel_expr_all: List[float] = self.condition_assembly["data"][gene_id]["expression_rel_all"][i]
                rel_expr_sum: float = sum(rel_expr_all)
                repl_count: int = self.condition_assembly["replicate_count"]

                # Calculate the average relative expression.
                rel_expr_avg: float = rel_expr_sum / repl_count
                self.condition_assembly["data"][gene_id]["expression_rel_avg"][i] = rel_expr_avg

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.library_pass_path["transcript_info"], "r") as f:
            info_dict: Dict[str, Dict[str, Any]] = json.load(f)
        with open(self.library_pass_path["transcript_seq"], "r") as f:
            seq_dict: Dict[str, Dict[str, Any]] = json.load(f)
        fas_dict: Dict[str, Dict[str, Any]] = dict()
        with open(self.library_pass_path["fas_index"], "r") as f1:
            fas_index: Dict[str, str] = json.load(f1)
            for path in set(fas_index.values()):
                with open(os.path.join(self.library_pass_path["fas_scores"], path), "r") as f2:
                    fas_sub_dict: Dict[str, Dict[str, Any]] = json.load(f2)
                    fas_dict.update(fas_sub_dict)
        gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(info_dict, seq_dict, fas_dict)
        return gene_assembly

    def cleanse_assembly(self):
        cleanse_dict: Dict[str, List[str]] = dict()
        for gene_id in self.condition_assembly["data"].keys():
            cleanse_dict[gene_id] = list()
            id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
            for index, transcript_id in enumerate(id_list):
                if self.condition_assembly["data"][gene_id]["expression_rel_avg"][index] == 0.0:
                    cleanse_dict[gene_id].append(transcript_id)
        for gene_id in cleanse_dict.keys():
            for transcript_id in cleanse_dict[gene_id]:
                id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
                index: int = id_list.index(transcript_id)
                self.condition_assembly["data"][gene_id]["ids"].pop(index)
                self.condition_assembly["data"][gene_id]["synonyms"].pop(index)
                self.condition_assembly["data"][gene_id]["biotypes"].pop(index)
                self.condition_assembly["data"][gene_id]["transcript_support_levels"].pop(index)
                self.condition_assembly["data"][gene_id]["tags"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_avg"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_all"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_all"].pop(index)

            if len(self.condition_assembly["data"][gene_id]["ids"]) == 0:
                del self.condition_assembly["data"][gene_id]
