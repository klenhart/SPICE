#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
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
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import math
import numpy as np

from typing import Dict, Any, List

from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ConditionAssembler:

    def __init__(self,
                 transcript_set_path: str,
                 condition_name: str,
                 initial_flag: bool = False):
        if initial_flag:
            self.transcript_set_path: str = transcript_set_path
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = transcript_set_path
            self.condition_assembly["normalization"]: str = ""
            self.condition_assembly["replicates"]: List[str] = list()
            self.condition_assembly["replicate_count"]: int = 0
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                self.condition_assembly["data"][gene_id]: Dict[str, List[Any]] = dict()
                self.condition_assembly["data"][gene_id]["ids"]: List[str] = list()
                self.condition_assembly["data"][gene_id]["biotypes"]: List[str] = list()
                self.condition_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = list()
                self.condition_assembly["data"][gene_id]["tags"]: List[List[str]] = list()
                self.condition_assembly["data"][gene_id]["expression"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_max"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_min"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_all"]: List[List[float]] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = list()
                self.condition_assembly["data"][gene_id]["expression_rel_std"]: List[float] = list()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.condition_assembly["data"][gene_id]["ids"].append(transcript.get_id())
                    biotype: str = transcript.get_biotype()
                    self.condition_assembly["data"][gene_id]["biotypes"].append(biotype)
                    tsl: int = transcript.get_transcript_support_level()
                    self.condition_assembly["data"][gene_id]["transcript_support_levels"].append(tsl)
                    self.condition_assembly["data"][gene_id]["tags"].append(transcript.get_tags())
                    self.condition_assembly["data"][gene_id]["expression"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_rel"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_rel_max"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_rel_min"].append(1.0)
                    self.condition_assembly["data"][gene_id]["expression_rel_avg"].append(1.0)
                    self.condition_assembly["data"][gene_id]["expression_all"].append([])
                    self.condition_assembly["data"][gene_id]["expression_rel_all"].append([])
                    self.condition_assembly["data"][gene_id]["expression_rel_std"].append(0.0)
        else:
            self.transcript_set_path: str = ""
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = transcript_set_path
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
                    expr_value = expr_assembly["data"][gene_id]["expression"][index]
                    rel_expr_value = expr_assembly["data"][gene_id]["expression_rel"][index]
                else:
                    expr_value = 0.0
                    rel_expr_value = 0.0

                # Add the expression to the total of the transcript.
                self.condition_assembly["data"][gene_id]["expression"][i] += expr_value

                # Append the expression to the list of all expressions.
                self.condition_assembly["data"][gene_id]["expression_all"][i].append(expr_value)

                # Append the relative expression to the list of all relative expressions.
                self.condition_assembly["data"][gene_id]["expression_rel_all"][i].append(rel_expr_value)

                rel_expr_all: List[float] = self.condition_assembly["data"][gene_id]["expression_rel_all"][i]
                rel_expr_sum: float = sum(rel_expr_all)
                repl_count: int = self.condition_assembly["replicate_count"]

                # Calculate the average relative expression.
                rel_expr_avg: float = rel_expr_sum / repl_count
                self.condition_assembly["data"][gene_id]["expression_rel_avg"][i] = rel_expr_avg

                current_min: float = self.condition_assembly["data"][gene_id]["expression_rel_min"][i]
                current_max: float = self.condition_assembly["data"][gene_id]["expression_rel_max"][i]

                # Replace the minimum and maximum relative expressions, if they change.
                if current_min > rel_expr_value:
                    self.condition_assembly["data"][gene_id]["expression_rel_min"][i] = rel_expr_value
                if current_max < rel_expr_value:
                    self.condition_assembly["data"][gene_id]["expression_rel_max"][i] = rel_expr_value

                expr_rel_std: float = math.sqrt(sum([(x - rel_expr_avg)**2 for x in rel_expr_all]) / repl_count)

                # Calculate the standard deviation for the transcript.
                self.condition_assembly["data"][gene_id]["expression_rel_std"][i] = expr_rel_std

            # Update all relative expressions.
            expr_sum: float = sum(self.condition_assembly["data"][gene_id]["expression"])
            for i, expr in enumerate(self.condition_assembly["data"][gene_id]["expression"]):
                if expr_sum == 0:
                    self.condition_assembly["data"][gene_id]["expression_rel"][i] = 0.0
                else:
                    self.condition_assembly["data"][gene_id]["expression_rel"][i] = expr / expr_sum

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.transcript_set_path, "r") as f:
            gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(json.load(f))
        return gene_assembly

    def cleanse_assembly(self):
        cleanse_dict: Dict[str, List[str]] = dict()
        for gene_id in self.condition_assembly["data"].keys():
            cleanse_dict[gene_id] = list()
            id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
            for index, transcript_id in enumerate(id_list):
                if self.condition_assembly["data"][gene_id]["expression"][index] == 0.0:
                    cleanse_dict[gene_id].append(transcript_id)
        for gene_id in cleanse_dict.keys():
            for transcript_id in cleanse_dict[gene_id]:
                id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
                index: int = id_list.index(transcript_id)
                self.condition_assembly["data"][gene_id]["ids"].pop(index)
                self.condition_assembly["data"][gene_id]["biotypes"].pop(index)
                self.condition_assembly["data"][gene_id]["transcript_support_levels"].pop(index)
                self.condition_assembly["data"][gene_id]["tags"].pop(index)
                self.condition_assembly["data"][gene_id]["expression"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_max"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_min"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_avg"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_all"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_all"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_std"].pop(index)

            if len(self.condition_assembly["data"][gene_id]["ids"]) <= 1:
                del self.condition_assembly["data"][gene_id]
