#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  ExpressionAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ExpressionAssembler is distributed in the hope that it will be useful,
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
            self.condition_assembly["expression_threshold"]: float = 0.0
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
                self.condition_assembly["data"][gene_id]["expression_max"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_min"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_avg"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_all"]: List[List[float]] = list()
                self.condition_assembly["data"][gene_id]["expression_std"]: List[float] = list()
                self.condition_assembly["data"][gene_id]["expression_rel"]: List[float] = list()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.condition_assembly["data"][gene_id]["ids"].append(transcript.get_id())
                    biotype: str = transcript.get_biotype()
                    self.condition_assembly["data"][gene_id]["biotypes"].append(biotype)
                    tsl: int = transcript.get_transcript_support_level()
                    self.condition_assembly["data"][gene_id]["transcript_support_levels"].append(tsl)
                    self.condition_assembly["data"][gene_id]["tags"].append(transcript.get_tags())
                    self.condition_assembly["data"][gene_id]["expression"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_rel"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_max"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_min"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_avg"].append(0.0)
                    self.condition_assembly["data"][gene_id]["expression_all"].append([])
                    self.condition_assembly["data"][gene_id]["expression_std"].append(0.0)
        else:
            self.transcript_set_path: str = ""
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = transcript_set_path
            self.condition_assembly["normalization"]: str = ""
            self.condition_assembly["expression_threshold"]: float = 0.0
            self.condition_assembly["replicates"]: List[str] = list()
            self.condition_assembly["replicate_count"]: int = 0
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

    def load(self, input_path: str):
        with open(input_path, "r") as f:
            self.condition_assembly = json.load(f)

    def insert_expression(self, expression_assembler: ExpressionAssembler):
        expr_assembly: Dict[str, Any] = expression_assembler.expression_assembly
        self.condition_assembly["replicates"].append(expr_assembly["name"])
        self.condition_assembly["replicate_count"] += 1
        for gene_id in self.condition_assembly["data"].keys():
            if gene_id in expr_assembly["data"].keys():
                for index, transcript_id in enumerate(self.condition_assembly["data"][gene_id]["ids"]):
                    if transcript_id in expr_assembly["data"][gene_id]["ids"]:
                        expr_value = expr_assembly["data"][gene_id]["expression"][index]
                        rel_expr_value = expr_assembly["data"][gene_id]["expression_rel"][index]
                        self.condition_assembly["data"][gene_id]["expression"][index] += expr_value
            else:
                pass

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.transcript_set_path, "r") as f:
            gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(json.load(f))
        return gene_assembly

          # if condition_flag:
          #     current_max: float = self.expression_assembly["data"][gene_id]["expression_max"][transcript_id]
          #     current_min: float = self.expression_assembly["data"][gene_id]["expression_min"][transcript_id]
          #     if expr_value > current_max:
          #         self.expression_assembly["data"][gene_id]["expression_max"][transcript_id] = expr_value
          #     if expr_value < current_min:
          #         self.expression_assembly["data"][gene_id]["expression_min"][transcript_id] = expr_value
          #     self.expression_assembly["data"][gene_id]["expression_all"][transcript_id].append(expr_value)
          #     expr_all: List[float]
          #     expr_all = self.expression_assembly["data"][gene_id]["expression_all"][transcript_id]
          #     expr_avg: float = sum(expr_all) / self.expression_assembly["replicate_count"]
          #     self.expression_assembly["data"][gene_id]["expression_avg"][transcript_id] = expr_avg
          #     expr_std: float = math.sqrt(sum([(x - expr_avg) ** 2 for x in expr_all]) / replicate_count)
          #     self.expression_assembly["data"][gene_id]["expression_std"][transcript_id] = expr_std

        # def update(self, assembly_update: Dict[str, Any]) -> None:
        #     self.expression_assembly["normalization"] = assembly_update["normalization"]
        #     for gene_id in self.expression_assembly["data"].keys():
        #         if gene_id in assembly_update["data"].keys():
        #             for transcript_id in self.expression_assembly["data"][gene_id]["biotypes"].keys():
        #                 if transcript_id in assembly_update["data"][gene_id]["biotypes"].keys():
        #                     expr_value = assembly_update["data"][gene_id]["expression"][transcript_id]
        #                 else:
        #                     expr_value = 0.0
        #                 self.expression_assembly["data"][gene_id]["expression"][transcript_id] += expr_value
