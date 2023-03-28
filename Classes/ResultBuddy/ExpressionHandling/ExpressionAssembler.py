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

from typing import Dict, Any, List

from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ExpressionAssembler:

    def __init__(self,
                 transcript_set_path: str,
                 expression_name: str,
                 origin_path: str = "",
                 normalization: str = "",
                 initial_flag: bool = False,
                 condition_flag: bool = False,
                 expression_threshold: float = 1.0):
        if initial_flag:
            self.transcript_set_path: str = transcript_set_path
            self.expression_assembly: Dict[str, Any] = dict()
            self.expression_assembly["name"]: str = expression_name
            self.expression_assembly["origin"]: str = origin_path
            self.expression_assembly["library"]: str = transcript_set_path
            self.expression_assembly["normalization"] = normalization
            self.expression_assembly["expression_threshold"] = expression_threshold
            self.expression_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            if condition_flag:
                self.expression_assembly["replicate_count"]: int = 0
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                self.expression_assembly["data"][gene_id]: Dict[str, Any] = dict()
                self.expression_assembly["data"][gene_id]["biotypes"]: Dict[str, str] = dict()
                self.expression_assembly["data"][gene_id]["transcript_support_levels"]: Dict[str, int] = dict()
                self.expression_assembly["data"][gene_id]["tags"]: Dict[str, List[str]] = dict()
                self.expression_assembly["data"][gene_id]["expression"]: Dict[str, float] = dict()
                if condition_flag:
                    self.expression_assembly["data"][gene_id]["expression_max"]: Dict[str, float] = dict()
                    self.expression_assembly["data"][gene_id]["expression_min"]: Dict[str, float] = dict()
                    self.expression_assembly["data"][gene_id]["expression_avg"]: Dict[str, float] = dict()
                    self.expression_assembly["data"][gene_id]["expression_all"]: Dict[str, List[float]] = dict()
                    self.expression_assembly["data"][gene_id]["expression_std"]: Dict[str, float] = dict()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    biotype: str = transcript.get_biotype()
                    self.expression_assembly["data"][gene_id]["biotypes"][transcript.get_id()] = biotype
                    tsl: int = transcript.get_transcript_support_level()
                    self.expression_assembly["data"][gene_id]["transcript_support_levels"][transcript.get_id()] = tsl
                    self.expression_assembly["data"][gene_id]["tags"][transcript.get_id()] = transcript.get_tags()
                    self.expression_assembly["data"][gene_id]["expression"][transcript.get_id()] = 0.0
                    if condition_flag:
                        self.expression_assembly["data"][gene_id]["expression_max"][transcript.get_id()] = 0.0
                        self.expression_assembly["data"][gene_id]["expression_min"][transcript.get_id()] = 0.0
                        self.expression_assembly["data"][gene_id]["expression_avg"][transcript.get_id()] = 0.0
                        self.expression_assembly["data"][gene_id]["expression_all"][transcript.get_id()] = []
                        self.expression_assembly["data"][gene_id]["expression_std"][transcript.get_id()] = 0.0
        else:
            self.transcript_set_path: str = ""
            self.expression_assembly: Dict[str, Any] = dict()
            self.expression_assembly["name"]: str = expression_name
            self.expression_assembly["origin"]: str = origin_path
            self.expression_assembly["library"]: str = transcript_set_path
            self.expression_assembly["normalization"]: str = ""
            self.expression_assembly["expression_threshold"] = expression_threshold
            self.expression_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

    def __len__(self) -> int:
        return len(self.expression_assembly["data"])

    def update(self, assembly_update: Dict[str, Any], condition_flag: bool = False) -> None:
        replicate_count: int
        if condition_flag:
            self.expression_assembly["replicate_count"] += 1
            replicate_count: int = self.expression_assembly["replicate_count"]
        else:
            replicate_count = 1
        self.expression_assembly["normalization"] = assembly_update["normalization"]
        for gene_id in self.expression_assembly["data"].keys():
            if gene_id in assembly_update["data"].keys():
                for transcript_id in self.expression_assembly["data"][gene_id]["biotypes"].keys():
                    if transcript_id in assembly_update["data"][gene_id]["biotypes"].keys():
                        expr_value = assembly_update["data"][gene_id]["expression"][transcript_id]
                    else:
                        expr_value = 0.0
                    self.expression_assembly["data"][gene_id]["expression"][transcript_id] += expr_value
                    if condition_flag:
                        current_max: float = self.expression_assembly["data"][gene_id]["expression_max"][transcript_id]
                        current_min: float = self.expression_assembly["data"][gene_id]["expression_min"][transcript_id]
                        if expr_value > current_max:
                            self.expression_assembly["data"][gene_id]["expression_max"][transcript_id] = expr_value
                        if expr_value < current_min:
                            self.expression_assembly["data"][gene_id]["expression_min"][transcript_id] = expr_value
                        self.expression_assembly["data"][gene_id]["expression_all"][transcript_id].append(expr_value)
                        expr_all: List[float]
                        expr_all = self.expression_assembly["data"][gene_id]["expression_all"][transcript_id]
                        expr_avg: float = sum(expr_all) / self.expression_assembly["replicate_count"]
                        self.expression_assembly["data"][gene_id]["expression_avg"][transcript_id] = expr_avg
                        expr_std: float = math.sqrt(sum([(x - expr_avg)**2 for x in expr_all]) / replicate_count)
                        self.expression_assembly["data"][gene_id]["expression_std"][transcript_id] = expr_std
            else:
                for transcript_id in self.expression_assembly["data"][gene_id]["biotypes"].keys():
                    expr_value = 0.0
                    self.expression_assembly["data"][gene_id]["expression"][transcript_id] += expr_value
                    if condition_flag:
                        current_max: float = self.expression_assembly["data"][gene_id]["expression_max"][transcript_id]
                        current_min: float = self.expression_assembly["data"][gene_id]["expression_min"][transcript_id]
                        if expr_value > current_max:
                            self.expression_assembly["data"][gene_id]["expression_max"][transcript_id] = expr_value
                        if expr_value < current_min:
                            self.expression_assembly["data"][gene_id]["expression_min"][transcript_id] = expr_value
                        self.expression_assembly["data"][gene_id]["expression_all"][transcript_id].append(expr_value)
                        expr_all: List[float]
                        expr_all = self.expression_assembly["data"][gene_id]["expression_all"][transcript_id]
                        expr_avg: float = sum(expr_all) / self.expression_assembly["replicate_count"]
                        self.expression_assembly["data"][gene_id]["expression_avg"][transcript_id] = expr_avg
                        expr_std: float = math.sqrt(sum([(x - expr_avg) ** 2 for x in expr_all]) / replicate_count)
                        self.expression_assembly["data"][gene_id]["expression_std"][transcript_id] = expr_std

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.transcript_set_path, "r") as f:
            gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(json.load(f))
        return gene_assembly

    def insert_expression_dict(self, insert_dict: Dict[str, Any]):
        gene_id: str = insert_dict["gene_id"]
        transcript_id: str = insert_dict["transcript_id"]
        protein_id: str = insert_dict["protein_id"]
        expression: float = float(insert_dict[self.expression_assembly["normalization"]])
        expression_threshold: float = self.expression_assembly["expression_threshold"]
        if len(protein_id) == 0:
            if expression >= expression_threshold:
                self.expression_assembly["data"][gene_id]["expression"][transcript_id] = expression
            else:
                self.expression_assembly["data"][gene_id]["expression"][transcript_id] = 0.0
        else:
            if expression >= expression_threshold:
                self.expression_assembly["data"][gene_id]["expression"][protein_id] = expression
            else:
                self.expression_assembly["data"][gene_id]["expression"][protein_id] = 0.0

    def cleanse_assembly(self, condition_flag: bool = False):
        cleanse_dict: Dict[str, List[str]] = dict()
        for gene_id in self.expression_assembly["data"].keys():
            cleanse_dict[gene_id] = list()
            for transcript_id in self.expression_assembly["data"][gene_id]["expression"].keys():
                if self.expression_assembly["data"][gene_id]["expression"][transcript_id] == 0.0:
                    cleanse_dict[gene_id].append(transcript_id)
        for gene_id in cleanse_dict.keys():
            for transcript_id in cleanse_dict[gene_id]:
                del self.expression_assembly["data"][gene_id]["biotypes"][transcript_id]
                del self.expression_assembly["data"][gene_id]["transcript_support_levels"][transcript_id]
                del self.expression_assembly["data"][gene_id]["tags"][transcript_id]
                del self.expression_assembly["data"][gene_id]["expression"][transcript_id]
                if condition_flag:
                    del self.expression_assembly["data"][gene_id]["expression_max"][transcript_id]
                    del self.expression_assembly["data"][gene_id]["expression_min"][transcript_id]
                    del self.expression_assembly["data"][gene_id]["expression_all"][transcript_id]
                    del self.expression_assembly["data"][gene_id]["expression_avg"][transcript_id]
                    del self.expression_assembly["data"][gene_id]["expression_std"][transcript_id]
            if len(self.expression_assembly["data"][gene_id]["biotypes"]) == 0:
                del self.expression_assembly["data"][gene_id]

    def load(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            self.expression_assembly = json.load(f)
            self.transcript_set_path = self.expression_assembly["library"]

    def save(self, output_path: str) -> None:
        with open(output_path, "w") as f:
            json.dump(self.expression_assembly, f, indent=4)

    def __str__(self) -> str:
        return str(self.expression_assembly)
