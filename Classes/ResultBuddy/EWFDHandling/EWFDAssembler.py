#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  EWFDAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EWFDAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import math
import numpy as np

from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class EWFDAssembler:

    def __init__(self,
                 species: str = "",
                 taxon_id: int = "",
                 library_pass_path: PassPath = PassPath(dict()),
                 expression_path: str = "",
                 initial_flag: bool = False,
                 condition_flag: bool = False):
        if initial_flag:
            if condition_flag:
                condition_assembler: ConditionAssembler = ConditionAssembler(library_pass_path)
                condition_assembler.load(expression_path)
                self.library_pass_path: PassPath = library_pass_path
                self.ewfd_assembly: Dict[str, Any] = dict()
                self.ewfd_assembly["name"]: str = condition_assembler.condition_assembly["name"]
                self.ewfd_assembly["library"]: str = condition_assembler.condition_assembly["library"]
                self.ewfd_assembly["normalization"] = condition_assembler.condition_assembly["normalization"]
                self.ewfd_assembly["replicates"] = condition_assembler.condition_assembly["replicates"]
                self.ewfd_assembly["replicate_count"] = condition_assembler.condition_assembly["replicate_count"]
                self.ewfd_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

                gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
                gene_assembler.load(self.library_pass_path)
                fas_dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = gene_assembler.get_fas_dist_matrix()
                condition_data: Dict[str, Any] = condition_assembler.condition_assembly["data"]

                for gene_id in condition_data.keys():
                    gene_dist_matrix: Dict[str, Dict[str, float]] = fas_dist_matrix[gene_id]
                    self.ewfd_assembly["data"][gene_id]: Dict[str, Any] = dict()
                    transcript_ids: List[str] = condition_data[gene_id]["ids"]
                    self.ewfd_assembly["data"][gene_id]["ids"]: List[str] = transcript_ids
                    transcript_count: int = len(self.ewfd_assembly["data"][gene_id]["ids"])
                    self.ewfd_assembly["data"][gene_id]["biotypes"]: List[str] = condition_data[gene_id]["biotypes"]
                    tsl_list: List[int] = condition_data[gene_id]["transcript_support_levels"]
                    self.ewfd_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = tsl_list
                    self.ewfd_assembly["data"][gene_id]["tags"]: List[List[str]] = condition_data[gene_id]["tags"]

                    expr_rel_avg_list: List[float] = condition_data[gene_id]["expression_rel_avg"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = expr_rel_avg_list

                    expr_rel_all_list: List[float] = condition_data[gene_id]["expression_rel_all"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel_all"]: List[float] = expr_rel_all_list

                    # These need to be calculated.
                    self.ewfd_assembly["data"][gene_id]["ewfd_avg_rel_expr"]: List[float]
                    self.ewfd_assembly["data"][gene_id]["ewfd_all"]: List[List[float]] = list()
                    self.ewfd_assembly["data"][gene_id]["ewfd_max"]: List[float] = [0.0] * transcript_count
                    self.ewfd_assembly["data"][gene_id]["ewfd_min"]: List[float] = [1.0] * transcript_count
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd"]: List[float] = list()
                    self.ewfd_assembly["data"][gene_id]["ewfd_std"]: List[float] = list()
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd+std"]: List[float] = list()
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd-std"]: List[float] = list()

                    # And here comes the calculation:
                    ewfd_avg_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                                  expr_rel_avg_list,
                                                                                  transcript_ids)
                    self.ewfd_assembly["data"][gene_id]["ewfd_avg_rel_expr"] = ewfd_avg_rel_expr

                    ewfd_all: List[List[float]] = list()
                    for i, repl_rel_expr_list in enumerate(np.array(expr_rel_all_list).transpose()):
                        ewfd_repl_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                                       list(repl_rel_expr_list),
                                                                                       transcript_ids)
                        ewfd_all.append(ewfd_repl_rel_expr)

                        # Check all EWFD values for this replicate if they exceed the current min max bounds.
                        for j, ewfd in enumerate(ewfd_repl_rel_expr):
                            if ewfd > self.ewfd_assembly["data"][gene_id]["ewfd_max"][j]:
                                self.ewfd_assembly["data"][gene_id]["ewfd_max"][j] = ewfd
                            if ewfd < self.ewfd_assembly["data"][gene_id]["ewfd_min"][j]:
                                self.ewfd_assembly["data"][gene_id]["ewfd_min"][j] = ewfd

                    ewfd_all = list(np.array(ewfd_all).transpose())
                    self.ewfd_assembly["data"][gene_id]["ewfd_all"] = [list(entry) for entry in ewfd_all]

                    # Calculate the average EWFD for each transcript.
                    for ewfd_list in self.ewfd_assembly["data"][gene_id]["ewfd_all"]:
                        count: int = len(ewfd_list)
                        ewfd_average: float = sum(ewfd_list) / count
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd"].append(ewfd_average)
                        std: float
                        std = math.sqrt(sum([(ewfd_value - ewfd_average)**2 for ewfd_value in ewfd_list]) / count)
                        self.ewfd_assembly["data"][gene_id]["ewfd_std"].append(std)
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd+std"].append(ewfd_average+std)
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd-std"].append(ewfd_average-std)
            else:
                expression_assembler: ExpressionAssembler = ExpressionAssembler(library_pass_path)
                expression_assembler.load(expression_path)
                self.library_pass_path: PassPath = library_pass_path
                self.ewfd_assembly: Dict[str, Any] = dict()
                self.ewfd_assembly["name"]: str = expression_assembler.expression_assembly["name"]
                self.ewfd_assembly["library"]: str = expression_assembler.expression_assembly["library"]
                self.ewfd_assembly["normalization"] = expression_assembler.expression_assembly["normalization"]
                self.ewfd_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

                gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
                gene_assembler.load(self.library_pass_path)
                fas_dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = gene_assembler.get_fas_dist_matrix()
                expression_data: Dict[str, Any] = expression_assembler.expression_assembly["data"]

                for gene_id in expression_data.keys():
                    gene_dist_matrix: Dict[str, Dict[str, float]] = fas_dist_matrix[gene_id]
                    self.ewfd_assembly["data"][gene_id]: Dict[str, Any] = dict()
                    transcript_ids: List[str] = expression_data[gene_id]["ids"]
                    self.ewfd_assembly["data"][gene_id]["ids"]: List[str] = transcript_ids
                    self.ewfd_assembly["data"][gene_id]["biotypes"]: List[str] = expression_data[gene_id]["biotypes"]
                    tsl_list: List[int] = expression_data[gene_id]["transcript_support_levels"]
                    self.ewfd_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = tsl_list
                    self.ewfd_assembly["data"][gene_id]["tags"]: List[List[str]] = expression_data[gene_id]["tags"]
                    expression_rel: List[float] = expression_data[gene_id]["expression_rel"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel"]: List[float] = expression_rel

                    # These need to be calculated.
                    self.ewfd_assembly["data"][gene_id]["ewfd_rel_expr"]: List[float]

                    # And here comes the calculation:
                    ewfd_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                              expression_rel,
                                                                              transcript_ids)
                    self.ewfd_assembly["data"][gene_id]["ewfd_rel_expr"] = ewfd_rel_expr
        else:
            self.ewfd_assembly: Dict[str, Any] = dict()

    def save(self, output_path) -> None:
        with open(output_path, "w") as f:
            json.dump(self.ewfd_assembly, f, indent=4)

    def load(self, input_path) -> None:
        with open(input_path, "r") as f:
            self.ewfd_assembly = json.load(f)

    @staticmethod
    def calculate_ewfd(gene_fas_dists: Dict[str, Dict[str, float]],
                       rel_expressions: List[float],
                       transcript_ids: List[str]) -> List[float]:
        ewfd_list: List[float] = [0.0] * len(transcript_ids)
        # This calculates the movement.
        for s, seed_id in enumerate(transcript_ids):
            for q, query_id in enumerate(transcript_ids):
                ewfd_list[s] += round(rel_expressions[q] * (1 - gene_fas_dists[seed_id][query_id]), 4)

        # ewfd_list = [round(1 - movement, 4) for movement in ewfd_list]

        return ewfd_list
