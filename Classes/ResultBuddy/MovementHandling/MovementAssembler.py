#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  MovementAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MovementAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json

from typing import Dict, Any, List

from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class MovementAssembler:

    def __init__(self,
                 species: str = "",
                 taxon_id: int = "",
                 transcript_set_path: str = "",
                 condition_path: str = "",
                 initial_flag: bool = False):
        if initial_flag:
            condition_assembler: ConditionAssembler = ConditionAssembler(transcript_set_path)
            condition_assembler.load(condition_path)
            self.transcript_set_path: str = transcript_set_path
            self.movement_assembly: Dict[str, Any] = dict()
            self.movement_assembly["name"]: str = condition_assembler.condition_assembly["name"]
            self.movement_assembly["library"]: str = condition_assembler.condition_assembly["library"]
            self.movement_assembly["normalization"] = condition_assembler.condition_assembly["normalization"]
            self.movement_assembly["replicates"] = condition_assembler.condition_assembly["replicates"]
            self.movement_assembly["replicate_count"] = condition_assembler.condition_assembly["replicate_count"]
            self.movement_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

            gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
            gene_assembler.load(transcript_set_path)
            fas_dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = gene_assembler.get_fas_dist_matrix()
            condition_data: Dict[str, Any] = condition_assembler.condition_assembly["data"]

            for gene_id in condition_data.keys():
                gene_dist_matrix: Dict[str, Dict[str, float]] = fas_dist_matrix[gene_id]
                self.movement_assembly["data"][gene_id]: Dict[str, Any] = dict()
                transcript_ids: List[str] = condition_data[gene_id]["ids"]
                self.movement_assembly["data"][gene_id]["ids"]: List[str] = transcript_ids
                self.movement_assembly["data"][gene_id]["biotypes"]: List[str] = condition_data[gene_id]["biotypes"]
                tsl_list: List[int] = condition_data[gene_id]["transcript_support_levels"]
                self.movement_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = tsl_list
                self.movement_assembly["data"][gene_id]["tags"]: List[List[str]] = condition_data[gene_id]["tags"]

                expr_rel_avg_list: List[float] = condition_data[gene_id]["expression_rel_avg"]
                expr_rel_all_list: List[List[float]] = condition_data[gene_id]["expression_rel_all"]
                expr_rel_stds: List[float] = condition_data[gene_id]["expression_rel_std"]
                expr_rel_max_list: List[float] = condition_data[gene_id]["expression_rel_max"]
                expr_rel_min_list: List[float] = condition_data[gene_id]["expression_rel_min"]
                self.movement_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = expr_rel_avg_list
                self.movement_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = expr_rel_all_list
                self.movement_assembly["data"][gene_id]["expression_rel_std"]: List[float] = expr_rel_stds
                self.movement_assembly["data"][gene_id]["expression_rel_max"]: List[float] = expr_rel_max_list
                self.movement_assembly["data"][gene_id]["expression_rel_min"]: List[float] = expr_rel_min_list

                # These need to be calculated.
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr"]: List[float]
                self.movement_assembly["data"][gene_id]["movement_max_rel_expr"]: List[float]
                self.movement_assembly["data"][gene_id]["movement_min_rel_expr"]: List[float]
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr+std"]: List[float]
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr-std"]: List[float]
                self.movement_assembly["data"][gene_id]["movement_all_rel_expr"]: List[List[float]] = list()

                # And here comes the calculation:
                mov_avg_rel_expr: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                     expr_rel_avg_list,
                                                                                     transcript_ids)
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr"] = mov_avg_rel_expr

                mov_max_rel_expr: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                     expr_rel_max_list,
                                                                                     transcript_ids)
                self.movement_assembly["data"][gene_id]["movement_max_rel_expr"] = mov_max_rel_expr

                mov_min_rel_expr: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                     expr_rel_min_list,
                                                                                     transcript_ids)
                self.movement_assembly["data"][gene_id]["movement_min_rel_expr"] = mov_min_rel_expr

                expr_avg_plus_stds: List[float] = [expr + expr_rel_stds[i] for i, expr in enumerate(expr_rel_avg_list)]
                mov_rel_expr_plus_std: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                          expr_avg_plus_stds,
                                                                                          transcript_ids)
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr+std"] = mov_rel_expr_plus_std

                expr_avg_minus_stds: List[float] = [expr - expr_rel_stds[i] for i, expr in enumerate(expr_rel_avg_list)]
                mov_rel_expr_minus_std: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                           expr_avg_minus_stds,
                                                                                           transcript_ids)
                self.movement_assembly["data"][gene_id]["movement_avg_rel_expr-std"] = mov_rel_expr_minus_std

                for i in range(self.movement_assembly["replicate_count"]):
                    repl_rel_expr_list: List[float] = list()
                    for rel_expr_list in expr_rel_all_list:
                        repl_rel_expr_list.append(rel_expr_list[i])
                    mov_repl_rel_expr: List[float] = MovementAssembler.calculate_movement(gene_dist_matrix,
                                                                                          repl_rel_expr_list,
                                                                                          transcript_ids)
                    self.movement_assembly["data"][gene_id]["movement_all_rel_expr"].append(mov_repl_rel_expr)
        else:
            self.movement_assembly: Dict[str, Any] = dict()

    def save(self, output_path) -> None:
        with open(output_path, "w") as f:
            json.dump(self.movement_assembly, f, indent=4)

    def load(self, input_path) -> None:
        with open(input_path, "r") as f:
            self.movement_assembly = json.load(f)

    @staticmethod
    def calculate_movement(gene_fas_dists: Dict[str, Dict[str, float]],
                           rel_expressions: List[float],
                           transcript_ids: List[str]) -> List[float]:
        movement_list: List[float] = [0.0] * len(transcript_ids)
        for s, seed_id in enumerate(transcript_ids):
            for q, query_id in enumerate(transcript_ids):
                movement_list[s] += rel_expressions[q] * gene_fas_dists[seed_id][query_id]

        return movement_list
