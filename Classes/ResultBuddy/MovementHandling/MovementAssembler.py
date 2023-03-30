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

from typing import Dict, Any, List

from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class MovementAssembler:

    def __init__(self,
                 species: str,
                 taxon_id: int,
                 transcript_set_path: str,
                 condition_path: str,
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
                self.movement_assembly["data"][gene_id]["ids"]: List[str] = list()
                self.movement_assembly["data"][gene_id]["biotypes"]: List[str] = list()
                self.movement_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = list()
                self.movement_assembly["data"][gene_id]["tags"]: List[List[str]] = list()
                self.movement_assembly["data"][gene_id]["movement"]: List[float] = list()
                self.movement_assembly["data"][gene_id]["movement_min"]: List[float] = list()
                self.movement_assembly["data"][gene_id]["movement_max"]: List[float] = list()
                self.movement_assembly["data"][gene_id]["movement_avg"]: List[float] = list()



                self.movement_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = list()
                self.movement_assembly["data"][gene_id]["expression_all"]: List[List[float]] = list()
                self.movement_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = list()
                self.movement_assembly["data"][gene_id]["expression_rel_std"]: List[float] = list()

    @staticmethod
    def calculate_movement(gene_fas_dists: Dict[str, Dict[str, float]],
                           rel_expressions: List[float],
                           transcript_ids: List[str]) -> List[float]:
        movement_list: List[float] = [0.0] * len(transcript_ids)

        for s, seed_id in enumerate(transcript_ids):
            for q, query_id in enumerate(transcript_ids):
                movement_list[s] += rel_expressions[q] * gene_fas_dists[seed_id][query_id]

        return movement_list


