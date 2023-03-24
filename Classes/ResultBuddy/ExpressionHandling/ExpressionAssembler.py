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

from typing import Dict, Any, List

from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ExpressionAssembler:

    def __init__(self,
                 transcript_set_path: str,
                 expression_name: str,
                 origin_path: str,
                 normalization: str,
                 initial_flag: bool = False):
        if initial_flag:
            self.transcript_set_path: str = transcript_set_path
            self.expression_assembly: Dict[str, Any] = dict()
            self.expression_assembly["name"]: str = expression_name
            self.expression_assembly["conditions"]: List[str] = []
            self.expression_assembly["origin"]: str = origin_path
            self.expression_assembly["library"]: str = transcript_set_path
            self.expression_assembly["normalization"] = normalization
            self.expression_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                self.expression_assembly[gene_id]: Dict[str, Any] = dict()
                self.expression_assembly[gene_id]["biotypes"]: Dict[str, str] = dict()
                self.expression_assembly[gene_id]["transcript_support_levels"]: Dict[str, int] = dict()
                self.expression_assembly[gene_id]["tags"]: Dict[str, List[str]] = dict()
                self.expression_assembly[gene_id]["expression"]: Dict[str, float] = dict()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.expression_assembly[gene_id]["biotypes"][transcript.get_id()] = transcript.get_biotype()
                    tsl: int = transcript.get_transcript_support_level()
                    self.expression_assembly[gene_id]["transcript_support_levels"][transcript.get_id()] = tsl
                    self.expression_assembly[gene_id]["tags"][transcript.get_id()] = transcript.get_tags()
                    self.expression_assembly[gene_id]["expression"][transcript.get_id()] = 0.0
        else:
            self.transcript_set_path: str = ""
            self.expression_assembly: Dict[str, Any] = dict()

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.transcript_set_path, "r") as f:
            gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(json.load(f))
        return gene_assembly

    def insert_expression_dict(self, insert_dict: Dict[str, Any]):
        gene_id: str = insert_dict["gene_id"]
        transcript_id: str = insert_dict["transcript_id"]
        protein_id: str = insert_dict["protein_id"]
        expression: float = insert_dict[self.expression_assembly["normalization"]]
        if len(protein_id) == 0:
            self.expression_assembly[gene_id]["expression"][transcript_id] = expression
        else:
            self.expression_assembly[gene_id]["expression"][protein_id] = expression

    def load(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            self.expression_assembly = json.load(f)
            self.transcript_set_path = self.expression_assembly["library"]

    def save(self, output_path: str) -> None:
        with open(output_path, "r") as f:
            json.dump(self.expression_assembly, f, indent=4)
