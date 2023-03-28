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

from typing import Dict, Any

from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class MovementAssembler:

    def __init__(self,
                 species: str,
                 taxon_id: int,
                 transcript_set_path: str,
                 expression_path: str,
                 condition_flag: bool = False):
        expr_assembler: ExpressionAssembler = ExpressionAssembler(transcript_set_path,
                                                                        "",
                                                                        "",
                                                                        "",
                                                                        False,
                                                                        condition_flag)
        expr_assembler.load(expression_path)
        self.transcript_set_path: str = transcript_set_path
        self.movement_assembly: Dict[str, Any] = dict()
        self.movement_assembly["name"]: str = expr_assembler.expression_assembly["name"]
        self.movement_assembly["origin"]: str = expr_assembler.expression_assembly["origin"]
        self.movement_assembly["library"]: str = expr_assembler.expression_assembly["library"]
        self.movement_assembly["normalization"] = expr_assembler.expression_assembly["normalization"]
        self.movement_assembly["expression_threshold"] = expr_assembler.expression_assembly["expression_threshold"]
        self.movement_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

        gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
        gene_assembler.load(transcript_set_path)


#
        #for expr_assembler.load(expression_path)

