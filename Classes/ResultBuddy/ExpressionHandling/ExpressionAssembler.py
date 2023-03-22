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

from typing import Dict

from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ExpressionAssembler:

    def __init__(self, transcript_set_path: str):
        self.transcript_set_path: str = transcript_set_path

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.transcript_set_path, "r") as f:
            gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(json.load(f))
        return gene_assembly

    def load(self) -> None:
        pass

    def save(self) -> None:
        pass

