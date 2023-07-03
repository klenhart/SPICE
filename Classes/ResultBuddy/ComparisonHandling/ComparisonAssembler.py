#!/bin/env python



#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ComparisonAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ComparisonAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict, Any

from Classes.PassPath.PassPath import PassPath


class ComparisonAssembler:

    def __init__(self,
                 condition_1: str,
                 condition_2: str,
                 result_pass_path: PassPath,
                 result_info: Dict[str, Any]):
        self.condition_1: str = condition_1
        self.condition_2: str = condition_2
        self.result_pass_path: PassPath = result_pass_path
        self.info: Dict[str, Any] = result_info

        self.condition_1_path: str = self.info["expression_import"]["conditions"][condition_1]["ewfd_path"]
        self.condition_2_path: str = self.info["expression_import"]["conditions"][condition_2]["ewfd_path"]


    def rmsd(self, gene_1, gene_2):


    def __str__(self):
        pass
