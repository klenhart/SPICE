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

from typing import Tuple

from Classes.PassPath.PassPath import PassPath


class ComparisonAssembler:

    def __init__(self,
                 condition_pair: Tuple[str, str],
                 result_pass_path: PassPath,
                 initial_flag: bool = False):
        result_pass_path[""]

    def __str__(self):
        pass
