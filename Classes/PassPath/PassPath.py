#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  PassPath is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PassPath is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
from typing import Dict


class PassPath:

    def __init__(self, path_dict: Dict[str, str]):
        self.path_dict: Dict[str, str] = path_dict

    def __getitem__(self, item: str) -> str:
        if item == "root":
            return self.path_dict["root"]
        else:
            return os.path.join(self.path_dict["root"], self.path_dict[item])
