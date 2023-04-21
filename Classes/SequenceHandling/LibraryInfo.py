#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  LibraryInfo is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  LibraryInfo is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import yaml
from typing import Any, Dict


class LibraryInfo:

    def __init__(self, info_path: str) -> None:
        self.path: str = info_path
        with open(info_path, "r") as f:
            self.info_dict: Dict[str, Any] = yaml.safe_load(f)
        if self.info_dict is None:
            self.info_dict = dict()

    def __getitem__(self, key: str) -> Any:
        return self.info_dict[key]

    def __setitem__(self, key: str, item: Any) -> None:
        self.info_dict[key] = item

    def __str__(self) -> str:
        output: str = LibraryInfo.format_dict(self.info_dict)
        return output

    def save(self):
        with open(self.path, "w") as f:
            yaml.dump(self.info_dict, f)

    def set_self_path(self, new_path: str):
        self.path = new_path

    @staticmethod
    def format_dict(dictionary: Dict[str, Any], layer: int = 0) -> str:
        output: str = ""
        for key in dictionary.keys():
            if isinstance(dictionary[key], dict):
                output += layer * "\t" + key + ":\n" + LibraryInfo.format_dict(dictionary[key], layer + 1)
            else:
                output += layer * "\t" + key + ":\t" + str(dictionary[key]) + "\n"
        return output
