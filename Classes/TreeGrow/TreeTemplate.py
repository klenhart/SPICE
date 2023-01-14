#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  TreeTemplate is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TreeTemplate is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict


class TreeTemplate:

    def __init__(self, paths_suffixes_dict: Dict[str, str]) -> None:
        self.paths_suffixes_dict: Dict[str, str] = paths_suffixes_dict
        self.paths_dict: Dict[str, str] = dict()
        self.prefix: str = ""

    def set_prefix(self, prefix: str) -> None:
        if prefix.endswith("/"):
            self.prefix = prefix
        else:
            self.prefix = prefix
        for key in list(self.paths_suffixes_dict.keys()):
            self.paths_dict[key] = self.prefix + self.paths_suffixes_dict[key]

    def get_paths_dict(self) -> Dict[str, str]:
        return self.paths_dict


def main() -> None:
    pass


if __name__ == "__main__":
    main()
