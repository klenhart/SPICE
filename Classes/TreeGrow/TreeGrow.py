#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  TreeGrow is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TreeGrow is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import os
#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  TreeGrow is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TreeGrow is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict

from pathlib import Path
import json


class TreeGrow:
    def __init__(self, paths_dict: Dict[str, str]) -> None:
        self.paths_dict: Dict[str, str] = paths_dict
        if "root" in self.paths_dict.keys():
            self.paths_dict["self"] = "paths.json"
        else:
            raise Exception("""Argument 'path_dict' passed to TreeGrow constructor must
            contain key-value-pair 'root': <path>. Exiting...""")

    def create_folders(self) -> None:
        root_path: str = self.paths_dict["root"]
        for path in self.paths_dict.values():
            if path == root_path:
                continue
            elif not TreeGrow.is_file(os.path.join(root_path, path)):
                Path(os.path.join(root_path, path)).mkdir(parents=True, exist_ok=True)
            else:
                Path(os.path.join(root_path, path)).parent.mkdir(parents=True, exist_ok=True)
                Path(os.path.join(root_path, path)).touch()

    def put_path_json(self) -> None:
        if "root" in self.paths_dict.keys():
            with open(os.path.join(self.paths_dict["root"], "paths.json"), "w") as json_file:
                json.dump(self.paths_dict, json_file, indent=4)

    @staticmethod
    def is_file(path: str) -> bool:
        for char in path[::-1]:
            if char == os.sep:
                return False
            elif char == ".":
                return True
        return False


def main() -> None:
    pass


if __name__ == "__main__":
    main()
