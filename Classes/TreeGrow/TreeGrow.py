#!/bin/env python
import os
#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
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
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict

from pathlib import Path
import json


class TreeGrow:
    def __init__(self, paths_dict: Dict[str, str]) -> None:
        self.paths_dict: Dict[str, str] = paths_dict
        if "root" in self.paths_dict.keys():
            self.paths_dict["self"] = os.path.join(self.paths_dict["root"], "paths.json")

    def create_folders(self) -> None:
        for path in self.paths_dict.values():
            if not TreeGrow.is_file(path):
                Path(path).mkdir(parents=True, exist_ok=True)
            else:
                Path(path).parent.mkdir(parents=True, exist_ok=True)
                Path(path).touch()

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
    paths_dict: Dict[str, str] = dict()
    paths_dict["root"] = "C:\\Users\\chris\\Desktop\\git\\root"
    paths_dict["info"] = "C:\\Users\\chris\\Desktop\\git\\root\\info.tsv"
    paths_dict["library"] = "C:\\Users\\chris\\Desktop\\git\\root\\library"
    paths_dict["FAS"] = "C:\\Users\\chris\\Desktop\\git\\root\\library\\FAS"
    paths_dict["fasta"] = "C:\\Users\\chris\\Desktop\\git\\root\\library\\FAS\\sequences.fasta"
    paths_dict["Disturbance"] = "C:\\Users\\chris\\Desktop\\git\\root\\library\\Disturbance\\disturb.txt"

    tree_grow = TreeGrow(paths_dict)
    tree_grow.create_folders()
    tree_grow.put_path_json()


if __name__ == "__main__":
    main()
