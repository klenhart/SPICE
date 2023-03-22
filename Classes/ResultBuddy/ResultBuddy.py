#!/bin/env python
import json
#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  ResultBuddy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ResultBuddy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os

from typing import Dict, Any, List

from Classes.SequenceHandling.LibraryInfo import LibraryInfo


class ResultBuddy:

    def __init__(self, library_path: str, result_path: str, initial_flag: bool = False):
        self.library_path: str = library_path
        self.result_path: str = result_path
        if initial_flag:
            library_info: LibraryInfo = LibraryInfo(os.path.join(self.library_path, "info.yaml"))
            self.result_info: Dict[str, Any] = dict()

            self.result_info["species"] = library_info["info"]["species"]
            self.result_info["taxon_id"] = library_info["info"]["taxon_id"]
            self.result_info["release"] = library_info["info"]["release"]
            self.result_info["library_integrity_flag"] = all(library_info["status"].values())

            self.result_info["expression_imports"]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["conditions"]: Dict[str, List[str]] = dict()
            self.result_info["expression_imports"]["replicates"]: Dict[str, Dict[str, str]] = dict()

            self.result_paths: Dict[str, Any] = dict()
            self.result_paths["library_path"] = self.library_path
            self.result_paths["result_path"] = self.result_path
            self.result_paths["result_info"] = os.path.join(self.result_path, "info.json")
            self.result_paths["result_paths"] = os.path.join(self.result_path, "paths.json")
            self.result_paths["expression"] = os.path.join(self.result_path, "expression")
            self.result_paths["movement"] = os.path.join(self.result_path, "movement")
            self.result_paths["comparison"] = os.path.join(self.result_path, "comparison")
            self.result_paths["expression_imports"]: Dict[str, Dict[str, Any]] = dict()
            self.result_paths["expression_imports"]["conditions"]: Dict[str, Dict[str, str]] = dict()
            self.result_paths["expression_imports"]["replicates"]: Dict[str, Dict[str, str]] = dict()

            with open(self.result_paths["result_info"], "w") as f:
                json.dump(self.result_info, f, indent=4)
            with open(self.result_paths["result_path"], "w") as f:
                json.dump(self.result_paths, f, indent=4)

        else:
            self.result_info: Dict[str, Any] = ResultBuddy.__load_info()
            self.result_paths: Dict[str, Any] = ResultBuddy.__load_paths()

    def __load_info(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "info.json"), "r") as f:
            result_info: Dict[str, Any] = json.load(f)
        return result_info

    def __load_paths(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "paths.json"), "r") as f:
            result_paths: Dict[str, Any] = json.load(f)
        return result_paths



