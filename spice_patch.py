#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  spice_path is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  spice_path is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import os
from pathlib import Path

from typing import Dict, Any

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo


def patch_fas_index(library_path: str):
    with open(os.path.join(library_path, "paths.json"), "r") as f:
        path_dict: Dict[str, Any] = json.load(f)
    path_dict["fas_scores"] = os.path.join("fas_data", "fas_scores")
    path_dict["fas_index"] = os.path.join("fas_data", "fas_index.json")
    with open(os.path.join(library_path, "paths.json"), "w") as f:
        json.dump(path_dict, f, indent=4)

    pass_path: PassPath = PassPath(path_dict)
    Path(pass_path["fas_scores"]).mkdir()

    index_dict: Dict[str, str] = dict()
    with open(os.path.join(library_path, "fas_data", "fas_scores.json"), "r") as f:
        fas_dict: Dict[str, Dict[str, Dict[str, float]]] = json.load(f)
    for key in fas_dict.keys():
        index_dict[key] = "fas_scores.json"

    with open(pass_path["fas_index"], "w") as f:
        json.dump(index_dict, f, indent=4)

    with open(os.path.join(pass_path["fas_scores"], "fas_scores.json"), "w") as f:
        json.dump(fas_dict, f, indent=4)

    lib_info: LibraryInfo = LibraryInfo(pass_path["info"])

    gene_assembler: GeneAssembler = GeneAssembler(lib_info["info"]["species"], str(lib_info["info"]["taxon_id"]))
    gene_assembler.load(pass_path)
    gene_assembler.save_fas(pass_path)

    os.remove(os.path.join(library_path, "fas_data", "fas_scores.json"))
    os.remove(os.path.join(library_path, "fas_data", "fas_scores", "fas_scores.json"))


def hotfix_fas_index(library_path: str):
    with open(os.path.join(library_path, "paths.json"), "r") as f:
        paths_dict = json.load(f)
    pass_path: PassPath = PassPath(paths_dict)
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens",
                                                  "9606")

    gene_assembler.load(pass_path)

    for gene in gene_assembler.get_genes():
        delete_list = list()
        for entry in gene.get_fas_dict().keys():
            if not entry.startswith("ENS"):
                delete_list.append(entry)
        for entry in delete_list:
            del gene.get_fas_dict()[entry]
        for entry in gene.get_fas_dict().keys():
            for del_entry in delete_list:
                del gene.get_fas_dict()[entry][del_entry]

    gene_assembler.save_fas(pass_path)


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--library", "--mode"],
                                                   [str, str],
                                                   ["store", "store"],
                                                   [None, None],
                                                   ["Directory of the library to be patched.",
                                                    """Name of the patch that shall be run:
                                                    fas_index: Splits the singles fas scores file into split up files
                                                    and an index."""])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()

    if argument_dict["mode"] == "fas_index":
        patch_fas_index(argument_dict["library"])
    elif argument_dict["mode"] == "hotfix_fas_index":
        hotfix_fas_index(argument_dict["library"])


if __name__ == "__main__":
    main()