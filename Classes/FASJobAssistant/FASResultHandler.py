#!/bin/env python
import json
import os
from typing import Dict, Any

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.WriteGuard.WriteGuard import WriteGuard


#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  FASResultHandler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASResultHandler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--pairings_path", "--gene_id", "--outdir", "--mode", "--anno_dir"],
                                                   [str, str, str, str, str],
                                                   ["store", "store", "store", "store", "store"],
                                                   [None, None, None, None, None],
                                                   ["Path to the pairings tsv.",
                                                    "Gene id to operate on.",
                                                    "Directory the FAS results will get stored in.",
                                                    """What operation shall be done on the results? 
                                                    'unpack', 'concat' or 'delete'""",
                                                    """Annotation directory that also contains
                                                     the concatenated FAS output."""])

    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()

    if argument_dict['mode'] == "unpack":
        with open(argument_dict["pairings_path"], "r") as f:
            gene_id_txt: str = json.load(f)[argument_dict['gene_id']]
        with open(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + ".txt"), "w") as f:
            f.write(gene_id_txt)
    elif argument_dict['mode'] == "delete":
        os.remove(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + ".txt"))
        os.remove(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + "_forward.domains"))
        os.remove(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + "_reverse.domains"))
        os.remove(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + ".phyloprofile"))
    elif argument_dict['mode'] == "concat":
        with WriteGuard(os.path.join(argument_dict["anno_dir"], "fas.phyloprofile"), argument_dict["anno_dir"]):
            with open(os.path.join(argument_dict["anno_dir"], "fas.phyloprofile", "a"):
        pass


if __name__ == "__main__":
    main()
