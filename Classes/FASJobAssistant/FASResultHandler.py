#!/bin/env python
import json
import os
from typing import Dict, Any

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse


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
    argument_parser: ReduxArgParse = ReduxArgParse(["--pairings_path", "--gene_id", "--outdir", "--mode"],
                                                   [str, str, str, str],
                                                   ["store", "store", "store", "store"],
                                                   [1, 1, 1, 1],
                                                   ["Path to the pairings tsv.",
                                                    "Gene id to operate on.",
                                                    "Directory the FAS results will get stored in.",
                                                    """What operation shall be done on the results? 
                                                    'unpack', 'concat' or 'delete'"""])

    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict['pairings_path'] = argument_dict['pairings_path'][0]
    argument_dict['gene_id'] = argument_dict['gene_id'][0]
    argument_dict['outdir'] = argument_dict['outdir'][0]
    argument_dict['mode'] = argument_dict['mode'][0]

    if argument_dict['mode'] == "unpack":
        with open(argument_dict["pairings_path"], "r") as f:
            gene_id_tsv: str = json.load(f)[argument_dict['gene_id']]
        with open(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + ".tsv"), "w") as f:
            f.write(gene_id_tsv)
    elif argument_dict['mode'] == "delete":
        os.remove(os.path.join(argument_dict["outdir"], argument_dict['gene_id'] + ".tsv"))
    elif argument_dict['mode'] == "concat":
        


if __name__ == "__main__":
    main()
