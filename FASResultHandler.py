#!/bin/env python
import json
import os
from typing import Dict, Any, List

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
    argument_parser: ReduxArgParse = ReduxArgParse(["--pairings_path", "--gene_id", "--out_dir", "--mode", "--anno_dir"],
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
        with open(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + ".tsv"), "w") as f:
            f.write(gene_id_txt)
    elif argument_dict['mode'] == "delete":
        os.remove(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + ".tsv"))
        os.remove(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + "_forward.domains"))
        os.remove(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + "_reverse.domains"))
        os.remove(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + "_config.yml"))
        os.remove(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + ".phyloprofile"))
    elif argument_dict['mode'] == "concat":
        with WriteGuard(os.path.join(argument_dict["anno_dir"], "fas.phyloprofile"), argument_dict["anno_dir"]):
            with open(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + ".phyloprofile"), "r") as f_in:
                fas_scores: str = "\n".join(f_in.read().split("\n")[1:])
            with open(os.path.join(argument_dict["anno_dir"], "fas.phyloprofile"), "a") as f_out:
                f_out.write(fas_scores)
        # Stopped constantly adding the the json dictionary since loading it into memory took too much time.
        # with WriteGuard(os.path.join(argument_dict["anno_dir"], "fas_scores.json"), argument_dict["anno_dir"]):
        #     with open(os.path.join(argument_dict["anno_dir"], "fas_scores.json"), "r") as f_in:
        #         fas_scores_dict: Dict[str, Dict[str, Dict[str, float]]] = json.load(f_in)
        #     with open(os.path.join(argument_dict["out_dir"], argument_dict['gene_id'] + ".phyloprofile")) as f_in:
        #         fas_score_list: List[str] = f_in.read().split("\n")[1:]
        #     for line in fas_score_list:
        #         if len(line) == 0:
        #             continue
        #         split_line: List[str] = line.split("\t")
        #         seed: str = split_line[0]
        #         query: str = split_line[2]
        #         fas_1: float = float(split_line[3])
        #         fas_2: float = float(split_line[4])
        #         split_seed: List[str] = seed.split("|")
        #         split_query: List[str] = query.split("|")
        #         gene_id: str = split_seed[0]
        #         seed_prot_id: str = split_seed[1]
        #         query_prot_id: str = split_query[1]
        #         fas_scores_dict[gene_id][seed_prot_id][query_prot_id] = fas_2
        #         fas_scores_dict[gene_id][query_prot_id][seed_prot_id] = fas_1
        #     with open(os.path.join(argument_dict["anno_dir"], "fas_scores.json"), "w") as f_out:
        #         json.dump(fas_scores_dict, f_out, indent=4)
        with WriteGuard(os.path.join(argument_dict["anno_dir"], "forward.domains"), argument_dict["anno_dir"]):
            with open(os.path.join(argument_dict["out_dir"],
                                   argument_dict['gene_id'] + "_forward.domains"), "r") as f_in:
                domains_input: str = f_in.read()
                with open(os.path.join(argument_dict["anno_dir"], "forward.domains"), "a") as f_out:
                    f_out.write(domains_input)
        with WriteGuard(os.path.join(argument_dict["anno_dir"], "reverse.domains"), argument_dict["anno_dir"]):
            with open(os.path.join(argument_dict["out_dir"],
                                   argument_dict['gene_id'] + "_reverse.domains"), "r") as f_in:
                domains_input: str = f_in.read()
                with open(os.path.join(argument_dict["anno_dir"], "reverse.domains"), "a") as f_out:
                    f_out.write(domains_input)


if __name__ == "__main__":
    main()
