#!/bin/env python

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

import json
import os
from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo
from Classes.WriteGuard.WriteGuard import WriteGuard


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
                                                     the concatenated FAS index JSON file."""])

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
    elif argument_dict['mode'] == "integrate":
        lib_path_dir = os.path.join("/".join(argument_dict["anno_dir"].split("/")[:-1]), "paths.json")
        with open(lib_path_dir, "r") as f:
            path_dict = json.load(f)
        pass_path: PassPath = PassPath(path_dict)
        lib_info: LibraryInfo = LibraryInfo(pass_path["info"])
        gene_assembler: GeneAssembler = GeneAssembler(lib_info["info"]["species"],
                                                      str(lib_info["info"]["taxon_id"]))
        gene_assembler.load(pass_path)
        with open(os.path.join(argument_dict["anno_dir"], "fas.phyloprofile")) as f_in:
            fas_score_list: List[str] = f_in.read().split("\n")[1:]
        fas_scores_dict: Dict[str, Dict[str, Dict[str, float]]] = dict()
        for line in fas_score_list:
            if len(line) == 0:
                continue
            split_line: List[str] = line.split("\t")
            seed: str = split_line[0]
            query: str = split_line[2]
            fas_1: float = float(split_line[3])
            fas_2: float = float(split_line[4])
            split_seed: List[str] = seed.split("|")
            split_query: List[str] = query.split("|")
            gene_id: str = split_seed[0]
            seed_prot_id: str = split_seed[1]
            query_prot_id: str = split_query[1]
            if gene_id not in fas_scores_dict.keys():
                fas_scores_dict[gene_id] = dict()
            if seed_prot_id not in fas_scores_dict[gene_id].keys():
                fas_scores_dict[gene_id][seed_prot_id] = dict()
            if query_prot_id not in fas_scores_dict[gene_id].keys():
                fas_scores_dict[gene_id][query_prot_id] = dict()
            fas_scores_dict[gene_id][seed_prot_id][query_prot_id] = fas_2
            fas_scores_dict[gene_id][query_prot_id][seed_prot_id] = fas_1

        for gene in gene_assembler.get_genes():
            if gene.get_id() in fas_scores_dict.keys():
                for prot_id_1 in fas_scores_dict[gene.get_id()].keys():
                    for prot_id_2 in fas_scores_dict[gene.get_id()][prot_id_1].keys():
                        gene.get_fas_dict()[prot_id_1][prot_id_2] = fas_scores_dict[gene.get_id()][prot_id_1][prot_id_2]

        gene_assembler.save_fas(pass_path)


if __name__ == "__main__":
    main()
