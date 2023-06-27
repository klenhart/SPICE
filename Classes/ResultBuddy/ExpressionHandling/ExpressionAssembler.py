#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ExpressionAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ExpressionAssembler is distributed in the hope that it will be useful,
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

from tqdm import tqdm

from Classes.PassPath.PassPath import PassPath
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.Transcript import Transcript


class ExpressionAssembler:

    def __init__(self,
                 library_pass_path: PassPath,
                 expression_name: str = "",
                 origin_path: str = "",
                 normalization: str = "",
                 initial_flag: bool = False,
                 expression_threshold: float = 1.0):
        if initial_flag:
            self.library_pass_path: PassPath = library_pass_path
            self.expression_assembly: Dict[str, Any] = dict()
            self.expression_assembly["name"]: str = expression_name
            self.expression_assembly["origin"]: str = origin_path
            self.expression_assembly["library"]: str = library_pass_path["root"]
            self.expression_assembly["normalization"] = normalization
            self.expression_assembly["expression_threshold"] = expression_threshold
            self.expression_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                self.expression_assembly["data"][gene_id]: Dict[str, Any] = dict()
                self.expression_assembly["data"][gene_id]["ids"]: List[str] = list()
                self.expression_assembly["data"][gene_id]["synonyms"]: List[List[str]] = list()
                self.expression_assembly["data"][gene_id]["biotypes"]: List[str] = list()
                self.expression_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = list()
                self.expression_assembly["data"][gene_id]["tags"]: List[List[str]] = list()
                self.expression_assembly["data"][gene_id]["expression"]: List[float] = list()
                self.expression_assembly["data"][gene_id]["expression_rel"]: List[float] = list()
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.expression_assembly["data"][gene_id]["ids"].append(transcript.get_id())
                    self.expression_assembly["data"][gene_id]["synonyms"].append(transcript.get_synonyms())
                    biotype: str = transcript.get_biotype()
                    self.expression_assembly["data"][gene_id]["biotypes"].append(biotype)
                    tsl: int = transcript.get_transcript_support_level()
                    self.expression_assembly["data"][gene_id]["transcript_support_levels"].append(tsl)
                    self.expression_assembly["data"][gene_id]["tags"].append(transcript.get_tags())
                    self.expression_assembly["data"][gene_id]["expression"].append(0.0)
                    self.expression_assembly["data"][gene_id]["expression_rel"].append(0.0)
        else:
            self.library_pass_path: PassPath = PassPath(dict())
            self.expression_assembly: Dict[str, Any] = dict()
            self.expression_assembly["name"]: str = expression_name
            self.expression_assembly["origin"]: str = origin_path
            self.expression_assembly["library"]: str = library_pass_path["root"]
            self.expression_assembly["normalization"]: str = ""
            self.expression_assembly["expression_threshold"] = expression_threshold
            self.expression_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

    def __len__(self) -> int:
        return len(self.expression_assembly["data"])

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        with open(self.library_pass_path["transcript_info"], "r") as f:
            info_dict: Dict[str, Dict[str, Any]] = json.load(f)
        with open(self.library_pass_path["transcript_seq"], "r") as f:
            seq_dict: Dict[str, Dict[str, Any]] = json.load(f)
        fas_dict: Dict[str, Dict[str, Any]] = dict()
        with open(self.library_pass_path["fas_index"], "r") as f1:
            fas_index: Dict[str, str] = json.load(f1)
            for path in set(fas_index.values()):
                with open(os.path.join(self.library_pass_path["fas_scores"], path), "r") as f2:
                    fas_sub_dict: Dict[str, Dict[str, Any]] = json.load(f2)
                    fas_dict.update(fas_sub_dict)
        gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(info_dict, seq_dict, fas_dict)
        return gene_assembly

    def insert_expression_dict(self, insert_dict: Dict[str, Any]):
        gene_id: str = insert_dict["gene_id"]
        transcript_id: str = insert_dict["transcript_id"]
        protein_id: str = insert_dict["protein_id"]
        expression: float = float(insert_dict[self.expression_assembly["normalization"]])
        expression_threshold: float = self.expression_assembly["expression_threshold"]
        id_list: List[str] = self.expression_assembly["data"][gene_id]["ids"]
        if len(protein_id) == 0:
            index: int = id_list.index(transcript_id)
            if expression >= expression_threshold:
                self.expression_assembly["data"][gene_id]["expression"][index] = expression
            else:
                self.expression_assembly["data"][gene_id]["expression"][index] = 0.0
        else:
            index: int = id_list.index(protein_id)
            if expression >= expression_threshold:
                self.expression_assembly["data"][gene_id]["expression"][index] = expression
            else:
                self.expression_assembly["data"][gene_id]["expression"][index] = 0.0

    def cleanse_assembly(self):
        cleanse_dict: Dict[str, List[str]] = dict()
        for gene_id in tqdm(self.expression_assembly["data"].keys(),
                            ncols=100,
                            total=len(self.expression_assembly["data"]),
                            desc=self.expression_assembly["name"] + ": extract cleanup progress"):
            cleanse_dict[gene_id] = list()
            id_list: List[str] = self.expression_assembly["data"][gene_id]["ids"]
            for index, transcript_id in enumerate(id_list):
                if self.expression_assembly["data"][gene_id]["expression"][index] == 0.0:
                    cleanse_dict[gene_id].append(transcript_id)
        for gene_id in cleanse_dict.keys():
            for transcript_id in cleanse_dict[gene_id]:
                id_list: List[str] = self.expression_assembly["data"][gene_id]["ids"]
                index: int = id_list.index(transcript_id)
                self.expression_assembly["data"][gene_id]["ids"].pop(index)
                self.expression_assembly["data"][gene_id]["synonyms"].pop(index)
                self.expression_assembly["data"][gene_id]["biotypes"].pop(index)
                self.expression_assembly["data"][gene_id]["transcript_support_levels"].pop(index)
                self.expression_assembly["data"][gene_id]["tags"].pop(index)
                self.expression_assembly["data"][gene_id]["expression"].pop(index)
                self.expression_assembly["data"][gene_id]["expression_rel"].pop(index)
            if len(self.expression_assembly["data"][gene_id]["ids"]) == 0:
                del self.expression_assembly["data"][gene_id]

    def calc_relative_expression(self):
        for gene_id in tqdm(self.expression_assembly["data"].keys(),
                            ncols=100,
                            total=len(self.expression_assembly["data"]),
                            desc=self.expression_assembly["name"] + ": relative expression calculation progress"):
            expressions: List[float] = self.expression_assembly["data"][gene_id]["expression"]
            total: float = sum(expressions)
            if total == 0.0:
                self.expression_assembly["data"][gene_id]["expression_rel"] = [0.0 for expr in expressions]
            else:
                self.expression_assembly["data"][gene_id]["expression_rel"] = [expr / total for expr in expressions]

    def load(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            self.expression_assembly = json.load(f)
        with open(os.path.join(self.expression_assembly["library"], "paths.json"), "r") as f:
            self.library_pass_path = PassPath(json.load(f))

    def save(self, output_path: str) -> None:
        with open(output_path, "w") as f:
            json.dump(self.expression_assembly, f, indent=4)

    def __str__(self) -> str:
        return str(self.expression_assembly)
