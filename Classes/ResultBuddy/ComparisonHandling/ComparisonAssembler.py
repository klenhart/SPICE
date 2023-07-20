#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ComparisonAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ComparisonAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict, Any, List, Set
import hashlib
import json
import math
import os
from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.EWFDHandling.EWFDAssembler import EWFDAssembler


class ComparisonGene:

    def __init__(self, gene_id: str,
                 data_dict_1: Dict[str, List[Any]],
                 data_dict_2: Dict[str, List[Any]],
                 biotype_filter: List[str], tag_filter: List[str],
                 fas_adjacency_matrix: Dict[str, Dict[str, float]]):
        self.gene_id = gene_id

        self.rmsd: float = 0.0

        data_dict_1 = ComparisonGene.apply_biotype_filter(data_dict_1, biotype_filter)
        data_dict_1 = ComparisonGene.apply_tag_filter(data_dict_1, tag_filter)

        data_dict_2 = ComparisonGene.apply_biotype_filter(data_dict_2, biotype_filter)
        data_dict_2 = ComparisonGene.apply_tag_filter(data_dict_2, tag_filter)

        data_dict_1, data_dict_2 = ComparisonGene.indices_expressed_in_one(data_dict_1,
                                                                           data_dict_2)

        ewfd_1: List[float] = ComparisonGene.recalculate_ewfd(data_dict_1, fas_adjacency_matrix)
        ewfd_2: List[float] = ComparisonGene.recalculate_ewfd(data_dict_2, fas_adjacency_matrix)

        if ComparisonGene.one_is_non_expressed(data_dict_1, data_dict_2):
            self.rmsd = 0.0
        else:
            self.calc_rmsd(ewfd_1, ewfd_2)

    @staticmethod
    def indices_expressed_in_one(data_dict_1, data_dict_2):
        delete_set: Set[int] = set()
        for i, entry in enumerate(data_dict_1["expression_rel_avg"]):
            if data_dict_1["expression_rel_avg"][i] == 0.0 and data_dict_2["expression_rel_avg"][i] == 0.0:
                delete_set.add(i)
        delete_list = list(delete_set)
        delete_list.sort(reverse=True)
        for index in delete_list:
            for key in data_dict_1.keys():
                data_dict_1[key].pop(index)
                data_dict_2[key].pop(index)
        return data_dict_1, data_dict_2

    @staticmethod
    def one_is_non_expressed(data_dict_1, data_dict_2):
        if all([x == 0.0 for x in data_dict_1["expression_rel_avg"]]):
            return True
        elif all([x == 0.0 for x in data_dict_2["expression_rel_avg"]]):
            return True
        else:
            return False

    def calc_rmsd(self, ewfd_1: List[float], ewfd_2: List[float]):
        squared_delta_list: List[float] = list()
        for i, _ in enumerate(ewfd_1):
            squared_delta_list.append((ewfd_1[i] - ewfd_2[i])**2)
        if len(ewfd_1) == 0:
            self.rmsd = 0.0
        else:
            self.rmsd = math.sqrt(sum(squared_delta_list)) / len(ewfd_1)

    def __str__(self):
        output: str = ",".join([self.gene_id, str(self.rmsd)])
        return output

    @staticmethod
    def recalculate_ewfd(data_dict: Dict[str, List[Any]], fas_adjacency_matrix: Dict[str, Dict[str, float]]):
        return EWFDAssembler.calculate_ewfd(fas_adjacency_matrix,
                                            data_dict["expression_rel_avg"],
                                            data_dict["ids"])

    @staticmethod
    def apply_biotype_filter(data_dict: Dict[str, List[Any]], biotype_filter: List[str]):
        delete_set = set()
        for entry in biotype_filter:
            for i, biotype in enumerate(data_dict["biotypes"]):
                if biotype == entry:
                    delete_set.add(i)
        delete_list = list(delete_set)
        delete_list.sort(reverse=True)
        for index in delete_list:
            for key in data_dict.keys():
                data_dict[key].pop(index)
        return data_dict

    @staticmethod
    def apply_tag_filter(data_dict: Dict[str, List[Any]], tag_filter: List[str]):
        delete_set = set()
        for entry in tag_filter:
            for i, tags in enumerate(data_dict["tags"]):
                if entry in tags:
                    delete_set.add(i)
        delete_list = list(delete_set)
        delete_list.sort(reverse=True)
        for index in delete_list:
            for key in data_dict.keys():
                data_dict[key].pop(index)
        return data_dict

    def __lt__(self, other):
        return self.rmsd < other.rmsd


class ComparisonAssembler:

    def __init__(self,
                 condition_1: str,
                 condition_2: str,
                 result_pass_path: PassPath,
                 result_info: Dict[str, Any]):
        self.condition_1: str = condition_1
        self.condition_2: str = condition_2
        self.result_pass_path: PassPath = result_pass_path
        self.info: Dict[str, Any] = result_info

        self.comparison_gene_list: List[ComparisonGene] = list()

        with open(os.path.join(result_pass_path["library_path"], "fas_data", "fas_index.json"), "r") as f:
            self.fas_index: Dict[str, str] = json.load(f)

        self.fas_scores_directory: str = os.path.join(result_pass_path["library_path"], "fas_data", "fas_scores")

        self.condition_1_path: str = self.info["expression_imports"]["conditions"][condition_1]["ewfd_path"]
        self.condition_2_path: str = self.info["expression_imports"]["conditions"][condition_2]["ewfd_path"]

        self.biotype_filter: List[str] = list()
        self.tag_filter: List[str] = list()
        self.filter_hash = md5_hash("".join(self.biotype_filter+self.tag_filter), 6)

    def add_biotype_filter(self, filter_out: str):
        self.biotype_filter.append(filter_out)
        self.filter_hash = md5_hash("".join(self.biotype_filter + self.tag_filter), 6)

    def add_tag_filter(self, filter_out: str):
        self.tag_filter.append(filter_out)
        self.filter_hash = md5_hash("".join(self.biotype_filter + self.tag_filter), 6)

    def compare_genes(self):
        with open(self.condition_1_path, "r") as f_1:
            data_cond_1: Dict[str, Dict[str, List[Any]]] = json.load(f_1)["data"]

        with open(self.condition_2_path, "r") as f_2:
            data_cond_2: Dict[str, Dict[str, List[Any]]] = json.load(f_2)["data"]

        for gene_id in data_cond_1.keys():
            with open(os.path.join(self.fas_scores_directory, self.fas_index[gene_id]), "r") as f:
                fas_adjacency_matrix: Dict[str, Dict[str, float]] = json.load(f)[gene_id]
            self.comparison_gene_list.append(ComparisonGene(gene_id,
                                                            data_cond_1[gene_id],
                                                            data_cond_2[gene_id],
                                                            self.biotype_filter,
                                                            self.tag_filter,
                                                            fas_adjacency_matrix))

    def sort_genes_by_rmsd(self):
        self.comparison_gene_list.sort(reverse=True)

    def delete_gene_below_rmsd(self, value: float):
        delete_list: List[int] = list()
        for i, comparison_gene in enumerate(self.comparison_gene_list):
            if comparison_gene.rmsd < value:
                delete_list.append(i)
        for i in sorted(delete_list, reverse=True):
            self.comparison_gene_list.pop(i)

    def __str__(self):
        return "\n".join([str(gene) for gene in self.comparison_gene_list])

    def save(self, out_dir):
        with open(os.path.join(out_dir,
                               self.condition_1 + "@" + self.condition_2 + "@" + self.filter_hash + ".tsv"), "w") as f:
            f.write(str(self))


def md5_hash(data, length: int = 32):
    hasher = hashlib.md5()
    hasher.update(data.encode("utf-8"))
    hash_string = hasher.hexdigest()
    return hash_string[:length]

