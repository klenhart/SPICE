#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  FASGeneAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASGeneAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import numpy as np
from typing import List, Dict, Any, Tuple

from tqdm import tqdm

from Classes.FASHandling.FASGene import FASGene
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class FASGeneAssembler(GeneAssembler):

    def __init__(self, species: str, taxon_id: str):
        super().__init__(species, taxon_id)
        self.gene_assembly: Dict[str, FASGene] = dict()

    def save(self, output_path: str) -> None:
        json_dict: Dict[str, Dict[str, Any]] = FASGeneAssembler.to_dict(self.gene_assembly)
        with open(output_path, "w") as f:
            json.dump(json_dict, f, indent=4)

    def load(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            self.gene_assembly = FASGeneAssembler.from_dict(json.load(f))

    def integrate_fas_distance_matrix_dict(self, input_path: str):
        with open(input_path, "r") as f:
            fas_distance_matrix_dict: Dict[str, Dict[str, Dict[str, float]]] = json.load(f)

        for gene_id in self.gene_assembly.keys():
            if gene_id in fas_distance_matrix_dict.keys():
                gene: FASGene = self.gene_assembly[gene_id]
                distance_matrix: Dict[str, Dict[str, float]] = fas_distance_matrix_dict[gene_id]
                gene.integrate_fas_dist_matrix(distance_matrix)
                self.gene_assembly[gene_id] = gene

    def get_fas_grouped_by_tsl(self):
        fas_values_grouped_by_tsl: List[List[float]] = [[], [], [], [], [], []]
        for gene_id in tqdm(self.gene_assembly.keys(),
                            ncols=100,
                            total=len(self.gene_assembly.keys()),
                            desc="Collecting FAS data grouped by TSL"):
            gene: FASGene = self.gene_assembly[gene_id]
            tsl_id_grouping_dict: Dict[str, List[str]] = {"0": [], "1": [], "2": [],
                                                          "3": [], "4": [], "5": [],
                                                          "6": [], "no_prot": []}
            for transcript in gene.get_transcripts():
                if transcript.get_biotype() != "protein_coding":
                    tsl_id_grouping_dict["no_prot"].append(transcript.get_id())
                else:
                    tsl_id_grouping_dict[str(transcript.get_transcript_support_level())].append(transcript.get_id())

            for i, tsl_list in enumerate([["0", "1"], ["2"], ["3"], ["4"], ["5"], ["6"]]):
                for tsl in tsl_list:
                    for transcript_id in tsl_id_grouping_dict[tsl]:
                        dist_row = gene.fas_distance_matrix_dict[transcript_id]
                        for key_id in dist_row.keys():
                            if key_id not in tsl_id_grouping_dict["no_prot"]:
                                fas_values_grouped_by_tsl[i].append(dist_row[key_id])
        return fas_values_grouped_by_tsl

    def get_protein_coding_transcript_counts(self, threshold: int = 7,
                                             include_nmd_flag: bool = False) -> List[List[int]]:
        if include_nmd_flag:
            possible_biotypes: List[str] = ["nonsense_mediated_decay", "protein_coding"]
        else:
            possible_biotypes: List[str] = ["protein_coding"]
        count_dict: Dict[int, int] = dict()
        for gene_id in tqdm(self.gene_assembly.keys(),
                            ncols=100,
                            total=len(self.gene_assembly.keys()),
                            desc="Counting splicing variants for each gene"):
            gene: FASGene = self.gene_assembly[gene_id]
            count: int = len([1 for transcript in gene.get_transcripts()
                              if transcript.get_transcript_support_level() < threshold and
                              transcript.get_biotype() in possible_biotypes])
            if count in count_dict.keys():
                count_dict[count] += 1
            else:
                count_dict[count] = 1

        count_list: List[int] = list(count_dict.items())
        count_list.sort()

        output_list: List[List[int]] = [[], []]
        for key, value in count_list:
            output_list[0].append(key)
            output_list[1].append(value)
        return output_list

    def get_all_fas_values(self, threshold: int = 7, include_nmd_flag: bool = False) -> List[float]:
        if include_nmd_flag:
            possible_biotypes: List[str] = ["nonsense_mediated_decay", "protein_coding"]
        else:
            possible_biotypes: List[str] = ["protein_coding"]
        fas_values: List[float] = []
        for gene_id in tqdm(self.gene_assembly.keys(),
                            ncols=100,
                            total=len(self.gene_assembly.keys()),
                            desc="Collecting FAS scores"):
            gene: FASGene = self.gene_assembly[gene_id]

            for protein_id in gene.fas_distance_matrix_dict.keys():
                protein_biotype: str = gene.transcripts[protein_id].get_biotype()
                protein_tsl: int = gene.transcripts[protein_id].get_transcript_support_level()
                if protein_tsl < threshold and protein_biotype in possible_biotypes:
                    dist_row = gene.fas_distance_matrix_dict[protein_id]
                    for key_id in dist_row.keys():
                        biotype: str = gene.transcripts[key_id].get_biotype()
                        tsl: int = gene.transcripts[key_id].get_transcript_support_level()
                        is_self: bool = key_id == protein_id
                        has_fas: bool = dist_row[key_id] != -1.0
                        if tsl < threshold and biotype in possible_biotypes and has_fas and not is_self:
                            fas_values.append(dist_row[key_id])
        return fas_values

    def get_nonsense_mediated_decay_gene_ratio(self) -> List[int]:
        gene_count_list: List[int] = [0, 0]
        for gene_id in tqdm(self.gene_assembly.keys(),
                            ncols=100,
                            total=len(self.gene_assembly.keys()),
                            desc="Collecting NMD ratio"):
            gene: FASGene = self.gene_assembly[gene_id]
            is_nmd_flag: bool = False
            for transcript in gene.get_transcripts():
                if transcript.get_biotype() == "nonsense_mediated_decay":
                    is_nmd_flag = True
                    break
            if is_nmd_flag:
                gene_count_list[0] += 1
            else:
                gene_count_list[1] += 1

        return gene_count_list

    def clear_empty_genes(self) -> None:
        gene_list: List[Gene] = self.get_genes()
        for gene in gene_list:
            if len(gene.get_transcripts()) == 0:
                del self.gene_assembly[gene.get_id()]

    @staticmethod
    def from_dict(json_dict: Dict[str, Dict[str, Any]]) -> Dict[str, FASGene]:
        output_dict: Dict[str, FASGene] = dict()
        for key in json_dict.keys():
            new_gene: FASGene = FASGene()
            new_gene.from_dict(json_dict[key])
            output_dict[new_gene.get_id()] = new_gene
        return output_dict

    @staticmethod
    def to_dict(gene_assembly: Dict[str, FASGene]) -> Dict[str, Dict[str, Any]]:
        json_dict: Dict[str, Dict[str, Any]] = dict()
        for key in gene_assembly.keys():
            json_dict[key] = gene_assembly[key].to_dict()
        return json_dict


def main():
    fas_gene_assembler = FASGeneAssembler("human", "NCBI9606")
    fas_gene_assembler.load("C:/Users/chris/Desktop/git/root/extract.json")
    fas_gene_assembler.integrate_fas_distance_matrix_dict("C:/Users/chris/Desktop/stats/old_fas.json")
    fas_gene_assembler.clear_empty_genes()
    fas_gene_assembler.save("C:/Users/chris/Desktop/git/root/extract_without_fas.json")


if __name__ == "__main__":
    main()
