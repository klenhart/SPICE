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
from typing import List, Dict, Any

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

    def cleanse_genes_without_fas(self):
        for gene_id in self.gene_assembly.keys():
            gene: FASGene = self.gene_assembly[gene_id]
            gene.clear_transcripts_without_fas()
            self.gene_assembly[gene_id] = gene

    def clear_empty_genes(self):
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
    fas_gene_assembler = FASGeneAssembler("human", "107")
    # fas_gene_assembler.load("C:/Users/chris/Desktop/git/root/extract.json")
    # fas_gene_assembler.integrate_fas_distance_matrix_dict("C:/Users/chris/Desktop/stats/old_fas.json")
    # fas_gene_assembler.cleanse_genes_without_fas()
    # fas_gene_assembler.save("C:/Users/chris/Desktop/git/root/extract_with_fas.json")

    fas_gene_assembler.load("C:/Users/chris/Desktop/git/root/extract_with_fas.json")
    fas_gene_assembler.clear_empty_genes()
    fas_gene_assembler.save("C:/Users/chris/Desktop/git/root/extract_with_fas_no_empty_genes.json")


if __name__ == "__main__":
    main()
