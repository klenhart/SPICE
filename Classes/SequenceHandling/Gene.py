#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Gene is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Gene is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

"""
Created on Wed Dec 21 10:03:16 2022

@author: Christian Bluemel
"""

from Classes.SequenceHandling import Protein


from typing import Type
from typing import List


class Gene:

    fasta_template = ">{0}|{1}|{2}\n{3}\n"

    def __init__(self) -> None:
        self.id_gene: str = ""
        self.id_taxon: str = ""
        self.species: str = ""
        self.expression_value: float = 0
        self.sequences: Type[Protein] = list()

    def set_id_gene(self, id_gene: str) -> None:
        self.id_gene = id_gene

    def set_id_taxon(self, id_taxon: str) -> None:
        """

        :type id_taxon: str
        """
        self.id_taxon = id_taxon

    def set_species(self, species: str):
        self.species = species

    def set_expression(self, expression: float) -> None:
        """

        :type expression: float
        """
        self.expression_value = expression

    def set_expression_value(self, expression: float) -> None:
        """

        :type expression: float
        """
        self.expression_value = expression

    def add_sequence(self, seq: Type[Protein]) -> None:
        """

        :type seq: Protein
        """
        self.sequences.append(seq)

    def get_sequences(self) -> List[Type[Protein]]:
        return self.sequences

    def get_expression_value(self) -> float:
        return self.expression_value

    def get_id_gene(self) -> str:
        return self.id_gene

    def get_id_taxon(self) -> str:
        return self.id_taxon

    def get_species(self) -> str:
        return self.species

    @property
    def to_fasta(self) -> str:
        output: str = ""
        for protein in self.sequences:
            output += self.fasta_template.format(self.get_id_gene(),
                                                 protein.get_id_transcript(),
                                                 protein.get_id_protein(),
                                                 protein.get_sequence())
        return output
