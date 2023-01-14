#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Protein is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Protein is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

"""
Created on Wed Dec 21 10:33:00 2022

@author: chris
"""

from Classes.SequenceHandling import Exon

from typing import Type
from typing import List


class Protein:

    def __init__(self) -> None:
        self.id_protein: str = ""
        self.id_gene: str = ""
        self.id_transcript: str = ""
        self.id_taxon: str = ""
        self.sequence: str = ""
        self.expression_value: float = 0
        self.annotation: List[str] = list()
        self.exons: List[Exon] = list()

    def set_id(self, id_protein: str) -> None:
        """

        :type id_protein: str
        """
        self.id_protein = id_protein

    def set_id_taxon(self, id_taxon):
        """

        :type id_taxon: str
        """
        self.id_taxon = id_taxon

    def set_id_gene(self, id_protein: str) -> None:
        """

        :type id_protein: str
        """
        self.id_gene = id_protein

    def set_id_transcript(self, id_transcript: str) -> None:
        """

        :type id_transcript: str
        """
        self.id_transcript = id_transcript

    def set_sequence(self, seq: str) -> None:
        """

        :type seq: str
        """
        self.sequence = seq

    def set_annotation(self, annotation: List[str]) -> None:
        """

        :type annotation: List[str]
        """
        self.annotation = annotation

    def set_expression_value(self, expression: float) -> None:
        """

        :type expression: float
        """
        self.expression_value = expression

    def add_exon(self, exon: Type[Exon]):
        self.exons.append(exon)

    def get_id_protein(self) -> str:
        return self.id_protein

    def get_id_gene(self) -> str:
        return self.id_gene

    def get_id_transcript(self) -> str:
        return self.id_transcript

    def get_sequence(self) -> str:
        return self.sequence

    def get_annotation(self) -> List[str]:
        return self.annotation

    def get_expression_value(self) -> float:
        return self.expression_value
