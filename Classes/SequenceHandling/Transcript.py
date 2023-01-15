#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Transcript is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Transcript is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


from Classes.SequenceHandling import Exon

from typing import Type
from typing import List


class Transcript:

    def __init__(self) -> None:
        self.id_transcript: str = ""
        self.id_gene: str = ""
        self.id_taxon: str = ""
        self.biotype: str = ""
        self.transcript_support_level: int = -1

    def set_id_transcript(self, id_transcript: str) -> None:
        """

        :type id_transcript: str
        """
        self.id_transcript = id_transcript

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

    def set_biotype(self, biotype: str) -> None:
        self.biotype = biotype

    def set_transcript_support_level(self, tsl: int):
        self.transcript_support_level = tsl

    def get_id_transcript(self) -> str:
        return self.id_transcript

    def get_id_taxon(self) -> str:
        return self.id_taxon

    def get_biotype(self) -> str:
        return self.biotype

    def get_id_gene(self) -> str:
        return self.id_gene

    def get_transcript_support_level(self) -> int:
        return self.transcript_support_level
