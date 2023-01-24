#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Exon is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Exon is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.SearchTree.AbstractSearchTreeEntry import AbstractSearchTreeEntry
from Classes.GTFBoy.GTFBoy import GTFBoy

from typing import Dict, Any, List


class Exon(AbstractSearchTreeEntry):

    GTF_MASK: List[str] = ["seqname", "source", "feature",
                           "start", "end", "score",
                           "strand", "frame", "attribute"]

    def __init__(self) -> None:
        self.id_exon: str = ""
        self.begin: int = 0
        self.end: int = 0
        self.id_protein: str = ""
        self.id_gene: str = ""

    def set_id(self, id_exon: str) -> None:
        """

        :type id_exon: str
        """
        self.id_exon = id_exon

    def set_id_gene(self, id_gene: str) -> None:
        self.id_gene = id_gene

    def set_id_protein(self, id_protein: str) -> None:
        self.id_protein = id_protein

    def set_begin(self, begin: int) -> None:
        self.begin = begin

    def set_end(self, end: int) -> None:
        self.end = end

    def get_id(self) -> str:
        return self.id_exon

    def get_begin(self) -> int:
        return self.begin

    def get_end(self) -> int:
        return self.end

    def get_id_gene(self) -> str:
        return self.id_gene

    def get_id_protein(self) -> str:
        return self.id_protein

    def from_dict(self, input_dict: Dict[str, Any]):
        self.set_id(input_dict["_id"])
        self.set_begin(input_dict["begin"])
        self.set_end(input_dict["end"])
        self.set_id_gene("gene_id")
        self.set_id_protein("protein_id")

    def to_dict(self) -> Dict[str, Any]:
        output: Dict[str, Any] = dict()
        output["_id"] = self.get_id()
        output["type"] = "exon"
        output["begin"] = self.get_begin()
        output["end"] = self.get_end()
        output["protein_id"] = self.get_id_protein()
        output["gene_id"] = self.get_id_gene()
        return output

    def from_gtf_line(self, gtf_split_line: List[str]):
        for i in range(len(self.GTF_MASK)):
            field_name: str = self.GTF_MASK[i]
            entry: str = gtf_split_line[i]

            if field_name == "start":
                self.set_begin(int(entry))
            elif field_name == "end":
                self.set_end(int(entry))
            elif field_name == "attribute":
                attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(entry)
                self.set_id(attribute_dict["exon_id"])

    def add_entry(self, entry_type: str, entry: Any) -> None:
        pass  # Exons have no entries yet.

    def __eq__(self, other):
        if isinstance(other, Exon):
            return self.get_id() == other.get_id()
        return False
