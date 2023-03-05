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


from Classes.SequenceHandling.Exon import Exon
from Classes.SequenceHandling.Transcript import Transcript
from Classes.GTFBoy.GTFBoy import GTFBoy

from typing import List, Dict, Any


class Protein(Transcript):

    def __init__(self) -> None:
        super().__init__()
        self.id_protein: str = ""
        self.sequence: str = ""
        self.feature: str = "protein"
        self.annotation: List[str] = list()

    def set_id(self, id_protein: str) -> None:
        """

        :type id_protein: str
        """
        self.id_protein = id_protein

    def set_id_transcript(self, id_transcript: str):
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

    # def add_exon(self, exon: Exon): # TODO Exons will be integrated in the future
    #     self.exons.insert_entry(exon)

    def get_id(self) -> str:
        return self.id_protein

    def get_id_transcript(self) -> str:
        return self.id_transcript

    def get_sequence(self) -> str:
        return self.sequence

    def get_annotation(self) -> List[str]:
        return self.annotation

    # def get_exons(self) -> List[AbstractSearchTreeEntry]: # TODO Exons will be integrated in the future
    #     return self.exons.flatten()

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_sequence(input_dict["sequence"])
        self.set_feature(input_dict["feature"])
        self.set_id_transcript(input_dict["transcript_id"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_id_gene(input_dict["gene_id"])
        self.set_annotation(input_dict["annotation"])
        self.set_biotype(input_dict["biotype"])
        self.set_transcript_support_level(input_dict["tsl"])
        # for exon_dict in input_dict["exons"]: # TODO Exons will be integrated in the future
        #     exon: Exon = Exon()
        #     exon.from_dict(exon_dict)
        #     self.add_exon(exon)

    def to_dict(self) -> Dict[str, Any]:
        output: Dict[str, Any] = dict()
        output["_id"] = self.get_id()
        output["type"] = "protein"
        output["gene_id"] = self.get_id_gene()
        output["transcript_id"] = self.get_id_transcript()
        output["taxon_id"] = self.get_id_taxon()
        output["sequence"] = self.get_sequence()
        output["annotation"] = self.get_annotation()
        output["biotype"] = self.get_biotype()
        output["tsl"] = self.get_transcript_support_level()
        # exon_list: List[Dict[str, Any]] = [] # TODO Exons will be integrated in the future.
        # for exon in self.get_exons():
        #     exon_list.append(exon.to_dict())
        # output["exons"] = exon_list
        return output

    def from_gtf_line(self, gtf_split_line: List[str]) -> None:
        for i in range(len(self.GTF_MASK)):
            field_name: str = self.GTF_MASK[i]
            entry: str = gtf_split_line[i]
            if field_name == "feature":
                self.set_feature("protein")
            elif field_name == "attribute":
                attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(entry)
                self.set_id(attribute_dict["protein_id"])
                try:
                    self.set_name(attribute_dict["transcript_name"])
                except KeyError:
                    self.set_name(".")
                self.set_biotype(attribute_dict["transcript_biotype"])
                self.set_id_gene(attribute_dict["gene_id"])
                self.set_id_transcript(attribute_dict["transcript_id"])
                self.set_transcript_support_level(int(attribute_dict["transcript_support_level"]))

    # def add_entry(self, entry_type: str, entry: Any) -> None: # TODO Exons will be integrated in the future
    #     if entry_type == "exon":
    #         self.add_exon(entry)

    def __eq__(self, other):
        if isinstance(other, Protein):
            return self.get_id() == other.get_id()
        return False
