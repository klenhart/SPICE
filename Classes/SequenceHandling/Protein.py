#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


from Classes.SequenceHandling.Transcript import Transcript
from Classes.GTFBoy.GTFBoy import GTFBoy

from typing import List, Dict, Any


class Protein(Transcript):

    def __init__(self) -> None:
        super().__init__()
        self.id_protein: str = ""
        self.sequence: str = ""
        self.feature: str = "protein"

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

    def get_id(self) -> str:
        return self.id_protein

    def get_id_transcript(self) -> str:
        return self.id_transcript

    def get_sequence(self) -> str:
        return self.sequence

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_sequence(input_dict["sequence"])
        self.set_feature(input_dict["feature"])
        self.set_id_transcript(input_dict["transcript_id"])
        self.set_name(input_dict["transcript_name"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_id_gene(input_dict["gene_id"])
        self.set_biotype(input_dict["biotype"])
        self.set_tags(input_dict["tags"])
        self.set_transcript_support_level(input_dict["tsl"])
        if "synonyms" in input_dict.keys():
            self.synonyms = input_dict["synonyms"]

    def to_dict(self) -> Dict[str, Any]:
        output: Dict[str, Any] = dict()
        output["_id"] = self.get_id()
        output["feature"] = self.get_feature()
        output["gene_id"] = self.get_id_gene()
        output["transcript_name"] = self.get_name()
        output["transcript_id"] = self.get_id_transcript()
        output["taxon_id"] = self.get_id_taxon()
        output["biotype"] = self.get_biotype()
        output["tags"] = self.get_tags()
        output["tsl"] = self.get_transcript_support_level()
        output["synonyms"] = self.get_synonyms()
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
                self.set_tags(attribute_dict["tag"].split(";"))
                try:
                    self.set_name(attribute_dict["transcript_name"])
                except KeyError:
                    self.set_name(".")
                self.set_biotype(attribute_dict["transcript_biotype"])
                self.set_id_gene(attribute_dict["gene_id"])
                self.set_id_transcript(attribute_dict["transcript_id"])
                self.set_transcript_support_level(int(attribute_dict["transcript_support_level"]))

    def make_header_pair(self, other: Transcript) -> str:
        return self.make_header() + "\t" + other.make_header()

    def __eq__(self, other) -> bool:
        return any([self.get_id() == other.get_id(),
                    self.get_id() in other.get_synonyms(),
                    other.get_id() in self.get_synonyms()])

    def __len__(self) -> int:
        return len(self.sequence)
