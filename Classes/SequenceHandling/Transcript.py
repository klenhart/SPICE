#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.GTFBoy.GTFBoy import GTFBoy

from typing import Dict, Any, List


class Transcript:

    GTF_MASK: List[str] = ["seqname", "source", "feature",
                           "start", "end", "score",
                           "strand", "frame", "attribute"]

    def __init__(self) -> None:
        self.id_transcript: str = ""
        self.name_transcript: str = ""
        self.feature = "transcript"
        self.id_gene: str = ""
        self.id_taxon: int = 0
        self.biotype: str = ""
        self.transcript_support_level: int = 6
        self.tags: List[str] = list()
        self.synonyms: List[str] = list()

    def __str__(self):
        output: str = self.get_id() + " " + self.get_biotype()
        return output

    def has_sequence(self) -> bool:
        if self.feature == "transcript":
            return True
        else:
            return len(self.get_sequence()) > 0

    def set_sequence(self, seq: str) -> None:
        pass

    def set_tags(self, tag_list: List[str]):
        self.tags = tag_list

    def set_id(self, id_transcript: str) -> None:
        """

        :type id_transcript: str
        """
        self.id_transcript = id_transcript

    def set_name(self, name_transcript: str) -> None:
        self.name_transcript = name_transcript

    def set_id_taxon(self, id_taxon: int) -> None:
        """

        :type id_taxon: int
        """
        self.id_taxon = id_taxon

    def set_id_gene(self, id_protein: str) -> None:
        """

        :type id_protein: str
        """
        self.id_gene = id_protein

    def set_feature(self, feature: str) -> None:
        self.feature = feature

    def set_biotype(self, biotype: str) -> None:
        self.biotype = biotype

    def set_transcript_support_level(self, tsl: int):
        self.transcript_support_level = tsl

    def add_synonym(self, synonym: str):
        self.synonyms.append(synonym)

    def get_id(self) -> str:
        return self.id_transcript

    def get_name(self) -> str:
        return self.name_transcript

    def get_feature(self) -> str:
        return self.feature

    def get_id_taxon(self) -> int:
        return self.id_taxon

    def get_biotype(self) -> str:
        return self.biotype

    def get_id_gene(self) -> str:
        return self.id_gene

    def get_transcript_support_level(self) -> int:
        return self.transcript_support_level

    def get_tags(self) -> List[str]:
        return self.tags

    def has_tag(self, tag: str) -> bool:
        return tag in self.tags

    def get_sequence(self):
        return ""

    def get_synonyms(self):
        return self.synonyms

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_name(input_dict["transcript_name"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_feature(input_dict["feature"])
        self.set_id_gene(input_dict["gene_id"])
        self.set_biotype(input_dict["biotype"])
        self.set_tags((input_dict["tags"]))
        self.set_transcript_support_level(input_dict["tsl"])
        if "synonyms" in input_dict.keys():
            self.synonyms = input_dict["synonyms"]

    def to_dict(self) -> Dict[str, Any]:
        output: Dict[str, Any] = dict()
        output["_id"] = self.get_id()
        output["transcript_name"] = self.get_name()
        output["feature"] = self.get_feature()
        output["gene_id"] = self.get_id_gene()
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
                self.set_feature(entry)
            elif field_name == "attribute":
                attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(entry)
                self.set_id(attribute_dict["transcript_id"])
                self.set_tags(attribute_dict["tag"].split(";"))
                try:
                    self.set_name(attribute_dict["transcript_name"])
                except KeyError:
                    self.set_name(".")
                self.set_biotype(attribute_dict["transcript_biotype"])
                self.set_id_gene(attribute_dict["gene_id"])
                self.set_transcript_support_level(int(attribute_dict["transcript_support_level"]))

    def make_header(self) -> str:
        return "|".join([self.get_id_gene(), self.get_id(), str(self.get_id_taxon())])

    def __eq__(self, other):
        return any([self.get_id() == other.get_id(),
                    self.get_id() in other.get_synonyms(),
                    other.get_id() in self.get_synonyms()])
