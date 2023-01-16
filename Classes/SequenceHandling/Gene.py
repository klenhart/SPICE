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


from typing import List, Dict, Any

from Classes.SearchTree.AbstractSearchTreeEntry import AbstractSearchTreeEntry
from Classes.SequenceHandling.Transcript import Transcript
from Classes.SequenceHandling.Protein import Protein


class Gene(AbstractSearchTreeEntry):

    fasta_template = ">{0}|{1}|{2}\n{3}"

    def __init__(self) -> None:
        self.id_gene: str = ""
        self.id_taxon: str = ""
        self.biotype: str = ""
        self.species: str = ""
        self.expression_value: float = 0
        self.transcripts: List[Transcript] = list()

    def set_id(self, id_gene: str) -> None:
        self.id_gene = id_gene

    def set_biotype(self, biotype: str) -> None:
        self.biotype = biotype

    def set_id_taxon(self, id_taxon: str) -> None:
        """

        :type id_taxon: str
        """
        self.id_taxon = id_taxon

    def set_species(self, species: str):
        self.species = species

    def set_expression_value(self, expression: float) -> None:
        """

        :type expression: float
        """
        self.expression_value = expression

    def add_transcript(self, transcript: Transcript) -> None:
        """

        :type transcript: Transcript
        """
        self.transcripts.append(transcript)

    def get_transcripts(self) -> List[Transcript]:
        return self.transcripts

    def get_expression_value(self) -> float:
        return self.expression_value

    def get_id(self) -> str:
        return self.id_gene

    def get_biotype(self) -> str:
        return self.biotype

    def get_id_taxon(self) -> str:
        return self.id_taxon

    def get_species(self) -> str:
        return self.species

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_species(input_dict["species"])
        self.set_expression_value(input_dict["expression_value"])
        self.set_biotype(input_dict["biotype"])
        self.set_expression_value(input_dict["expression_value"])
        for transcript_dict in input_dict["transcripts"]:
            if transcript_dict["type"] == "transcript":
                transcript: Transcript = Transcript()
                transcript.from_dict(transcript_dict)
                self.add_transcript(transcript)
            else:
                protein: Protein = Protein()
                protein.from_dict(transcript_dict)
                self.add_transcript(protein)

    def to_dict(self) -> Dict[str, Any]:
        output: Dict[str, Any]
        output: Dict[str, Any] = dict()
        output["_id"] = self.get_id()
        output["type"] = "gene"
        output["taxon_id"] = self.get_id_taxon()
        output["expression_value"] = self.get_expression_value()
        output["biotype"] = self.get_biotype()
        output["species"] = self.get_species()
        transcript_list: List[Dict[str, Any]] = []
        for transcript in self.get_transcripts():
            transcript_list.append(transcript.to_dict())
        output["transcripts"] = transcript_list
        return output

    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.get_id() == other.get_id()
        return False

    @property
    def fasta(self) -> str:
        proteins: List[Protein] = [transcript for transcript in self.transcripts if isinstance(transcript, Protein)]
        output: str = "\n".join([self.fasta_template.format(self.get_id(),
                                                            protein.get_id_transcript(),
                                                            protein.get_id(),
                                                            protein.get_sequence()) for protein in proteins])
        return output


def main():
    gene = Gene()
    gene.set_id("standin_gene_id")
    gene.set_species("human")
    gene.set_expression_value(12.9)
    gene.set_id_taxon("9606")

    protein1: Protein = Protein()
    protein1.set_id("p1")
    protein1.set_id_transcript("tp1")
    protein1.set_sequence("AAA")

    protein2: Protein = Protein()
    protein2.set_id("p2")
    protein2.set_id_transcript("tp2")
    protein2.set_sequence("CCC")

    protein3: Protein = Protein()
    protein3.set_id("p3")
    protein3.set_id_transcript("tp3")
    protein3.set_sequence("LLL")

    transcript1: Transcript = Transcript()
    transcript1.set_id("t1")

    gene.add_transcript(protein1)
    gene.add_transcript(protein2)
    gene.add_transcript(protein3)
    gene.add_transcript(transcript1)

    print(gene.fasta)


if __name__ == "__main__":
    main()
