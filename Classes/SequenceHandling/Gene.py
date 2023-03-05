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

from Classes.SequenceHandling.Transcript import Transcript
from Classes.SequenceHandling.Protein import Protein
from Classes.GTFBoy.GTFBoy import GTFBoy


class Gene:

    GTF_MASK: List[str] = ["seqname", "source", "feature",
                           "start", "end", "score",
                           "strand", "frame", "attribute"]

    fasta_template = ">{0}|{1}|{2}\n{3}"

    def __init__(self) -> None:
        self.id_gene: str = ""
        self.name_gene: str = ""
        self.id_taxon: str = ""
        self.feature: str = "gene"
        self.chromosome: str = ""
        self.biotype: str = ""
        self.species: str = ""
        self.transcripts: Dict[str, Transcript] = dict()
        self.sequences_complete_flag = False
        self.check_sequence_status()

    def set_id(self, id_gene: str) -> None:
        self.id_gene = id_gene

    def set_name(self, name_gene: str) -> None:
        self.name_gene = name_gene

    def set_feature(self, feature: str) -> None:
        self.feature = feature

    def set_biotype(self, biotype: str) -> None:
        self.biotype = biotype

    def set_id_taxon(self, id_taxon: str) -> None:
        """

        :type id_taxon: str
        """
        self.id_taxon = id_taxon

    def set_species(self, species: str):
        self.species = species

    def set_chromosome(self, chromosome: str) -> None:
        self.chromosome = chromosome

    def add_transcript(self, transcript: Transcript) -> None:
        """

        :type transcript: Transcript
        """
        self.transcripts[transcript.get_id()] = transcript
        self.check_sequence_status()

    def get_transcripts(self) -> List[Transcript]:
        return list(self.transcripts.values())

    def get_proteins(self) -> List[Protein]:
        return [protein for protein in self.transcripts.values() if isinstance(protein, Protein)]

    def get_id(self) -> str:
        return self.id_gene

    def get_name(self) -> str:
        return self.name_gene

    def get_feature(self) -> str:
        return self.feature

    def get_biotype(self) -> str:
        return self.biotype

    def get_id_taxon(self) -> str:
        return self.id_taxon

    def get_species(self) -> str:
        return self.species

    def get_chromosome(self) -> str:
        return self.chromosome

    def is_sequence_complete(self) -> bool:
        return self.sequences_complete_flag

    def set_sequence_of_transcript(self, transcript_id: str, sequence: str) -> None:
        protein: Protein = self.transcripts[transcript_id]
        protein.set_sequence(sequence)
        self.check_sequence_status()

    def check_sequence_status(self) -> None:
        self.sequences_complete_flag = all([len(protein.get_sequence()) > 0
                                            for protein in self.get_transcripts() if isinstance(protein, Protein)])

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_name(input_dict["name"])
        self.set_feature(input_dict["feature"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_chromosome(input_dict["chromosome"])
        self.set_species(input_dict["species"])
        self.set_biotype(input_dict["biotype"])
        for key in input_dict["transcripts"].keys():
            transcript_dict = input_dict["transcripts"][key]
            if transcript_dict["feature"] == "transcript":
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
        output["name"] = self.get_name()
        output["feature"] = self.get_feature()
        output["taxon_id"] = self.get_id_taxon()
        output["chromosome"] = self.get_chromosome()
        output["biotype"] = self.get_biotype()
        output["species"] = self.get_species()
        transcript_dict: Dict[str, Dict[str, Any]] = dict()
        for transcript in self.get_transcripts():
            transcript_dict[transcript.get_id()] = transcript.to_dict()
        output["transcripts"] = transcript_dict
        return output

    def from_gtf_line(self, gtf_split_line: List[str]):
        for i in range(len(self.GTF_MASK)):
            field_name: str = self.GTF_MASK[i]
            entry: str = gtf_split_line[i]

            if field_name == "seqname":
                self.set_chromosome(entry)
            elif field_name == "feature":
                self.set_feature(entry)
            elif field_name == "attribute":
                attribute_dict: Dict[str, str] = GTFBoy.build_attribute_dict(entry)
                self.set_id(attribute_dict["gene_id"])
                try:
                    self.set_name(attribute_dict["gene_name"])
                except KeyError:
                    self.set_name(".")
                self.set_biotype(attribute_dict["gene_biotype"])

    def add_entry(self, entry_type: str, entry: Any):
        if entry_type == "transcript":
            self.add_transcript(entry)
        # elif entry_type == "exon": # TODO Exons will be integrated in the future
        #     self.transcripts.find(entry.get_id_protein()).add_entry("exon", entry)

    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.get_id() == other.get_id()
        return False

    def __getitem__(self, item: str) -> Any:
        key_list: List[str] = ["_id", "name", "feature",
                               "taxon_id", "chromosome", "species", "biotype", "transcripts"]
        if item in key_list:
            return self.attributes[item]

    @property
    def attributes(self) -> Dict[str, Any]:
        output_dict: Dict[str, Any] = {"_id": self.get_id(),
                                       "name": self.get_name(),
                                       "feature": self.get_feature(),
                                       "taxon_id": self.get_id_taxon(),
                                       "chromosome": self.get_chromosome,
                                       "species": self.get_species(),
                                       "biotype": self.get_biotype(),
                                       "transcripts": self.get_transcripts(),
                                       "proteins": self.get_proteins()}
        return output_dict

    @property
    def fasta(self) -> str:
        proteins: List[Protein] = [transcript for transcript in self.get_transcripts() if isinstance(transcript,
                                                                                                     Protein)]
        output: str = "\n".join([self.fasta_template.format(self.get_id(),
                                                            protein.get_id_transcript(),
                                                            protein.get_id(),
                                                            protein.get_sequence()) for protein in proteins])
        return output


def main():
    gene = Gene()
    gene.set_id("standin_gene_id")
    gene.set_species("human")
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
