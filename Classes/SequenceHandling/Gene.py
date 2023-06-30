#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
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
        self.fas_complete_flag = False
        self.fas_dict: Dict[str, Dict[str, float]] = dict()
        self.check_sequence_status()
        self.check_fas_status()

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

    def add_transcript(self, transcript: Transcript, initial_add: bool = False) -> None:
        """
        :param initial_add: True if this transcript is being newly added to the gene or if it comes from an old load.
        :type transcript: Transcript
        """
        self.transcripts[transcript.get_id()] = transcript
        if transcript.get_biotype() in ["protein_coding", "nonsense_mediated_decay", "non_coding"] and initial_add:
            if transcript.get_id() not in self.fas_dict.keys():
                self.fas_dict[transcript.get_id()] = dict()
                for transcript_id in self.fas_dict.keys():
                    if transcript.get_biotype() in ["nonsense_mediated_decay", "non_coding"]:
                        self.fas_dict[transcript.get_id()][transcript_id] = 0.0
                        self.fas_dict[transcript_id][transcript.get_id()] = 0.0
                    elif self.transcripts[transcript_id].get_biotype() in ["nonsense_mediated_decay", "non_coding"]:
                        self.fas_dict[transcript.get_id()][transcript_id] = 0.0
                        self.fas_dict[transcript_id][transcript.get_id()] = 0.0
                    else:
                        self.fas_dict[transcript.get_id()][transcript_id] = -1.0
                        self.fas_dict[transcript_id][transcript.get_id()] = -1.0

        self.check_sequence_status()
        self.check_fas_status()

    def get_transcripts(self, no_sequence_flag: bool = False) -> List[Transcript]:
        if no_sequence_flag:
            return [transcript for transcript in list(self.transcripts.values()) if transcript.has_sequence()]
        else:
            return list(self.transcripts.values())

    def get_proteins(self,
                     no_sequence_flag: bool = False,
                     no_fas_flag: bool = False) -> List[Protein]:
        if no_sequence_flag:
            return [protein for protein in self.transcripts.values() if
                    isinstance(protein, Protein) and not protein.has_sequence()]
        elif no_fas_flag:
            output_list: List[Protein] = list()
            for protein_id in self.fas_dict.keys():
                if -1 in self.fas_dict[protein_id].values():
                    output_list.append(self.transcripts[protein_id])
            return output_list
        else:
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

    def get_fas_dict(self) -> Dict[str, Dict[str, float]]:
        return self.fas_dict

    def reset_fas(self) -> None:
        for transcript_1 in self.get_transcripts():
            for transcript_2 in self.get_transcripts():
                if transcript_1.get_biotype() == "protein_coding" and transcript_1 != transcript_2:
                    self.fas_dict[transcript_1.get_id()][transcript_2.get_id()] = -1.0

    def is_sequence_complete(self) -> bool:
        return self.sequences_complete_flag

    def is_fas_complete(self) -> bool:
        return self.fas_complete_flag

    def set_fas_dict(self, fas_dict: Dict[str, Dict[str, float]]):
        self.fas_dict = fas_dict

    def set_sequence_of_transcript(self, transcript_id: str, sequence: str) -> None:
        protein: Protein = self.transcripts[transcript_id]
        protein.set_sequence(sequence)
        self.check_sequence_status()

    def check_sequence_status(self) -> None:
        self.sequences_complete_flag = all([len(protein.get_sequence()) > 0
                                            for protein in self.get_proteins() if isinstance(protein, Protein)])

    def check_fas_status(self) -> None:
        self.fas_complete_flag = all([-1 not in list(row_dict.values()) for row_dict in self.get_fas_dict().values()])

    def from_dict(self,
                  info_dict: Dict[str, Any],
                  seq_dict: Dict[str, Any],
                  fas_dict: Dict[str, Any]) -> None:
        self.set_id(info_dict["_id"])
        self.set_name(info_dict["name"])
        self.set_feature(info_dict["feature"])
        self.set_id_taxon(info_dict["taxon_id"])
        self.set_chromosome(info_dict["chromosome"])
        self.set_species(info_dict["species"])
        self.set_biotype(info_dict["biotype"])
        self.set_fas_dict(fas_dict)
        for key in info_dict["transcripts"].keys():
            transcript_dict = info_dict["transcripts"][key]
            if transcript_dict["biotype"] != "protein_coding":
                transcript: Transcript = Transcript()
                transcript.from_dict(transcript_dict)
                self.add_transcript(transcript)
            else:
                transcript_dict["sequence"] = seq_dict[key]
                protein: Protein = Protein()
                protein.from_dict(transcript_dict)
                self.add_transcript(protein)

    def to_dict(self, mode: str) -> Dict[str, Any]:
        output: Dict[str, Any]
        output: Dict[str, Any] = dict()
        if mode == "fas":
            output = self.get_fas_dict()
        elif mode == "seq":
            for protein in self.get_proteins():
                output[protein.get_id()] = protein.get_sequence()
        elif mode == "info":
            output["_id"] = self.get_id()
            output["name"] = self.get_name()
            output["feature"] = self.get_feature()
            output["taxon_id"] = self.get_id_taxon()
            output["chromosome"] = self.get_chromosome()
            output["species"] = self.get_species()
            output["biotype"] = self.get_biotype()
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

    def get_transcript_count(self, no_sequence_flag: bool = False) -> int:
        return len(self.get_transcripts(no_sequence_flag))

    def get_protein_count(self,
                          no_sequence_flag: bool = False,
                          no_fas_flag: bool = False) -> int:
        return len(self.get_proteins(no_sequence_flag, no_fas_flag))

    def delete_transcript(self, transcript_id: str):
        print("\tDeleting ", transcript_id)
        del self.transcripts[transcript_id]
        del self.fas_dict[transcript_id]
        for key in self.fas_dict.keys():
            del self.fas_dict[key][transcript_id]
        self.check_sequence_status()

    def calculate_implicit_fas_scores(self):
        for key_x in self.fas_dict.keys():
            for key_y in self.fas_dict.keys():
                if key_x == key_y:
                    self.fas_dict[key_x][key_y] = 1.0

    def make_pairings(self) -> str:
        protein_list: List[Protein] = self.get_proteins(False, True)
        row_list: List[str] = list()
        for protein_1 in protein_list:
            for protein_2 in protein_list:
                if protein_1 != protein_2 and protein_2.make_header_pair(protein_1) not in row_list:
                    row_list.append(protein_1.make_header_pair(protein_2))
        return "\n".join(row_list)

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
                                                            protein.get_id(),
                                                            protein.get_id_taxon(),
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
