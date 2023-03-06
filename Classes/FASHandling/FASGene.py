#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  FASGene is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASGene is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import List, Dict, Any

from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript


class FASGene(Gene):

    FAS_BIOTYPES: List[str] = ["protein_coding", "nonsense_mediated_decay"]

    def __init__(self) -> None:
        super().__init__()
        self.fas_distance_matrix_dict: Dict[str, Dict[str, float]] = dict()
        self.fas_complete_flag = True
        self.check_fas_status()

    def integrate_fas_dist_matrix(self, dist_matrix: Dict[str, Dict[str, float]]):
        for transcript_id_1 in dist_matrix.keys():
            if transcript_id_1 in self.fas_distance_matrix_dict.keys():
                matrix_row: Dict[str, float] = dist_matrix[transcript_id_1]
                for transcript_id_2 in matrix_row.keys():
                    if transcript_id_2 in self.fas_distance_matrix_dict[transcript_id_1].keys():
                        self.fas_distance_matrix_dict[transcript_id_1][transcript_id_2] = matrix_row[transcript_id_2]
        self.check_fas_status()

    def add_transcript(self, transcript: Transcript) -> None:
        """

        :type transcript: Transcript
        """
        self.transcripts[transcript.get_id()] = transcript
        if transcript.get_biotype() in FASGene.FAS_BIOTYPES:
            if transcript.get_id() not in self.fas_distance_matrix_dict.keys():
                self.fas_distance_matrix_dict[transcript.get_id()] = dict()
                for transcript_id in self.fas_distance_matrix_dict.keys():
                    if transcript_id == transcript.get_id():
                        self.fas_distance_matrix_dict[transcript.get_id()][transcript_id] = 1.0
                        self.fas_distance_matrix_dict[transcript_id][transcript.get_id()] = 1.0
                    elif transcript.get_biotype() == "nonsense_mediated_decay":
                        self.fas_distance_matrix_dict[transcript.get_id()][transcript_id] = 0.0
                        self.fas_distance_matrix_dict[transcript_id][transcript.get_id()] = 0.0
                    elif self.transcripts[transcript_id].get_biotype() == "nonsense_mediated_decay":
                        self.fas_distance_matrix_dict[transcript.get_id()][transcript_id] = 0.0
                        self.fas_distance_matrix_dict[transcript_id][transcript.get_id()] = 0.0
                    else:
                        self.fas_distance_matrix_dict[transcript.get_id()][transcript_id] = -1.0
                        self.fas_distance_matrix_dict[transcript_id][transcript.get_id()] = -1.0
        self.check_sequence_status()
        self.check_fas_status()

    def remove_transcript(self, transcript_id: str, from_dist_matrix_flag: bool = True) -> None:
        del self.transcripts[transcript_id]
        if from_dist_matrix_flag:
            del self.fas_distance_matrix_dict[transcript_id]
            for other_id in self.fas_distance_matrix_dict.keys():
                del self.fas_distance_matrix_dict[other_id][transcript_id]
        self.check_sequence_status()
        self.check_fas_status()

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
        output["fas_distance_matrix"] = self.fas_distance_matrix_dict
        return output

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.set_id(input_dict["_id"])
        self.set_name(input_dict["name"])
        self.set_feature(input_dict["feature"])
        self.set_id_taxon(input_dict["taxon_id"])
        self.set_chromosome(input_dict["chromosome"])
        self.set_species(input_dict["species"])
        self.set_biotype(input_dict["biotype"])
        if "fas_distance_matrix" in input_dict.keys():
            self.fas_distance_matrix_dict = input_dict["fas_distance_matrix"]
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

    def check_fas_status(self) -> None:
        self.fas_complete_flag = True
        for fas_dicts in self.fas_distance_matrix_dict.values():
            for fas_value in fas_dicts.values():
                if fas_value == -1.0:
                    self.fas_complete_flag = False
                    break

    def is_fas_complete(self) -> bool:
        return self.fas_complete_flag

