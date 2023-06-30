#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  FastaBoy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FastaBoy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import sys
from typing import Dict, List, Tuple, Iterator, Any
import argparse
import json

from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript


class SpiceFastaBoy:

    def __init__(self, fasta_path: str, taxon_id: int):
        self.fasta_path: str = fasta_path
        self.fasta_dict: Dict[str, List[Transcript]] = dict()
        self.taxon_id: int = taxon_id

    def __iter__(self) -> Iterator[str]:
        with open(self.fasta_path, "r") as f:
            for line in f:
                yield line

    def get_fasta_dict(self):
        return self.fasta_dict

    def parse_fasta(self):
        gene_id: str = ""
        sequence: str = ""
        for line in self:
            if line.startswith(">"):
                if gene_id != "" and isinstance(self.fasta_dict[gene_id][-1], Protein):
                    self.fasta_dict[gene_id][-1].set_sequence(sequence)
                sequence = ""
                fasta_header_dict: Dict[str, str] = self.make_fasta_header_dict(line)
                gene_id: str = fasta_header_dict["gene_id"]
                transcript_id: str = fasta_header_dict["transcript_id"]
                if gene_id not in self.fasta_dict.keys():
                    self.fasta_dict[gene_id] = list()
                if fasta_header_dict["biotype"] in ["non_coding", "nonsense_mediated_decay"]:
                    transcript = Transcript()
                elif fasta_header_dict["biotype"] in ["protein_coding"]:
                    fasta_header_dict["sequence"] = ""
                    transcript = Protein()
                else:
                    print("Biotype " + fasta_header_dict["biotype"] + " of transcript " + transcript_id + "unknown.")
                    sys.exit(1)
                transcript.from_dict(fasta_header_dict)
                self.fasta_dict[gene_id].append(transcript)
                continue
            sequence += line.strip()

    def make_fasta_header_dict(self, header: str) -> Dict[str, Any]:
        header_dict: Dict[str, Any] = dict()
        gene_id, transcript_id, tags, synonyms = header.split("|")

        header_dict["gene_id"] = gene_id[1:]
        header_dict["_id"] = transcript_id
        header_dict["transcript_id"] = transcript_id
        header_dict["transcript_name"] = transcript_id

        tag_list: List[str] = tags.split(" ")
        delete_index = -1
        for i, tag in enumerate(tag_list):
            if tag.startswith("biotype:"):
                header_dict["biotype"] = tag[8:]
                delete_index: int = i
                break
        if delete_index == -1:
            print("Biotype not defined in tags for transcript", transcript_id,
                  "; must include biotype:<BIOTYPE> in header tags.")
            sys.exit(1)
        else:
            tag_list.pop(delete_index)
        header_dict["tags"] = tag_list

        header_dict["synonyms"] = [entry.strip() for entry in synonyms.split(" ")]
        if len(header_dict["synonyms"]) == 1:
            if header_dict["synonyms"][0] == "":
                header_dict["synonyms"] = list()

        header_dict["feature"] = "novel_transcript"
        header_dict["taxon_id"] = self.taxon_id
        header_dict["tsl"] = 6
        return header_dict


class TransDecoderFastaBoy:

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.fasta_dict: Dict[str, List[Dict[str, str]]] = dict()

    def __iter__(self) -> Iterator[str]:
        with open(self.fasta_path, "r") as f:
            for line in f:
                yield line

    def get_fasta_dict(self) -> Dict[str, List[Dict[str, str]]]:
        return self.fasta_dict

    def parse_fasta(self):
        transcript_id: str = ""
        for line in self:
            if line.startswith(">"):
                fasta_header_dict: Dict[str, str] = TransDecoderFastaBoy.make_fasta_header_dict(line)
                transcript_id: str = fasta_header_dict["transcript_id"]
                del fasta_header_dict["transcript_id"]
                fasta_header_dict["sequence"] = ""
                if transcript_id not in self.fasta_dict.keys():
                    self.fasta_dict[transcript_id] = list()
                self.fasta_dict[transcript_id].append(fasta_header_dict)
                continue
            self.fasta_dict[transcript_id][-1]["sequence"] += line.strip().replace("*", "")

    @staticmethod
    def make_fasta_header_dict(fasta_header: str) -> Dict[str, str]:
        entry_list: List[str] = [entry for entry in fasta_header.split(" ")]
        fasta_header_dict: Dict[str, str] = dict()

        for entry in entry_list:
            if entry.startswith(">"):
                fasta_header_dict["transcript_id"] = entry[1:].split(".")[0]
                fasta_header_dict["transcript_tag"] = entry[1:].split(".")[1]
            elif entry.startswith("type:"):
                fasta_header_dict["orf_type"] = entry.split(":")[1]
            elif entry.startswith("gc:"):
                fasta_header_dict["orf_type"] = entry.split(":")[1]
            elif entry.startswith(fasta_header_dict["transcript_id"]):
                indices: List[str] = entry.split(":")[1].split("-")
                fasta_header_dict["start_orf"] = indices[0]
                fasta_header_dict["end_orf"] = indices[1].split("(")[0]
                fasta_header_dict["strand"] = entry.split(":")[1].split("(")[1][0]
        return fasta_header_dict


class EnsemblFastaBoy:

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.fasta_dict: Dict[str, Dict[str, str]] = dict()
        self.filter_list: List[Tuple[str, str]] = list()

    def get_fasta_dict(self) -> Dict[str, Dict[str, str]]:
        return self.fasta_dict

    def set_filter(self, key: str, value: str):
        self.filter_list.append((key, value))

    def parse_fasta(self):
        gene_id: str = ""
        protein_id: str = ""
        found_flag: bool = False
        for line in self:
            if line.startswith(">"):
                found_flag = False
                fasta_header_dict: Dict[str, str] = EnsemblFastaBoy.make_fasta_header_dict(line)
                if self.__apply_filter(fasta_header_dict):
                    found_flag = True
                    gene_id: str = fasta_header_dict["gene_id"]
                    protein_id: str = fasta_header_dict["protein_id"]
                    if gene_id not in self.fasta_dict.keys():
                        self.fasta_dict[gene_id] = dict()
                    self.fasta_dict[gene_id][protein_id] = ""
            elif found_flag:
                self.fasta_dict[gene_id][protein_id] += line.strip()

    def __apply_filter(self, fasta_header_dict: Dict[str, str]) -> bool:
        for key, value in self.filter_list:
            if fasta_header_dict[key] != value:
                return False
        return True

    @staticmethod
    def make_fasta_header_dict(fasta_header: str) -> Dict[str, str]:
        entry_list: List[str] = [entry for entry in fasta_header.split(" ") if EnsemblFastaBoy.discriminate_header(entry)]
        fasta_header_dict: Dict[str, str] = dict()
        for entry in entry_list:
            if entry.startswith(">"):
                fasta_header_dict["protein_id"] = entry[1:].split(".")[0]
            elif entry.startswith("gene:"):
                fasta_header_dict["gene_id"] = entry.split(":")[1].split(".")[0]
            elif entry.startswith("transcript:"):
                fasta_header_dict["transcript_id"] = entry.split(":")[1].split(".")[0]
            elif entry.startswith("gene_biotype:"):
                fasta_header_dict["gene_biotype"] = entry.split(":")[1]
            elif entry.startswith("transcript_biotype:"):
                fasta_header_dict["transcript_biotype"] = entry.split(":")[1]
        return fasta_header_dict

    @staticmethod
    def discriminate_header(entry: str):
        if entry.startswith(">"):
            return True
        elif entry.startswith("gene_biotype"):
            return True
        elif entry.startswith("transcript_biotype"):
            return True
        elif entry.startswith("gene"):
            return True
        elif entry.startswith("transcript"):
            return True
        else:
            return False

    def __iter__(self) -> Iterator[str]:
        with open(self.fasta_path, "r") as f:
            for line in f:
                yield line


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        action="store",
                        help="Path to a fasta file.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        action="store",
                        help="Path to json output file.")
    parser.add_argument("-m",
                        "--mode",
                        type=str,
                        action="store",
                        help="Either ensembl or transdecoder.")

    argument_dict: Dict[str, str] = vars(parser.parse_args())

    if argument_dict["mode"] == "ensembl":
        fasta_iterator: EnsemblFastaBoy = EnsemblFastaBoy(argument_dict["input"])
        fasta_iterator.set_filter("transcript_biotype", "protein_coding")
        fasta_iterator.set_filter("gene_biotype", "protein_coding")

        fasta_iterator.parse_fasta()

        with open(argument_dict["output"], "w") as f:
            json.dump(fasta_iterator.get_fasta_dict(), f, indent=4)
    elif argument_dict["mode"] == "transdecoder":
        fasta_iterator: TransDecoderFastaBoy = TransDecoderFastaBoy(argument_dict["input"])
        fasta_iterator.parse_fasta()

        with open(argument_dict["output"], "w") as f:
            json.dump(fasta_iterator.get_fasta_dict(), f, indent=4)


if __name__ == "__main__":
    main()