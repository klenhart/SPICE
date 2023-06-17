#!/bin/env python


#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  NaiveTranslator is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  NaiveTranslator is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


from typing import List, Dict

STANDARD: Dict[str, str] = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S",
                            "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C",
                            "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P",
                            "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
                            "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
                            "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V",
                            "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
                            "GGG": "G"}


INVERT: Dict[str, str] = {"G": "C", "C": "G", "A": "T", "T": "A"}


def seq_inverter(sequence: str):
    new_seq: str = ""
    for char in sequence:
        new_seq += INVERT[char]
    return new_seq[::-1]


class NaiveTranslator:

    def __init__(self, sequence: str,
                 start_codons=None,
                 stop_codons=None,
                 translation_dict=None):

        if start_codons is None:
            start_codons = ["ATG"]
        if stop_codons is None:
            stop_codons = ["TAA", "TGA", "TAG"]
        if translation_dict is None:
            translation_dict = STANDARD

        self.start_flag: bool = True
        self.stop_flag: bool = False
        self.non_coding_flag: bool = False

        self.translation_dict = translation_dict

        self.start_codons = start_codons
        self.stop_codons = stop_codons
        self.raw_sequence = sequence
        self.orf = ""
        self.peptide = ""

        self.start_index: int = -1
        self.stop_index: int = -1

        self.__identify_frame()

        if self.stop_flag and self.start_flag:
            self.__translate()
        else:
            print("NO ORF COULD BE PREDICTED.")
            self.non_coding_flag = True

    def get_orf(self):
        return self.orf

    def get_peptide(self):
        return self.peptide

    def get_orf_indices(self):
        return self.start_index, self.stop_index

    def get_raw_sequence(self):
        return self.raw_sequence

    def is_coding(self):
        return not self.non_coding_flag

    def __translate(self):
        codon_count: int = len(self.orf) // 3
        for codon in range(codon_count):
            i: int = codon * 3
            self.peptide += self.translation_dict[self.orf[i:i+3]]

    def __identify_frame(self) -> None:
        start_index_list: List[int] = list()
        for codon in self.start_codons:
            try:
                self.raw_sequence.index(codon)
                index: int = self.raw_sequence.index(codon)
                start_index_list.append(index)
            except ValueError:
                print(codon, "not found.")

        if len(start_index_list) == 0:
            self.start_flag = False
        else:
            self.start_index = min(start_index_list)

            subsequence: str = self.raw_sequence[self.start_index:]
            for i, letter in enumerate(subsequence):
                if i % 3 == 0 and subsequence[i:i+3] in self.stop_codons:
                    self.stop_index = i + self.start_index
                    self.stop_flag = True
                    self.orf = self.raw_sequence[self.start_index: self.stop_index]
                    break


def main():
    print("Hello World!")

if __name__ == "__main__":
    main()
