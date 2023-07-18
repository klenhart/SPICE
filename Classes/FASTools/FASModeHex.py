#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  FASModeHex is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASModeHex is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import List, Tuple, Dict


class FASModeHex:

    hex_dict: Dict[str, str] = {"0000": "0", "0001": "1", "0010": "2", "0011": "3", "0100": "4",
                                "0101": "5", "0110": "6", "0111": "7", "1000": "8", "1001": "9",
                                "1010": "a", "1100": "b", "1101": "c", "1110": "d", "1111": "e"}

    def __init__(self,
                 mode_byte_count: int = 3,
                 mode_tuple: Tuple[str] = ("#linearized", "Pfam", "#normal",
                                           "fLPS", "COILS2", "SEG",
                                           "SignalP", "TMHMM", "#checked")):
        self.mode_byte_count: int = mode_byte_count
        self.mode_hex: str = "0" * mode_byte_count
        self.mode_bits: List[bool] = [False] * mode_byte_count * 4
        self.blocked_bits: int = mode_byte_count*4-len(mode_tuple)
        self.mode_tuple: Tuple[str] = mode_tuple

    def get_mode_hex(self) -> str:
        return self.mode_hex

    def __update_mode_hex__(self) -> None:
        self.mode_hex = FASModeHex.bool_list_to_hex(self.mode_bits, self.mode_byte_count)

    def activate_all(self):
        self.activate_modes(list(self.mode_tuple))

    def activate_modes(self, mode_name_list: List[str]) -> None:
        for mode_name in mode_name_list:
            self.activate_mode(mode_name)

    def activate_mode(self, mode_name: str) -> None:
        try:
            index: int = self.mode_tuple.index(mode_name) + self.blocked_bits
            self.mode_bits[index] = True
            self.__update_mode_hex__()
        except ValueError:
            print("Given mode name", mode_name, "not in mode list: ", self.mode_tuple)

    def deactivate_modes(self, mode_name_list: List[str]) -> None:
        for mode_name in mode_name_list:
            self.deactivate_mode(mode_name)

    def deactivate_mode(self, mode_name: str) -> None:
        try:
            index: int = self.mode_tuple.index(mode_name) + self.blocked_bits
            self.mode_bits[index] = False
            self.__update_mode_hex__()
        except ValueError:
            print("Given mode name", mode_name, "not in mode list: ", self.mode_tuple)

    def __str__(self) -> str:
        output_list: List[str] = list()
        for i, bit in enumerate(self.mode_bits):
            if bit:
                output_list.append(self.mode_tuple[i-self.blocked_bits])
        return "\n".join(output_list)

    @staticmethod
    def bool_list_to_hex(bool_list: List[bool], byte_limit: int) -> str:
        byte_count: int = 0
        hex_string: str = ""
        while byte_count < byte_limit:
            byte_start: int = byte_count * 4
            byte_end: int = byte_start + 4
            bit_string: str = ""
            for entry in bool_list[byte_start:byte_end]:
                if entry:
                    bit_string += "1"
                else:
                    bit_string += "0"
            hex_string += FASModeHex.hex_dict[bit_string]
            byte_count += 1
        return hex_string


def main():
    fas_mode_hex: FASModeHex = FASModeHex()
    fas_mode_hex.activate_modes(["#linearized", "#normal", "TMHMM", "#checked"])
    print(fas_mode_hex)
    print(fas_mode_hex.get_mode_hex())


if __name__ == "__main__":
    main()
