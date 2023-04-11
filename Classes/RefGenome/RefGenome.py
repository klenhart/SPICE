#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  RefGenome is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RefGenome is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Iterator, List


class RefGenome:

    def __init__(self, genome_path: str):
        self.path = genome_path

    def __iter__(self) -> Iterator[str]:
        with open(self.path, "r") as f:
            for line in f:
                yield line

    def get_sequence(self, chromosome: str, start: int, end: int) -> str:
        sequence: str = ""
        start_found: bool = False
        end_found: bool = False
        found_flag: bool = False
        for line in self:
            if end_found:
                break
            elif start_found:
                pass
            if found_flag:
                if start - len(line) < 0:
                    pass

            else:
                if line.startswith(">"):
                    split_line: List[str] = line[1:].split(" ")
                    if chromosome in split_line:
                        found_flag: bool = True
        return sequence


def main():
    path: str = "C:\\Users\\chris\\Downloads\\GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    ref_genome: RefGenome = RefGenome(path)
    print(ref_genome.get_sequence("chr1", 1598010, 1599476))


if __name__ == "__main__":
    main()
