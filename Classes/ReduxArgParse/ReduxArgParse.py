#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  ReduxArgParse is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ReduxArgParse is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import argparse
from typing import List
from typing import Any
from typing import Dict


class ReduxArgParse:

    def __init__(self,
                 parse_tags: List[str],
                 parse_type: List[type],
                 parse_nargs: List[str],
                 parse_help: List[str]) -> None:
        self.parse_tags: List[str] = parse_tags
        self.parse_type: List[type] = parse_type
        self.parse_nargs: List[str] = parse_nargs  # ? zero or 1, * zero or more, + one or more
        self.parse_help: List[str] = parse_help
        self.parser: argparse.ArgumentParser = argparse.ArgumentParser()

    def generate_parser(self) -> None:
        for i in range(len(self.parse_tags)):
            shortcut_tag: str = "-" + self.parse_tags[i][2]
            self.parser.add_argument(shortcut_tag,
                                     self.parse_tags[i],
                                     type=self.parse_type[i],
                                     nargs=self.parse_nargs[i],
                                     help=self.parse_help[i])

    def execute(self) -> None:
        args: Any = self.parser.parse_args()
        return args

    def get_args(self) -> Dict[str, Any]:
        return vars(self.execute())


def main() -> None:
    argument_parser: ReduxArgParse = ReduxArgParse(["--input", "--output"], [str, str], ["?", "+"],  ["The input we want", "The output we want"])
    argument_parser.generate_parser()
    argument_parser.execute()
    print(argument_parser.get_args())


if __name__ == "__main__":
    main()
