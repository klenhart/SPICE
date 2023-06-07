#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import argparse
from typing import List
from typing import Any
from typing import Dict


class ReduxArgParse:

    def __init__(self,
                 parse_tags: List[str],
                 parse_types: List[type],
                 parse_actions: List[str],
                 parse_nargs: List[Any],
                 parse_help: List[str]) -> None:
        self.parse_tags: List[str] = parse_tags
        self.parse_types: List[type] = parse_types
        self.parse_actions: List[str] = parse_actions  # store, append, store_true, store_false, count, help, version
        self.parse_nargs: List[str] = parse_nargs  # ? 0 or 1, * 0 or more, + 1 or more, <number> that exact amount
        self.parse_help: List[str] = parse_help
        self.parser: argparse.ArgumentParser = argparse.ArgumentParser()

    def generate_parser(self) -> None:
        for i in range(len(self.parse_tags)):
            shortcut_tag: str = "-" + self.parse_tags[i][2]
            if self.parse_actions[i] in ["store_true", "store_false"]:
                self.parser.add_argument(shortcut_tag,
                                         self.parse_tags[i],
                                         action=self.parse_actions[i],
                                         help=self.parse_help[i])
            else:
                self.parser.add_argument(shortcut_tag,
                                         self.parse_tags[i],
                                         type=self.parse_types[i],
                                         action=self.parse_actions[i],
                                         nargs=self.parse_nargs[i],
                                         help=self.parse_help[i])

    def execute(self):
        args: Any = self.parser.parse_args()
        return args

    def get_args(self) -> Dict[str, Any]:
        return vars(self.execute())


def main() -> None:
    argument_parser: ReduxArgParse = ReduxArgParse(["--flag", "--list", "--integer", "--multilist"],
                                                   [None, str, int, str],
                                                   ["store_true", "store", "store", "append"],
                                                   [None, 1, "?", "*"],
                                                   ["help1", "help2", "help3", "help4"])
    argument_parser.generate_parser()
    argument_parser.execute()
    print(argument_parser.get_args())


if __name__ == "__main__":
    main()
