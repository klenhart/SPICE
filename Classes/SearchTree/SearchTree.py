# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  SearchTree is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SearchTree is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict, Tuple, List, Any

from Classes.SearchTree.AbstractSearchTreeEntry import AbstractSearchTreeEntry


class SearchTree(AbstractSearchTreeEntry):

    def __init__(self, key: str, leaf_flag: bool = False) -> None:
        self.leaf_flag: bool = leaf_flag
        self.key: str = key
        if self.leaf_flag:
            self.entry: AbstractSearchTreeEntry
        else:
            self.children: Dict[str, SearchTree] = dict()

    def get_id(self) -> str:
        return self.key

    def pass_down_entry(self, entry_pair: Tuple[str, AbstractSearchTreeEntry]) -> None:
        key: str = entry_pair[0][0]
        if len(entry_pair[0]) == 1:
            if key not in self.children.keys():
                self.children[key] = SearchTree(key, True)
                self.children[key].insert_entry(entry_pair[1])
        else:
            entry_pair: Tuple[str, AbstractSearchTreeEntry] = (entry_pair[0][1:], entry_pair[1])
            if entry_pair[0] not in self.children.keys():
                self.children[key] = SearchTree(key)
                self.children[key].pass_down_entry(entry_pair)

    def insert_entry(self, entry: AbstractSearchTreeEntry) -> None:
        if self.leaf_flag:
            self.entry: AbstractSearchTreeEntry = entry
        else:
            entry_pair: Tuple[str, AbstractSearchTreeEntry] = (entry.get_id(), entry)
            self.pass_down_entry(entry_pair)

    def get_entry(self) -> List[AbstractSearchTreeEntry]:
        if self.leaf_flag:
            return [self.entry]
        else:
            output: List[AbstractSearchTreeEntry] = list()
            for key in self.children.keys():
                output += self.children[key].get_entry()

    def find(self, find_id: str, depth: int = 0) -> AbstractSearchTreeEntry:
        if self.leaf_flag:
            return self.entry
        else:
            try:
                key: str = find_id[0]
                find_id = find_id[1:]
                return self.children[key].find(find_id, depth+1)
            except KeyError:
                if depth == 0:
                    print(find_id, "is not a key in this SearchTree.")

    def flatten(self) -> List[AbstractSearchTreeEntry]:
        if self.leaf_flag:
            return [self.entry]
        else:
            flattened_tree: List[AbstractSearchTreeEntry] = []
            for key in self.children.keys():
                flattened_tree += self.children[key].flatten()

    def to_dict(self) -> Dict[str, Any]:
        all_entries: List[AbstractSearchTreeEntry] = self.flatten()
        output_dict: Dict[str, Dict[str, Any]] = dict()
        for entry in all_entries:
            output_dict[entry.get_id()] = entry.to_dict()

    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        self.children: Dict[str, SearchTree] = dict()
        for key in input_dict.keys():
            input_dict[key]

    def add_entry(self, entry_type: str, entry: Any) -> None:  # TODO Implement add_entry of SearchTree
        pass
