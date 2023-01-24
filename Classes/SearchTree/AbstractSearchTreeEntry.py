# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  AbstractSearchTreeNode is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AbstractSearchTreeNode is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from abc import ABC
from abc import abstractmethod
from typing import Dict, Any


class AbstractSearchTreeEntry(ABC):

    @abstractmethod
    def get_id(self) -> str:
        pass

    @abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        pass

    @abstractmethod
    def from_dict(self, input_dict: Dict[str, Any]) -> None:
        pass

    @abstractmethod
    def add_entry(self, entry_type: str, entry: Any) -> None:
        pass
