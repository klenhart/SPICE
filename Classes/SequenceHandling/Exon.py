#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Exon is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Exon is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

"""
Created on Wed Dec 21 10:33:00 2022

@author: chris
"""

from typing import List


class Exon:

    def __init__(self) -> None:
        self.id: str = ""
        self.begin: int = 0
        self.end: int = 0

    def set_id(self, id_exon: str) -> None:
        """

        :type id_exon: str
        """
        self.id = id_exon

    def set_begin(self, begin: int) -> None:
        self.begin = begin

    def set_end(self, end: int) -> None:
        self.end = end

    def get_id(self) -> str:
        return self.id

    def get_begin(self) -> int:
        return self.begin

    def get_end(self) -> int:
        return self.end
