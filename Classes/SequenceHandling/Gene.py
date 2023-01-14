#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Gene is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Gene is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

"""
Created on Wed Dec 21 10:03:16 2022

@author: Christian Bluemel
"""

from Classes.SequenceHandling import Protein


from typing import Type
from typing import List


class Gene:

    def __init__(self) -> None:
        self.expression_value: float
        self.sequences: Type[Protein] = list()
        self.expression_value: float = 0

    def set_expression(self, expression: float) -> None:
        pass

    def set_expression_value(self, expression: float) -> None:
        """

        :type expression: float
        """
        self.expression_value = expression

    def add_sequence(self, seq: Type[Protein]) -> None:
        pass

    def get_sequences(self) -> List[Type[Protein]]:
        pass

    def get_expression_value(self) -> float:
        return self.expression_value
