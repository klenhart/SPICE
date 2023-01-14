#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  get_domain_importance is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  get_domain_importance is distributed in the hope that it will be useful,
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

from typing import Type
from typing import List


class Sequence():

    def __init__(self) -> None:
        pass

    def setID(self, seq : str) -> None:
        pass

    def setSequence(self, seq : str) -> None:
        pass

    def setAnnotation(self, annotation : List[str]) -> None:
        pass

    def setExpression(self, expression : float) -> None:
        pass

    def getID(self) -> str:
        pass

    def getSequence(self) -> str:
        pass

    def getAnnotation(self) -> List[str]:
        pass

    def getExpression(self) -> float:
        pass
