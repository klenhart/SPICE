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
Created on Wed Dec 21 10:03:16 2022

@author: Christian Bluemel
"""

import Sequence

import abc
from typing import Type
from typing import List

class SequenceContainer(metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def __init__(self) -> None:
        pass
    
    @abc.abstractmethod
    def setExpression(self, expression : float) -> None:
        pass
    
    @abc.abstractmethod
    def addSequence(self, seq : Type[Sequence]) -> None:
        pass
    
    @abc.abstractmethod
    def getSequences(self) -> List[Type[Sequence]]:
        pass
    
    @abc.abstractmethod
    def getExpression(self) -> float:
        pass
        
    