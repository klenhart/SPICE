#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################################################
# Copyright (C) 2022 Christian, Blümel, Julian Dosch
#
# This file is part of grand-trumpet.
#
#  grand-trumpet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  grand-trumpet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with grand-trumpet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
"""
Created on Thu Oct  6 10:44:50 2022

@author: chrisbl
"""

from valves.fas_utility import longest_common_prefix

class SearchTree:
    def __init__(self, tuple_list):
        words = [ entry[0] for entry in tuple_list ]
        values = [ entry[1] for entry in tuple_list ]

        self.flag_leaf = len(words[0]) == 1
        self.prefix = ""
        
        if not self.flag_leaf:
            self.prefix = longest_common_prefix(words)

        self.node_dict = dict()
                
        words = [ entry[len(self.prefix):] for entry in words ]
        tuple_list = zip(words, values)
        
        

        if self.flag_leaf:
            for letter, value in tuple_list:
                self.node_dict[letter] = value 
        else:  
            starting_letters = set([ entry[0] for entry in words ])
            for letter in starting_letters:
                self.node_dict[letter] = SearchTree([ (entry[0][1:], entry[1]) for entry in tuple_list if entry[0][0] == letter ])
    
    def search_word(self, word):
        if self.flag_leaf:
            return self.node_dict[word]
        else:
            if word.startswith(self.prefix):
                word = word[len(self.prefix):]
                letter = word[0]
                word = word[1:]
                return self.node_dict[letter].search_word(word)
            else:
                raise "The searched word is not part of the tree."
            
        
        