#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:34:52 2022

@author: chrisbl
"""

"""
Move all minor tools from all modules to this.
"""

def start_stop_range(length, n):
    """
    Splits length into equal chunks of size n.
    
    Parameters
    ----------
    length : int
        length to be split apart.
    n : int
        size of chunks.

    Yields
    ------
    generator of (int, int)
        Splits length into equal chunks of size n.
    """
    for i in range(0, length, n):
        yield (i, min(i+n-1, length))

def get_name(path):
    start = path.index("polygonFAS_") + 11
    end = path.index(".tsv")
    return path[start:end]