#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import argparse
import json
import os
from pathlib import Path


def option_parse():
    parser = argparse.ArgumentParser(epilog="This script restructures the annotation file format into a mapped version where each feature instance of a protein gets an id. This is to make saving linearized architectures less data heavy.")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--inPath", default='.', type=str, required=True,
                          help="path to input json")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory. Name will be based on input file.")
    args = parser.parse_args()
    main(args.inPath, args.outPath)


def main(inpath, outpath):
    output = read_input(inpath)
    name = '.'.join(inpath.split('/')[-1].split('.')[0:-1])
    save2json(output, name, outpath)


    
def save2json(dict2save, name, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    out = open(directory + '/' + name + '_map.json', 'w')
    out.write(jsonOut)
    out.close()


def read_input(inpath):    # reads input json that contains the isoform annotations and restructures the data in a gene centric fashion
    fa_map = {}
    with open(inpath, 'r') as infile:
        in_dict = json.loads(infile.read())
        features = in_dict['feature']
        for protid in features:
            gid, pid, tid = protid.split('|')
            if not gid in fa_map:
                fa_map[gid] = {}
            i = 0
            fa_map[gid][pid] = {'fmap': {}}
            for tool in features[protid]:
                if tool == 'length':
                    fa_map[gid][pid]['length'] = features[protid]['length']
                else:
                    for feature in features[protid][tool]:
                        for instance in features[protid][tool][feature]['instance']:
                            fa_map[gid][pid]['fmap'][i] = (feature, instance[0], instance[1])
                            i += 1
    return fa_map
                        


if __name__ == '__main__':
    option_parse()

