#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
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
    optional.add_argument("-c", "--genesPerFile", default=100, type=int, required=False,
                          help="Number of genes in one file. Higher number result in less file with bigger size. ")
    args = parser.parse_args()
    main(args.inPath, args.outPath, args.genesPerFile)


def main(inpath, outpath, filesize):
    arc = read_input(inpath)
    index = 0
    count = 0
    output = {}
    indexout = {'genes': {}}
    for gene in arc:
        output[gene] = arc[gene]
        indexout['genes'][gene] = index
        count += 1
        if count > filesize:
            name = str(index).rjust(9, '0')
            save2json(output, name, outpath)
            output = {}
            index += 1
            count = 0
    name = str(index).rjust(9, '0')
    save2json(output, name, outpath)
    indexout['#files'] = index
    save2json(indexout, 'index', outpath)

    
def save2json(dict2save, name, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    out = open(directory + '/' + name + '.json', 'w')
    out.write(jsonOut)
    out.close()


def read_input(inpath):    # reads input json that contains the isoform annotations and restructures the data in a gene centric fashion
    fa_map = {}
    with open(inpath, 'r') as infile:
        features = json.loads(infile.read())['feature']
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

