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
    required.add_argument("-f", "--forwardPath", default='.', type=str, required=True,
                          help="path to _forward.domains")
    required.add_argument("-r", "--reversePath", default='.', type=str, required=True,
                          help="path to _reverse.domains")
    required.add_argument("-m", "--mapPath", default='.', type=str, required=True,
                          help="path to feature mapping json")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory.")
    args = parser.parse_args()
    main(args.forwardPath, args.reversePath, args.mapPath, args.outPath)


def main(forwardpath, reversepath, mappath, outpath):
    with open(mappath + '/index.json', 'r') as infile:
        files = json.loads(infile.read())['#files']
    read_input((forwardpath, reversepath), mappath, files)


    
def save2json(dict2save, name, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    out = open(directory + '/' + name + '.json', 'w')
    out.write(jsonOut)
    out.close()


def read_input(inpaths, mappath, files):
    for i in [0, 1]:
        with open(inpaths[i], 'r') as infile:
            lines = infile.readlines()
        index2 = {}
        y = 0
        current = None
        for line in lines:
            if line == '\n' or line.startswith('#'):
                y += 1
                continue
            gid = line.split('|')[0]
            if not gid in index2:
                index2[gid] = [y]
            if current and not gid == current:
                index2[current].append(y)
                current = gid
            if not current:
                current = gid
            y += 1
        index2[current].append(y)
        for index in range(files+1):
            with open(mappath + '/' + str(index).rjust(9, '0') + '.json', 'r') as infile2:     
                mapfile = json.loads(infile2.read())
            if i == 0:
                lin = {}
            else:
                with open(mappath + '/' + str(index).rjust(9, '0') + '_paths.json', 'r') as pathfile:
                    lin = json.loads(pathfile.read())
            for gene in mapfile:
                if gene in index2:
                    for z in range(index2[gene][0], index2[gene][1]):
                        line = lines[z]
                        cells = line.rstrip('\n').split('\t')
                        if line == '\n' or line.startswith('#'):
                            continue
                        if i == 0:
                            p1, p2 = cells[0].split('#')
                        else:
                            p2, p1 = cells[0].split('#')
                        if p1 == p2:
                            continue
                        gid, p1, tax = p1.split('|')
                        p2 = p2.split('|')[1]
                        if gid not in lin:
                            lin[gid] = {}
                        pkey = '@'.join((p1, p2))
                        if pkey not in lin[gid]:
                            lin[gid][pkey] = [[], []]
                        rp = cells[1].split('|')[1]
                        if cells[7] == 'Y':
                            fid = None
                            x = 0
                            while fid == None:
                                if mapfile[gid][rp]['fmap'][str(x)][0] == cells[3] and mapfile[gid][rp]['fmap'][str(x)][1] == int(cells[4]) and mapfile[gid][rp]['fmap'][str(x)][2] == int(cells[5]):
                                    fid = str(x)
                                x += 1
                            if rp == p1:
                                lin[gid][pkey][0].append(fid)
                            elif rp == p2:
                                lin[gid][pkey][1].append(fid)
            save2json(lin, str(index).rjust(9, '0') + '_paths', mappath)



if __name__ == '__main__':
    option_parse()

