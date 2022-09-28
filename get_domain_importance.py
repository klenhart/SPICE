#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of main.
#
#  XXXX is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  XXXX is distributed in the hope that it will be useful,
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
    parser = argparse.ArgumentParser(epilog="This script uses the annotation file of the different protein isoforms of each gene in a proteome to generate a dictionary of features that have presence/absence changes due to AS.")
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--inPath", default='.', type=str, required=True,
                          help="path to input json")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory. Name will be based on input file.")
    args = parser.parse_args()
    main(args.inPath, args.outPath)


def main(inpath, outpath):
    genes = read_input(inpath)
    output = get_fdict(genes)
    name = '.'.join(inpath.split('/')[-1].split('.')[0:-1])
    save2json(output, name, outpath)
    
    
def save2json(dict2save, name, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    jsonOut = json.dumps(dict2save, ensure_ascii=False)
    out = open(directory + '/' + name + '_important_features.json', 'w')
    out.write(jsonOut)
    out.close()
    

def get_fdict(genes):    # generates dict {feature: [gene1, gene2...]}
    f_dict = {}
    for gene in genes:
        f_list = []
        for isoform in genes[gene]['isoforms']:
            for feature in genes[gene]['isoforms'][isoform]['types']:
                if not feature in f_list:
                    if genes[gene]['isoforms'][isoform]['types'][feature] < genes[gene]['max'][feature]:
                        f_list.append(feature)
            for clan in genes[gene]['isoforms'][isoform]['clans']:
                if not clan in f_list:
                    if genes[gene]['isoforms'][isoform]['clans'][clan] < genes[gene]['max'][clan]:
                        f_list.append(clan)
        for feature in f_list:
            if not feature in f_dict:
                f_dict[feature] = []
            f_dict[feature].append(gene)
    return f_dict


def read_input(inpath):    # reads input json that contains the isoform annotations and restructures the data in a gene centric fashion
    with open(inpath, 'r') as infile:
        in_dict = json.loads(infile.read())
        features, clans = in_dict['feature'], in_dict['clan']
    genes = {}
    for protein in features:
        gid, pid, tid = protein.split('|')    # isoform id structure: geneid|proteinid|taxid
        if not gid in genes:
            genes[gid] = {'max': {}, 'isoforms': {}}
        genes[gid]['isoforms'][pid] = {'types': {}, 'clans': {}}
        for tool in features[protein]:
            if not tool=="length":
                for ftype in features[protein][tool]:
                    count = len(features[protein][tool][ftype]["instance"])
                    if not ftype in genes[gid]['max']:
                        genes[gid]['max'][ftype] = count
                    elif count > genes[gid]['max'][ftype]:
                        genes[gid]['max'][ftype] = count
                    genes[gid]['isoforms'][pid]['types'][ftype] = count
                    if ftype in clans:
                        clan = clans[ftype]
                        if not clan in genes[gid]['isoforms'][pid]['clans']:
                            genes[gid]['isoforms'][pid]['clans'][clan] = count
                        else:
                            genes[gid]['isoforms'][pid]['clans'][clan] += count
        for clan in genes[gid]['isoforms'][pid]['clans']:
            count = genes[gid]['isoforms'][pid]['clans'][clan]
            if not clan in genes[gid]['max']:
                genes[gid]['max'][clan] = count
            elif count > genes[gid]['max'][clan]:
                genes[gid]['max'][clan] = count
    return genes
                        


if __name__ == '__main__':
    option_parse()

