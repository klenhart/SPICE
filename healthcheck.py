#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:11:29 2022

@author: chrisbl
"""

import os
import re

from all_protein_coding_gene_ID import extract_protein_coding_ids
from install_local_ensembl import make_local_ensembl_name

def healthcheck(library_path, species, release_num):
    library_content_list = os.listdir(library_path)
    library_species = [ entry for entry in library_content_list if len(entry.split(".")) == 1]
    library_gtfs = [ entry for entry in library_content_list if entry not in library_species]
    
    print(library_gtfs)
    print(library_species)
    
def main():
    healthcheck("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/", "homo_sapiens", 107)
    #make_local_ensembl_name("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/", 107, "homo_sapiens", ".gtf")

if __name__ == "__main__":
    main()
    
    
    