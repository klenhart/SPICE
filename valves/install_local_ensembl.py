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
Created on Tue Jul 19 13:05:13 2022

@author: chrisbl
"""

import requests
import os
import sys
import shutil
import gzip
import urllib.request as request
from contextlib import closing

from valves.library_class import Library

TMHMM_SIGNALP = """#linearized
#normal
TMHMM
SignalP
#checked
"""

LCR = """#linearized
#normal
flPS
SEG
#checked
"""



def ping_ensembl():
    import requests, sys
    
    
    server = "https://rest.ensembl.org"
    ext = "/info/ping?"
    
    for i in range(3):
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
        if r.ok:
            break
        elif i >= 2:
            print("Ensembl is currently down. Can't download sequences. Aborting...")
            r.raise_for_status()
            sys.exit() 
    decoded = r.json()
    return bool(decoded["ping"])

def make_local_ensembl_name(path, release_num, species, suffix, assembly_name, url_species):
        release_num = str(release_num)        
        ensembl_path = path + url_species + "." + assembly_name + "." + release_num + suffix
        return ensembl_path

def get_species_info(species):
    server = "https://rest.ensembl.org"
    ext = "/info/genomes/" + species + "?"
     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      print("Could not find species!")
      r.raise_for_status()
      sys.exit()
     
    decoded = r.json()
    url_name = decoded["url_name"]
    name = decoded["name"]
    assembly_default = decoded["assembly_default"]
    return name, url_name, assembly_default


def get_release():
        print("Checking release number...")
        server = "https://rest.ensembl.org"
        ext = "/info/data/?"
 
        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not current release number...")
                r.raise_for_status()
                sys.exit()
            else:
                release_num = max(r.json()["releases"])
        return str(release_num)

def make_rootpath(library_path, species, assembly_num):
    return library_path + species + "/release-" + str(assembly_num) + "/"

def get_taxon_id(species): 
    server = "https://rest.ensembl.org"
    ext = "/info/genomes/taxonomy/" + species +"?"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    decoded = r.json()
    return decoded[0]["species_taxonomy_id"]

def make_folders_and_files(root_path):
    if not os.path.exists(root_path):
        os.makedirs(root_path)
    
    tsv_buffer_path = root_path + "temp/tsv_buffer/"
    if not os.path.exists(tsv_buffer_path):
        os.makedirs(tsv_buffer_path)
    
    slurm_path = root_path + "temp/SLURM/"
    if not os.path.exists(slurm_path):
        os.makedirs(slurm_path)
    
    fas_buffer_path = root_path + "temp/FAS_buffer/"
    if not os.path.exists(fas_buffer_path):
        os.makedirs(fas_buffer_path)
    
    fas_buffer_tmhmm_path = root_path + "temp/FAS_buffer_tmhmm/"
    if not os.path.exists(fas_buffer_tmhmm_path):
        os.makedirs(fas_buffer_tmhmm_path)
    
    fas_buffer_lcr_path = root_path + "temp/FAS_buffer_lcr/"
    if not os.path.exists(fas_buffer_lcr_path):
        os.makedirs(fas_buffer_lcr_path)
        
    annotation_path = root_path + "src/annotation/"
    if not os.path.exists(annotation_path):
        os.makedirs(annotation_path)
    
    isoforms_path = root_path + "src/isoforms.fasta"
    if not os.path.isfile(isoforms_path):
        with open(isoforms_path, "w") as fp:
            pass
    
    distance_master_path = root_path + "src/fas_all.phyloprofile"
    if not os.path.isfile(distance_master_path):
        with open(distance_master_path, "w") as fp:
            fp.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    
    distance_master_lcr_path = root_path + "src/fas_lcr.phyloprofile"
    if not os.path.isfile(distance_master_lcr_path):
        with open(distance_master_lcr_path, "w") as fp:
            fp.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
            
    distance_master_tmhmm_path = root_path + "src/fas_tmhmm.phyloprofile"
    if not os.path.isfile(distance_master_tmhmm_path):
        with open(distance_master_tmhmm_path, "w") as fp:
            fp.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")
    
    phyloprofile_ids_path = root_path + "src/phyloprofile_ids.tsv"
    if not os.path.isfile(phyloprofile_ids_path):
        with open(phyloprofile_ids_path, "w") as fp:
            pass
    
    gene_ids_path = root_path + "src/gene_ids.txt"
    if not os.path.isfile(gene_ids_path):
        with open(gene_ids_path, "w") as fp:
            pass
        
    protein_coding_ids_path = root_path + "src/protein_coding_genes.tsv"
    if not os.path.isfile(protein_coding_ids_path):
        with open(protein_coding_ids_path, "w") as fp:
            pass
    
    canonical_path = root_path + "src/canonical_transcripts.fasta"
    if not os.path.isfile(canonical_path):
        with open(canonical_path, "w") as fp:
            pass
    
    tmhmm_path = root_path + "src/featuretypes_tmhmm_SignalP.txt"
    if not os.path.isfile(tmhmm_path):
        with open(tmhmm_path, "w") as fp:
            fp.write(TMHMM_SIGNALP)
    
    lcr_path = root_path + "src/featuretypes_lcr.txt"
    if not os.path.isfile(lcr_path):
        with open(lcr_path, "w") as fp:
            fp.write(LCR)
    
    result_path = root_path + "result/"
    if not os.path.exists(result_path):
        os.makedirs(result_path)
        
    movement_path = root_path + "movement/"
    if not os.path.exists(movement_path):
        os.makedirs(movement_path)
        
    expression_path = root_path + "expression/"
    if not os.path.exists(expression_path):
        os.makedirs(expression_path)

    main_comparison_path = root_path + "main_comparison/"
    if not os.path.exists(main_comparison_path):
        os.makedirs(main_comparison_path)

    result_config_path = root_path + "result_config.tsv"
    if not os.path.isfile(result_config_path):
        with open(lcr_path, "w") as fp:
            fp.write("condition\tsamples\tFAS_modes\tLibrary\texpr_normalization\n")
    
    path_dict = dict()
    path_dict["fas_all_path"] = distance_master_path
    path_dict["fas_lcr_path"] = distance_master_lcr_path
    path_dict["fas_tmhmm_path"] = distance_master_tmhmm_path
    path_dict["tmhmm_path"] = tmhmm_path
    path_dict["lcr_path"] = lcr_path
    path_dict["root_path"] = root_path
    path_dict["gene_ids_path"] = gene_ids_path
    path_dict["isoforms_path"] = isoforms_path
    path_dict["protein_coding_ids_path"] = protein_coding_ids_path
    path_dict["phyloprofile_ids_path"] = phyloprofile_ids_path
    path_dict["slurm_path"] = slurm_path
    path_dict["tsv_buffer_path"] = tsv_buffer_path
    path_dict["fas_buffer_path"] = fas_buffer_path
    path_dict["fas_buffer_tmhmm_path"] = fas_buffer_tmhmm_path
    path_dict["fas_buffer_lcr_path"] = fas_buffer_lcr_path
    path_dict["annotation_path"] = annotation_path
    path_dict["canonical_path"] = canonical_path
    path_dict["result"] = result_path
    path_dict["result_config"] = result_config_path
    path_dict["main_comparison"] = main_comparison_path
    path_dict["expression"] = expression_path
    path_dict["movement"] = movement_path
    
    return path_dict

def install_local_ensembl(species, output_dir):
        ping_ensembl()
        library_path = output_dir + "/FAS_library/"
        if not os.path.exists(library_path):
            os.makedirs(library_path)
            os.makedirs(library_path + "/annoTools/")
        release_num = get_release()
        species, url_name, assembly_default = get_species_info(species)
        taxon_id = get_taxon_id(species)
        root_path = make_rootpath(library_path, species, release_num)
        
        path_dict = make_folders_and_files(root_path)
        pairings_tsv_json_path = root_path + "src/pairings_tsv.json"

        ftp_prefix = "http://ftp.ensembl.org/pub/release-"
        ftp_suffix = ".gtf.gz"        
        ftp_infix = "/gtf/" + species + "/"
        
        print("The current release is", release_num + "...")
        
        ### Assemble FTP address according to standard /release-{release_num}/{species}.{assembly_dafault}.{assembly_num}.gtf.gz
        file_name = url_name + "." + assembly_default + "." + release_num + ftp_suffix 
        ftp_address = ftp_prefix + release_num + ftp_infix + file_name 
        local_assembly_path = root_path + "src/" + file_name
        print("Downloading", ftp_address, "to", local_assembly_path + "...")
        with closing(request.urlopen(ftp_address)) as r:
            with open(local_assembly_path, 'wb') as f:
                shutil.copyfileobj(r, f)
        
        ### Unpacking file
        with gzip.open(local_assembly_path, 'rb') as f_in:
            with open(local_assembly_path[:-3], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(local_assembly_path)
        
        config_path = root_path + "config.tsv"
        fas_lib = Library(None, True)
        fas_lib.set_config("self_path", config_path)
        fas_lib.set_config("pairings_tsv_json_path", pairings_tsv_json_path)

        fas_lib.set_config("species", species)
        fas_lib.set_config("release_num", release_num)
        fas_lib.set_config("url_name", url_name)
        fas_lib.set_config("assembly_default", assembly_default)
        fas_lib.set_config("taxon_id", taxon_id)
    
        for key in path_dict.keys():
            fas_lib.set_config(key, path_dict[key])
        
        fas_lib.save_config()
        
        print("Library scaffold constructed. Now run main.py using the option -c", config_path)


def main():
    pass

if __name__ == "__main__":
    main()
        
                
                        
                    