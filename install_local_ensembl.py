#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

def make_local_ensembl_name(path, release_num, species, suffix):
        release_num = str(release_num)    
    
        server = "https://rest.ensembl.org"
        ext = "/info/assembly/" + species + "?"
        ####
        ## Getting current assembly accession
        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get assembly info for " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                assembly_accession = r.json()["assembly_accession"]

        ###########
        ### Getting current release NUM
        print("Getting url_name of species and default assembly name...")
        ext = "/info/genomes/assembly/" + assembly_accession + "?"

        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get genome info for assembly" + assembly_accession + " of " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                decoded = r.json()
                url_species = decoded["url_name"]
                assembly_name = decoded["assembly_default"]
        
        ensembl_path = path + url_species + "." + assembly_name + "." + release_num + suffix
        return ensembl_path

def get_species_info(species):
        server = "https://rest.ensembl.org"
        ext = "/info/assembly/" + species + "?"
        
        ####
        ## Getting current assembly accession
        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get assembly info for " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                assembly_accession = r.json()["assembly_accession"]

        ###########
        ### Getting current release NUM
        ext = "/info/genomes/assembly/" + assembly_accession + "?"

        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get genome info for assembly" + assembly_accession + " of " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                decoded = r.json()
                name = decoded["name"]
        return name

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
        return release_num

def install_local_ensembl(species, library_path):
        if not os.path.exists(library_path):
            os.makedirs(library_path)
    
        ftp_prefix = "http://ftp.ensembl.org/pub/release-"
        ftp_suffix = ".gtf.gz"

        print("Getting assembly accession id...")
        server = "https://rest.ensembl.org"
        ext = "/info/assembly/" + species + "?"
        
        ####
        ## Getting current assembly accession
        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get assembly info for " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                assembly_accession = r.json()["assembly_accession"]

        ###########
        ### Getting current release NUM
        print("Getting url_name of species and default assembly name...")
        ext = "/info/genomes/assembly/" + assembly_accession + "?"

        for x in range(3):
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

            if not r.ok and x < 2:
                continue
            elif x > 2:
                print("Could not get genome info for assembly" + assembly_accession + " of " + species + "...")
                r.raise_for_status()
                sys.exit()
            else:
                decoded = r.json()
                url_species = decoded["url_name"]
                url_infix_species = decoded["name"]
                assembly_name = decoded["assembly_default"]
        
        ftp_infix = "/gtf/" + url_infix_species+ "/"
        ###########
        ### Getting current release NUM
        release_num = get_release()    
        print("The current release is", release_num + "...")
        
        ### Assemble FTP address according to standard /release-{release_num}/{species}.{assembly_name}.{assembly_num}.gtf.gz
        file_name = url_species + "." + assembly_name + "." + release_num + ftp_suffix 
        ftp_address = ftp_prefix + release_num + ftp_infix + file_name 
        file_output_path = library_path + "/" + file_name
        print("Downloading", ftp_address, "to", file_output_path + "...")
        with closing(request.urlopen(ftp_address)) as r:
            with open(file_output_path, 'wb') as f:
                shutil.copyfileobj(r, f)
        
        ### Unpacking file
        print("Unpacking file...")
        with gzip.open(file_output_path, 'rb') as f_in:
            with open(file_output_path[:-3], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        os.remove(file_output_path) 
        print("Download complete. Local Ensembl Assembly can now be used.")

def main():
    #install_local_ensembl("mouse", 107, "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_Pipe/")
    print(make_local_ensembl_name("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/", 107, "homo_sapiens", ".gtf"))
if __name__ == "__main__":
    main()
        
                
                        
                    