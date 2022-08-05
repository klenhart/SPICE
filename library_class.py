# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 09:30:21 2022

@author: chris
"""

import os
import sys

class Library:
    def __init__(self, config_path, flag_init):
        #config dict
        self.config = dict()
        
        #stats
        self.config["species"] = None
        self.config["release_num"] = None
        self.config["url_name"] = None
        self.config["assembly_default"] = None
        self.config["taxon_id"] = None
        
        #paths
        self.config["library_path"] = None
        self.config["root_path"] = None
        self.config["self_path"] = None
        self.config["pairings_tsv_json_path"] = None
        self.config["gene_ids_path"] = None
        self.config["isoforms_path"] = None
        self.config["phyloprofile_ids_path"] = None
        self.config["protein_coding_ids_path"] = None
        self.config["local_assembly_path"] = None
        self.config["slurm_path"] = None
        self.config["tsv_buffer_path"] = None
        self.config["fas_buffer_path"] = None
        self.config["annotation_path"] = None
        
        #counts
        self.config["acquired_seq_count"] = "0"
        self.config["total_seq_count"] = "0"
        self.config["max_isoform_num"] = "0"
        self.config["gene_count"] = "0"
        
        #flags
        self.config["flag_protein_coding_genes"] = "False"
        self.config["flag_sequence_collection"] = "False"
        self.config["flag_gene_ids_collection"] = "False"
        self.config["flag_made_pairings"] = "False"

        if not flag_init:
            print("Reached it!")
            if os.path.isfile(config_path):
                with open(config_path, "r") as f:
                    config = f.read()
                    config = config.split("\n")
                    config = [ entry.split("\t") for entry in config if len(entry) > 0 ]
                    config = [ (category, entry if entry != "." else None) for category, entry in config ]
                self.extract_config(config)            

    def increment_acquired_seq_count(self):
        self.config["acquired_seq_count"] = str(int(self.config["acquired_seq_count"]) + 1)

    def extract_config(self, config):
        for category, entry in config:
            self.config[category] = entry
    
    def __str__(self):
        config_tsv = ""
        for key in self.config.keys():
            config_tsv += key + "\t" + str(self.config[key]) + "\n"
        return config_tsv
            
    
    def save_config(self):
        with open(self.config["self_path"], "w") as f:
            f.write(str(self))
            
            
    def set_config(self, category, value):
        if category in self.config.keys():
            self.config[category] = value
        else:
            print(category, "is not a valid key in the config.")
            sys.exit()
    
    def get_config(self, category):
        if category in self.config.keys():
            return self.config[category]
        else:
            print(category, "is not a valid key in the config.")
            sys.exit()

        