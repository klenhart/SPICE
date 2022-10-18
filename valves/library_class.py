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
        self.config["distance_master_path"] = None
        self.config["distance_master_lcr_path"] = None
        self.config["distance_master_tmhmm_path"] = None
        self.config["lcr_path"] = None
        self.config["tmhmm_path"] = None
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
        self.config["fas_buffer_tmhmm_path"] = None
        self.config["fas_buffer_lcr_path"] = None
        self.config["annotation_path"] = None
        self.config["canonical_path"] = None
        
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
        self.config["flag_assemble_canonical_transcripts"] = "False"

        if not flag_init:
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

        