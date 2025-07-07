#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ConditionAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ConditionAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import os

from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class ConditionAssembler:

    def __init__(self,
                 library_pass_path: PassPath,
                 condition_name: str = "",
                 initial_flag: bool = False):
        if initial_flag:
            # Initialize ConditionAssembler datastructure condition_assembly
            self.library_pass_path: PassPath = library_pass_path
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = library_pass_path["root"]
            self.condition_assembly["normalization"]: str = ""
            # Will store names of replicates
            self.condition_assembly["replicates"]: List[str] = list()
            # Will store number of replicates
            self.condition_assembly["replicate_count"]: int = 0
            # Condition data
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
            # gene_assembly loaded from library path
            gene_assembly: Dict[str, Gene] = self.__load_gene_assembly()
            for gene_id in gene_assembly.keys():
                # Initiize condition_assembly["data"] structure -> gene-top-level structure
                self.condition_assembly["data"][gene_id]: Dict[str, List[Any]] = dict()
                # Will store Transcript and Protein ids for that gene
                self.condition_assembly["data"][gene_id]["ids"]: List[str] = list()
                # Will store synonyms of transcripts for that gene
                self.condition_assembly["data"][gene_id]["synonyms"]: List[List[str]] = list()
                # Will store biotypes
                self.condition_assembly["data"][gene_id]["biotypes"]: List[str] = list()
                # Will store tsl
                self.condition_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = list()
                self.condition_assembly["data"][gene_id]["tags"]: List[List[str]] = list()
                # Stores average relative expression of each transcript/protein across all replicates in condition
                self.condition_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = list()
                # Stores relative expression of each transcript in each replicate in
                # list of lists: Each transcript will represented here by a list of length = number of replicates
                # storing the relative expression in each of them
                self.condition_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = list()
                # Stores total (normalized) expression of each transcript/protein across all replicates in condition
                # Same format as expression_rel_all but with absolute values
                self.condition_assembly["data"][gene_id]["expression_all"]: List[List[float]] = list()
                # Loop over all tanscripts in gene_assembly of library and extend information for each gene
                for transcript in gene_assembly[gene_id].get_transcripts():
                    self.condition_assembly["data"][gene_id]["ids"].append(transcript.get_id())
                    self.condition_assembly["data"][gene_id]["synonyms"].append(transcript.get_synonyms())
                    biotype: str = transcript.get_biotype()
                    self.condition_assembly["data"][gene_id]["biotypes"].append(biotype)
                    tsl: int = transcript.get_transcript_support_level()
                    self.condition_assembly["data"][gene_id]["transcript_support_levels"].append(tsl)
                    self.condition_assembly["data"][gene_id]["tags"].append(transcript.get_tags())
                    # Initially, set average relative expression of each transcript/protein to zero
                    self.condition_assembly["data"][gene_id]["expression_rel_avg"].append(1.0)
                    # Append list for each transcript
                    self.condition_assembly["data"][gene_id]["expression_rel_all"].append([])
                    self.condition_assembly["data"][gene_id]["expression_all"].append([])
        else:
            # initialize attributes that will later be loaded using internal func load.
            self.library_pass_path: PassPath = library_pass_path
            self.condition_assembly: Dict[str, Any] = dict()
            self.condition_assembly["name"]: str = condition_name
            self.condition_assembly["library"]: str = library_pass_path["root"]
            self.condition_assembly["normalization"]: str = ""
            self.condition_assembly["replicates"]: List[str] = list()
            self.condition_assembly["replicate_count"]: int = 0
            self.condition_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

    def save(self, output_path: str):
        with open(output_path, "w") as f:
            json.dump(self.condition_assembly, f, indent=4)

    def load(self, input_path: str):
        """
        Loads condition data from disk into an existing ConditionAssembler instance.

        Parameters:
            input_path: Usually the condition paths (/expression/conditions/expression_[cond_name].json
                        containing information (e.g. relative expression of transcripts
                        of each gene in each replicate and average relative expression across)
                        about a previously populated ConditionAssembler instance which data was saved
                        on disk.

        Behavior:
            - Loads information and sets condition_assembly attribute
            - Sets library_pass_path attribute to direct to paths.json in library dir

        Use Case:
            Used for condition data reuse without recalcuation during result
            mode "condition". Will typically be used when initial_flag = False after creating
            a ConditionAssembler object in EWFDAssembler, created with condition_flag = True.
        """
        with open(input_path, "r") as f:
            self.condition_assembly = json.load(f)

    def insert_expression(self, expression_assembler: ExpressionAssembler):
        """
        Populates the condition_assembly "data" structure through saving original and relative
        expression of each transcript for that replicate and calculating the average relative
        expression over all replicates.

        New behavior:
            - For each gene:
                1. Collect per-replicate expression for each transcript.
                2. Compute the mean expression per transcript.
                3. Normalize to get final relative expression.
        Note:
            - Applied on each recreated replicate ExpressionAssembler object of a condition.
              This will be called in ResultBuddy's build_condition()
        """
        # Load expression_assembly dict of a replicate belonging to this condition
        expr_assembly: Dict[str, Any] = expression_assembler.expression_assembly
        # Add replicate name
        self.condition_assembly["replicates"].append(expr_assembly["name"])
        # Increase replicate_count_number
        self.condition_assembly["replicate_count"] += 1
        # Loop over all genes
        for gene_id in self.condition_assembly["data"].keys():
            # Checks if gene can be found in the expression_assembly of that replicate
            if gene_id in expr_assembly["data"].keys():
                gene_found_flag: bool = True
            else:
                gene_found_flag: bool = False
            # Loop over transcript/protein ids of that gene in condition assembly
            for i, transcript_id in enumerate(self.condition_assembly["data"][gene_id]["ids"]):
                # This is the case if the expression gtf file only contained genes with non-zero expression values
                # Deprecated: ExpressionAssembler's cleanse_assembly() is not used in ResultBuddy
                # (commented out) which means transcripts with zero expression are maintained in the expression
                # assembly
                if not gene_found_flag:
                    expr_value = 0.0
                    rel_expr_value = 0.0
                # Gene is present in expression_assembly of that replicate as well as transcript
                elif transcript_id in expr_assembly["data"][gene_id]["ids"]:
                    index: int = expr_assembly["data"][gene_id]["ids"].index(transcript_id)
                    rel_expr_value = expr_assembly["data"][gene_id]["expression_rel"][index]
                    expr_value = expr_assembly["data"][gene_id]["expression"][index]
                # Gene is present in expression_assembly but transcript is not 
                # This is a deprecated option as ExpressionAssembler's cleanse_assembly() is not used in ResultBuddy
                # (commented out) which means transcripts with zero expression are maintained in the expression
                # assembly
                else:
                    expr_value = 0.0
                    rel_expr_value = 0.0
                # Add relative expression and expression values of replicate to datastructure
                self.condition_assembly["data"][gene_id]["expression_rel_all"][i].append(rel_expr_value)
                self.condition_assembly["data"][gene_id]["expression_all"][i].append(expr_value)

            # === NEW: Calculate averaged relative expression vector from averaged absolute expressions ===
            expr_all = self.condition_assembly["data"][gene_id]["expression_all"]
            replicate_count = self.condition_assembly["replicate_count"]

            # Average expression values across replicates per transcript
            expr_means = [sum(vals) / replicate_count for vals in expr_all]

            # Normalize to get average relative expression
            expr_total = sum(expr_means)
            if expr_total > 0.0:
                expr_rel_avg = [v / expr_total for v in expr_means]
            else:
                expr_rel_avg = [0.0 for _ in expr_means]

            self.condition_assembly["data"][gene_id]["expression_rel_avg"] = expr_rel_avg

    def __load_gene_assembly(self) -> Dict[str, Gene]:
        """
        Loads gene, transcript, protein, and FAS data from Disk and creates GeneAssembler.gene_assembly
        like dictionary

        Behavior:
            - Loads transcript_info.json, sequences.json, fas_index.json and FAS score JSON files.
            - Employs GeneAssembler's from_dict() function to 
        
        Returns:
            gene_assembly: Dict[str, Gene]: Dictionary mapping gene IDs to Gene instances, fully restored with
                           their transcripts, proteins, sequences, and FAS scores.
        """
        with open(self.library_pass_path["transcript_info"], "r") as f:
            info_dict: Dict[str, Dict[str, Any]] = json.load(f)
        with open(self.library_pass_path["transcript_seq"], "r") as f:
            seq_dict: Dict[str, Dict[str, Any]] = json.load(f)
        fas_dict: Dict[str, Dict[str, Any]] = dict()
        with open(self.library_pass_path["fas_index"], "r") as f1:
            fas_index: Dict[str, str] = json.load(f1)
            for path in set(fas_index.values()):
                with open(os.path.join(self.library_pass_path["fas_scores"], path), "r") as f2:
                    fas_sub_dict: Dict[str, Dict[str, Any]] = json.load(f2)
                    fas_dict.update(fas_sub_dict)
        gene_assembly: Dict[str, Gene] = GeneAssembler.from_dict(info_dict, seq_dict, fas_dict)
        return gene_assembly

    def cleanse_assembly(self):
        """ DEPRECATED """
        cleanse_dict: Dict[str, List[str]] = dict()
        for gene_id in self.condition_assembly["data"].keys():
            cleanse_dict[gene_id] = list()
            id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
            for index, transcript_id in enumerate(id_list):
                if self.condition_assembly["data"][gene_id]["expression_rel_avg"][index] == 0.0:
                    cleanse_dict[gene_id].append(transcript_id)
        for gene_id in cleanse_dict.keys():
            for transcript_id in cleanse_dict[gene_id]:
                id_list: List[str] = self.condition_assembly["data"][gene_id]["ids"]
                index: int = id_list.index(transcript_id)
                self.condition_assembly["data"][gene_id]["ids"].pop(index)
                self.condition_assembly["data"][gene_id]["synonyms"].pop(index)
                self.condition_assembly["data"][gene_id]["biotypes"].pop(index)
                self.condition_assembly["data"][gene_id]["transcript_support_levels"].pop(index)
                self.condition_assembly["data"][gene_id]["tags"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_avg"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_rel_all"].pop(index)
                self.condition_assembly["data"][gene_id]["expression_all"].pop(index)

            if len(self.condition_assembly["data"][gene_id]["ids"]) == 0:
                del self.condition_assembly["data"][gene_id]
