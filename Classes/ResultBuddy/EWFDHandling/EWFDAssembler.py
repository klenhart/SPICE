#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  EWFDAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EWFDAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import json
import math
import numpy as np

from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


class EWFDAssembler:

    def __init__(self,
                 species: str = "",
                 taxon_id: int = "",
                 library_pass_path: PassPath = PassPath(dict()),
                 expression_path: str = "",
                 initial_flag: bool = False,
                 condition_flag: bool = False):
        if initial_flag:
            # If mode condition
            if condition_flag:
                # Initialize ConditionAssembler object
                condition_assembler: ConditionAssembler = ConditionAssembler(library_pass_path)
                # Populate ConditionAssembler datastructure from disk
                condition_assembler.load(expression_path)

                self.library_pass_path: PassPath = library_pass_path
                self.ewfd_assembly: Dict[str, Any] = dict()
                self.ewfd_assembly["name"]: str = condition_assembler.condition_assembly["name"]
                self.ewfd_assembly["library"]: str = condition_assembler.condition_assembly["library"]
                self.ewfd_assembly["normalization"] = condition_assembler.condition_assembly["normalization"]
                self.ewfd_assembly["replicates"] = condition_assembler.condition_assembly["replicates"]
                self.ewfd_assembly["replicate_count"] = condition_assembler.condition_assembly["replicate_count"]
                self.ewfd_assembly["data"]: Dict[str, Dict[str, Any]] = dict()

                gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
                # Populate GeneAssembler datastructure from disk
                gene_assembler.load(self.library_pass_path)
                # Extract FAS score matrices
                fas_dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = gene_assembler.get_fas_dist_matrix()
                # Extract data in condition_assembler.condition_assembly
                condition_data: Dict[str, Any] = condition_assembler.condition_assembly["data"]

                for gene_id in condition_data.keys():
                    # Extract FAS score matrix for gene
                    gene_dist_matrix: Dict[str, Dict[str, float]] = fas_dist_matrix[gene_id]
                    # Initiate emtpy dict in ewfd_assembly["data"] for that gene
                    self.ewfd_assembly["data"][gene_id]: Dict[str, Any] = dict()
                    # Extract transcript ids and save in ewfd_assembly
                    transcript_ids: List[str] = condition_data[gene_id]["ids"]
                    self.ewfd_assembly["data"][gene_id]["ids"]: List[str] = transcript_ids
                    # Calculate number of transcripts for gene
                    transcript_count: int = len(self.ewfd_assembly["data"][gene_id]["ids"])
                    # Extract biotypes, tsl and tags of transcripts and save in ewfd_assembly
                    self.ewfd_assembly["data"][gene_id]["biotypes"]: List[str] = condition_data[gene_id]["biotypes"]
                    tsl_list: List[int] = condition_data[gene_id]["transcript_support_levels"]
                    self.ewfd_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = tsl_list
                    self.ewfd_assembly["data"][gene_id]["tags"]: List[List[str]] = condition_data[gene_id]["tags"]
                    # Extract the average relative expression values of that gene over replicates and save in ewfd_assembly
                    expr_rel_avg_list: List[float] = condition_data[gene_id]["expression_rel_avg"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel_avg"]: List[float] = expr_rel_avg_list
                    # Extract the relative expression values of that gene in all replicates and save in ewfd_assembly
                    # Note: condition_data[gene_id]["expression_rel_all"] is a nested list of length = number of transcripts
                    # whereas each transcript is represented by a list of length replicate_count, storing the absolute 
                    # expression values in each replicate for that transcript, shape:  n_transcripts × n_replicates
                    expr_rel_all_list: List[List[float]] = condition_data[gene_id]["expression_rel_all"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel_all"]: List[List[float]] = expr_rel_all_list

                    # These need to be calculated.
                    # Initialize list to save result of ewfd calculation using average relative expression vectors (condition-level)
                    self.ewfd_assembly["data"][gene_id]["ewfd_avg_rel_expr"]: List[float]
                    # Initialize list to save result of ewfd calculation in each replicate for each transcript (replicate_level)
                    # Shape: n_transcripts × n_replicates
                    self.ewfd_assembly["data"][gene_id]["ewfd_all"]: List[List[float]] = list()
                    # Initialize list to save result of maximum ewfd value for each transcript across replicates, initially zero
                    self.ewfd_assembly["data"][gene_id]["ewfd_max"]: List[float] = [0.0] * transcript_count
                    # Initialize list to save result of minimum ewfd value for each transcript across replicates, initially one
                    self.ewfd_assembly["data"][gene_id]["ewfd_min"]: List[float] = [1.0] * transcript_count
                    # Initialize list to save average ewfd for transcripts across replicates (length = num transcrpts)
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd"]: List[float] = list()
                    # Initialize list to save std of ewfd for transcripts across replicates (length = num transcrpts)
                    self.ewfd_assembly["data"][gene_id]["ewfd_std"]: List[float] = list()
                    # Initialize list to save error margins/confidence of ewfd for transcripts across replicates (length = num transcrpts)
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd+std"]: List[float] = list()
                    self.ewfd_assembly["data"][gene_id]["avg_ewfd-std"]: List[float] = list()

                    # Calculate "average" ewfd for gene using average relative expression values across replicates
                    # generated in ConditionAssembler
                    ewfd_avg_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                                  expr_rel_avg_list,
                                                                                  transcript_ids)
                    # Save in ewfd_assembly
                    self.ewfd_assembly["data"][gene_id]["ewfd_avg_rel_expr"] = ewfd_avg_rel_expr
                    # Note: Except of recalcalculating ewfd for each replicate again, the internal function load
                    # could have been used to extract ewfds from /expression/replicates/expression_[replicate_name].json
                    ewfd_all: List[List[float]] = list()
                    # turn shape  n_transcripts × n_replicates to shape n_replicates x n_transcripts
                    for i, repl_rel_expr_list in enumerate(np.array(expr_rel_all_list).transpose()):
                        ewfd_repl_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                                       list(repl_rel_expr_list),
                                                                                       transcript_ids)
                        ewfd_all.append(ewfd_repl_rel_expr)

                        # Set max and min ewfd value for each transcript across replicates
                        for j, ewfd in enumerate(ewfd_repl_rel_expr):
                            if ewfd > self.ewfd_assembly["data"][gene_id]["ewfd_max"][j]:
                                self.ewfd_assembly["data"][gene_id]["ewfd_max"][j] = ewfd
                            if ewfd < self.ewfd_assembly["data"][gene_id]["ewfd_min"][j]:
                                self.ewfd_assembly["data"][gene_id]["ewfd_min"][j] = ewfd
                    # turn shape  n_replicates × n_transcripts to shape n_transcripts x n_replicates
                    ewfd_all = list(np.array(ewfd_all).transpose())
                    # Add ewfd calculated across all replicates to ewfd_assembly
                    self.ewfd_assembly["data"][gene_id]["ewfd_all"] = [list(entry) for entry in ewfd_all]

                    # Calculate the average EWFD for each transcript.
                    # For each transcript ewfd list across replicate
                    for ewfd_list in self.ewfd_assembly["data"][gene_id]["ewfd_all"]:
                        # count: num replicates
                        count: int = len(ewfd_list)
                        # Calculate average ewfd across replicates
                        ewfd_average: float = sum(ewfd_list) / count
                        # Add to ewfd_assembly
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd"].append(ewfd_average)
                        std: float
                        std = math.sqrt(sum([(ewfd_value - ewfd_average)**2 for ewfd_value in ewfd_list]) / count)
                        self.ewfd_assembly["data"][gene_id]["ewfd_std"].append(std)
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd+std"].append(ewfd_average+std)
                        self.ewfd_assembly["data"][gene_id]["avg_ewfd-std"].append(ewfd_average-std)
            # mode expression
            else:
                # Initiate ExpressionAssembler instance
                expression_assembler: ExpressionAssembler = ExpressionAssembler(library_pass_path)
                # Load previously (spice_result mode expression) generated expression data 
                expression_assembler.load(expression_path)
                self.library_pass_path: PassPath = library_pass_path
                # Initiate EWFDAssembler datastructure
                self.ewfd_assembly: Dict[str, Any] = dict()
                self.ewfd_assembly["name"]: str = expression_assembler.expression_assembly["name"]
                self.ewfd_assembly["library"]: str = expression_assembler.expression_assembly["library"]
                self.ewfd_assembly["normalization"] = expression_assembler.expression_assembly["normalization"]
                self.ewfd_assembly["data"]: Dict[str, Dict[str, Any]] = dict()
                # Initiate GeneAssembler instance
                gene_assembler: GeneAssembler = GeneAssembler(species, str(taxon_id))
                # Load previously GeneAssembler data created in spice library generation
                gene_assembler.load(self.library_pass_path)
                # Obtain FAS distance matrix calculated during library generation
                fas_dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = gene_assembler.get_fas_dist_matrix()
                # Obtain expression data (contains expression values and relative expression values)
                expression_data: Dict[str, Any] = expression_assembler.expression_assembly["data"]
                # Loop over all genes
                for gene_id in expression_data.keys():
                    # Obtain FAS matrix of gene
                    gene_dist_matrix: Dict[str, Dict[str, float]] = fas_dist_matrix[gene_id]
                    # Start poupulating ewfd_assembly["data"]
                    self.ewfd_assembly["data"][gene_id]: Dict[str, Any] = dict()
                    transcript_ids: List[str] = expression_data[gene_id]["ids"]
                    self.ewfd_assembly["data"][gene_id]["ids"]: List[str] = transcript_ids
                    self.ewfd_assembly["data"][gene_id]["biotypes"]: List[str] = expression_data[gene_id]["biotypes"]
                    tsl_list: List[int] = expression_data[gene_id]["transcript_support_levels"]
                    self.ewfd_assembly["data"][gene_id]["transcript_support_levels"]: List[int] = tsl_list
                    self.ewfd_assembly["data"][gene_id]["tags"]: List[List[str]] = expression_data[gene_id]["tags"]
                    expression_rel: List[float] = expression_data[gene_id]["expression_rel"]
                    self.ewfd_assembly["data"][gene_id]["expression_rel"]: List[float] = expression_rel

                    # Initiate empty ewfd list for each gene. These need to be calculated.
                    self.ewfd_assembly["data"][gene_id]["ewfd_rel_expr"]: List[float]

                    # Calculatino of ewfd for each gene:
                    ewfd_rel_expr: List[float] = EWFDAssembler.calculate_ewfd(gene_dist_matrix,
                                                                              expression_rel,
                                                                              transcript_ids)
                    # Save result list in datastructure
                    self.ewfd_assembly["data"][gene_id]["ewfd_rel_expr"] = ewfd_rel_expr
        else:
            self.ewfd_assembly: Dict[str, Any] = dict()

    def save(self, output_path) -> None:
        """
        Saves attribute ewfd_assembly as JSON file in defined path.

        Note:
            Usually /ewfd/replicates or ewfd/conditions.

        Use Case:
            Part of result files.
        """
        with open(output_path, "w") as f:
            json.dump(self.ewfd_assembly, f, indent=4)

    def load(self, input_path) -> None:
        with open(input_path, "r") as f:
            self.ewfd_assembly = json.load(f)

    @staticmethod
    def calculate_ewfd(gene_fas_dists: Dict[str, Dict[str, float]],
                       rel_expressions: List[float],
                       transcript_ids: List[str]) -> List[float]:
        """
        Calculates expression-weighted-functional-disturbance (EWFD) vectors for each gene.

        Parameters:
            - gene_fas_dists: Dict[str, Dict[str, float]]: Ordered FAS score matrix for a Gene
                              in form of a dictionary. Usually obtained after intialization of
                              EWFDAssembler instance.
            - rel_expressions: List[float]: List containing relative expression values of tran-
                               scripts of that gene.
            - transcript_ids: List[str]): List containing transcript ids.

        Returns:
            ewfd_list: List[float]: List containing ewfd values for each transcript, rounded to
                                    4 digits.
        
        Behavior:
            - Effectively performs a Matrix-Vector multiplication with the Matrix containing pairwise
              FAS scores between selected transcripts/proteins of a gene and the Vector containing the
              relative expression values for each transcript.
            - Uses 1-FAS score which reflects FA dissimilarity (FAS-score represents FA similarity)
        
        Note
            Calculation is either done for each sample (each replicate in any condition) or on a
            condition level. On which level is determined by the condition_flag which is passed 
            to build the EWFDAssembler instance and depends on the mode selected when runnning 
            spice_result.py.
        """
        # Initiate ewfd list/vector with zeros
        ewfd_list: List[float] = [0.0] * len(transcript_ids)
        # Performs effectively a Matrix-Vector multiplication with complement FAS scores (dissimilarity score)
        for s, seed_id in enumerate(transcript_ids):
            for q, query_id in enumerate(transcript_ids):
                ewfd_list[s] += round(rel_expressions[q] * (1 - gene_fas_dists[seed_id][query_id]), 4)

        # ewfd_list = [round(1 - movement, 4) for movement in ewfd_list]

        return ewfd_list
