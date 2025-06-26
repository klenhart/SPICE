#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ComparisonAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ComparisonAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict, Any, List, Set
import hashlib
import json
import math
import os
import numpy as np
from scipy.spatial.distance import jensenshannon
from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.EWFDHandling.EWFDAssembler import EWFDAssembler


class ComparisonGene:

    def __init__(self, gene_id: str,
                 data_dict_1: Dict[str, List[Any]],
                 data_dict_2: Dict[str, List[Any]],
                 biotype_filter: List[str], tag_filter: List[str],
                 fas_adjacency_matrix: Dict[str, Dict[str, float]]):
        self.gene_id = gene_id
        # Euclidean distance as absolute change measure
        self.edist: float = 0.0
        # rmsd and rmsd_max used for calculation of maximum change possible for that gene (rmsd_rel)
        self.rmsd: float = 0.0
        self.rmsd_max: float = 0.0
        self.rmsd_rel: float = 0.0
        # MPD as measure for isoform diversity for included transcripts given expression data (mpd_expr)
        # It is the diversity actually used
        self.mpd_included: float = 0.0
        self.sd_included: float = 0.0
        # MPD as measure for isoform diversity for all transcripts of that gene that are in the
        # library (mpd_lib)
        self.mpd_all: float = 0.0
        self.sd_all: float = 0.0
        # JSD as emeasure for strength of change
        self.expressed_div_ratio: float = 0.0
        self.jsd: float = 0.0
        self.score: float = 0.0
        self.num_considered_transcripts: int = 0
        if not ComparisonGene.check_gene_not_expressed(data_dict_1, data_dict_2):

            # Remove transcripts with zero expression in both conditions
            data_cond1, data_cond2 = ComparisonGene.indices_expressed_in_one(data_dict_1,
                                                                            data_dict_2)
            # Remove dissimiarity transcripts which show zero expressioin in both conditions
            # from the FAS adjacency matrix -> excludes these entries for rmsd_max calculation
            fas_data_updated = ComparisonGene.update_fas_adjacency_matrix(data_cond1, fas_adjacency_matrix)

            # Extract ewfd vectors calculated using the average relative expression
            # calculated across replicates
            ewfd_1 : List[float] = data_cond1["ewfd_avg_rel_expr"]
            ewfd_2 : List[float] = data_cond2["ewfd_avg_rel_expr"]

            avg_rel_expr_v1 : List[float] = data_cond1["expression_rel_avg"]
            avg_rel_expr_v2 : List[float] = data_cond2["expression_rel_avg"]

            inverted_fas_dict = ComparisonGene.invert_fas_adjacency_matrix(fas_adjacency_matrix)
            updated_inverted_fas = ComparisonGene.invert_fas_adjacency_matrix(fas_data_updated)
            self.num_considered_transcripts: int = len(ewfd_1)

            self.calc_rmsd(ewfd_1, ewfd_2)
            self.calc_max_rmsd(updated_inverted_fas, avg_rel_expr_v1, avg_rel_expr_v2)
            self.rmsd_rel = self.rmsd / self.rmsd_max if self.rmsd_max > 0 else 0.0
            self.calc_euclidean(self.num_considered_transcripts)
            self.calc_jsd(avg_rel_expr_v1, avg_rel_expr_v2)
            self.mpd_included, self.sd_included = ComparisonGene.estimate_diversity(updated_inverted_fas)
            self.mpd_all, self.sd_all = ComparisonGene.estimate_diversity(inverted_fas_dict)
            self.expressed_div_ratio = self.mpd_included/self.mpd_all if self.mpd_all > 0 else 0.0
            self.score = self.rmsd_rel * self.mpd_included

    @staticmethod
    def update_fas_adjacency_matrix(data_dict_1, fas_adjacency_matrix):
        ids_to_keep = set(data_dict_1["ids"])
        for outer_transcript in list(fas_adjacency_matrix.keys()):
            if outer_transcript not in ids_to_keep:
                for inner_transcript in list(fas_adjacency_matrix.keys()):
                    fas_adjacency_matrix[inner_transcript].pop(outer_transcript, None)
                fas_adjacency_matrix.pop(outer_transcript)
        return fas_adjacency_matrix

    @staticmethod
    def indices_expressed_in_one(data_dict_1, data_dict_2):
        delete_set: Set[int] = set()
        for i, entry in enumerate(data_dict_1["expression_rel_avg"]):
            # If there is any transcript that is not expressed in both conditions it is removed from the
            # vector
            if data_dict_1["expression_rel_avg"][i] == 0.0 and data_dict_2["expression_rel_avg"][i] == 0.0:
                delete_set.add(i)
        delete_list = list(delete_set)
        delete_list.sort(reverse=True)
        for index in delete_list:
            for key in data_dict_1.keys():
                data_dict_1[key].pop(index)
                data_dict_2[key].pop(index)
        return data_dict_1, data_dict_2

    @staticmethod
    def check_gene_not_expressed(data_dict_1, data_dict_2):
        if all([x == 0.0 for x in data_dict_1["expression_rel_avg"]]):
            return True
        elif all([x == 0.0 for x in data_dict_2["expression_rel_avg"]]):
            return True
        else:
            return False

    @staticmethod
    def invert_fas_adjacency_matrix(fas_adjacency_matrix: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
        """
        Calculates the invert FAS adjacency matrix in nested dict format to
        turn it into a dissimilarity_matrix.

        Behavior:
            row: Source protein/transcript which is the current row of the matrix.
            col: The target protein which is the current column of the matrix.
            inner: The dictionary mapping from row to all its coll values.
            value: The FAS score from pair row; col.
        """
        return {
            row: {col: 1.0 - value for col, value in inner.items()}
            for row, inner in fas_adjacency_matrix.items()
        }

    @staticmethod
    def dict_to_matrix(inverted_fas_dict: Dict[str, Dict[str, float]]) -> np.ndarray:
        """Convert nested dict of inverted FAS adjacency matrix to numpy matrix."""
        keys_outer = list(inverted_fas_dict.keys()) # row names 
        keys_inner = list(next(iter(inverted_fas_dict.values())).keys()) # column names
        matrix = np.array([[inverted_fas_dict[row][col] for col in keys_inner] for row in keys_outer])
        return matrix
    

    def calc_rmsd(self, ewfd_1: List[float], ewfd_2: List[float]):
        squared_delta_list = [(x - y)**2 for x, y in zip(ewfd_1, ewfd_2)]
        if len(ewfd_1) == 0:
            self.rmsd = 0.0
        else:
            self.rmsd = math.sqrt(sum(squared_delta_list) / len(ewfd_1))

    def calc_euclidean(self, num_transcripts):
        self.edist = self.rmsd * math.sqrt(num_transcripts)

    
    def calc_max_rmsd(self,
                      inverted_fas: Dict[str, Dict[str, float]],
                      avg_rel_expr_v1: List[float], 
                      avg_rel_expr_v2: List[float]):
        """
        Note:
            RMSD_max is the maximum possible RMSD that could result from v,
            if v were perfectly aligned with the most "sensitive" or
            "amplifying" direction defined by the dissimilarity structure in M.
        """

        # Transform average relative expression vectors into numpy arrays
        a = np.array(avg_rel_expr_v1)
        b = np.array(avg_rel_expr_v2)
        # Calculate v as the change vector (a-b) in relative
        # expression between both conditions
        v = a - b

        diss_matrix = ComparisonGene.dict_to_matrix(inverted_fas)
        # M is dissimilarity-induced metric, representing how pairwise
        # expression differences are weighted according to FA dissimilarity
        M = diss_matrix.T @ diss_matrix 
        # Compute the maximum eigenvalue of M
        lambda_max = np.max(np.linalg.eigvalsh(M))

        n = len(v)
        # Calculate squared norm of vector v (||v||**2) -> strength of change
        norm_v_squared = np.linalg.norm(v)**2
        self.rmsd_max = np.sqrt((1 / n) * lambda_max * norm_v_squared) if n > 0 else 0.0

    @staticmethod
    def estimate_diversity(inverted_fas: Dict[str, Dict[str, float]]) -> tuple[float, float]:
        keys = list(inverted_fas.keys())
        num_isoforms = len(keys)
        if num_isoforms < 2:
            return 0.0, 0.0

        values = []
        for i in keys:
            for j in keys:
                if i != j:
                    values.append(inverted_fas[i][j])

        n = len(values)
        if n == 0:
            return 0.0, 0.0

        mean = sum(values) / n
        variance = sum((x - mean) ** 2 for x in values) / n
        sd = math.sqrt(variance)

        return mean, sd

    def calc_jsd(self, avg_rel_expr_v1: List[float], avg_rel_expr_v2: List[float]):
        """
        Calculates Jensen-Shannon distance between two relative expression vectors.
        Stores the result in self.jsd.

        Both input vectors should be non-negative and sum to 1 (i.e., relative expression).
        If they don't, the function will normalize them automatically.
        """
        a = np.array(avg_rel_expr_v1, dtype=np.float64)
        b = np.array(avg_rel_expr_v2, dtype=np.float64)

        # Normalize to ensure valid probability distributions -> due to floating point errors
        # might not be exactly 1
        if a.sum() != 0:
            a /= a.sum()
        if b.sum() != 0:
            b /= b.sum()

        # Compute JSD using scipy (returns distance, not divergence)
        self.jsd = jensenshannon(a, b, base=2)

    def __str__(self):
        output: str = ",".join([self.gene_id, str(self.rmsd)])
        return output

    # @staticmethod
    # def recalculate_ewfd(data_dict: Dict[str, List[Any]], fas_adjacency_matrix: Dict[str, Dict[str, float]]):
    #     return EWFDAssembler.calculate_ewfd(fas_adjacency_matrix,
    #                                         data_dict["expression_rel_avg"],
    #                                         data_dict["ids"])

    # @staticmethod
    # def apply_biotype_filter(data_dict: Dict[str, List[Any]], biotype_filter: List[str]):
    #     """
    #     Removes all transcript entries from the data_dict that match any biotype in the biotype_filter list.

    #     This function operates by identifying indices of transcripts whose biotype is listed in the 
    #     biotype_filter. For each such index, it removes corresponding entries across all keys in 
    #     the data_dict (e.g., 'ids', 'biotypes', 'expression_rel_avg', etc.), preserving structural alignment.

    #     Parameters:
    #         data_dict (Dict[str, List[Any]]): A dictionary containing transcript-associated lists 
    #                                         for a gene (e.g., 'ids', 'biotypes', 'expression_rel_avg').
    #         biotype_filter (List[str]): A list of biotypes to filter out (e.g., ['retained_intron']).

    #     Returns:
    #         Dict[str, List[Any]]: The modified data_dict with filtered transcripts removed.
    #     """
    #     delete_set = set()
    #     for entry in biotype_filter:
    #         for i, biotype in enumerate(data_dict["biotypes"]):
    #             if biotype == entry:
    #                 delete_set.add(i)
    #     delete_list = list(delete_set)
    #     delete_list.sort(reverse=True)
    #     for index in delete_list:
    #         for key in data_dict.keys():
    #             data_dict[key].pop(index)
    #     return data_dict

    # @staticmethod
    # def apply_tag_filter(data_dict: Dict[str, List[Any]], tag_filter: List[str]):
    #     """
    #     Removes all transcript entries from the data_dict that contain any tag in the tag_filter list.

    #     Similar to apply_biotype_filter, this function scans the 'tags' field for each transcript.
    #     If any tag matches the tag_filter, it deletes that transcript's information across all fields
    #     (e.g., 'ids', 'expression_rel_avg', etc.) to ensure consistent indexing.

    #     Parameters:
    #         data_dict (Dict[str, List[Any]]): A dictionary containing transcript-associated lists 
    #                                         for a gene (e.g., 'ids', 'tags', 'expression_rel_avg').
    #         tag_filter (List[str]): A list of tag values to filter out (e.g., ['incomplete', 'cds_start_NF']).

    #     Returns:
    #         Dict[str, List[Any]]: The modified data_dict with filtered transcripts removed.
    #     """
    #     delete_set = set()
    #     for entry in tag_filter:
    #         for i, tags in enumerate(data_dict["tags"]):
    #             if entry in tags:
    #                 delete_set.add(i)
    #     delete_list = list(delete_set)
    #     delete_list.sort(reverse=True)
    #     for index in delete_list:
    #         for key in data_dict.keys():
    #             data_dict[key].pop(index)
        # return data_dict

    def __lt__(self, other):
        return self.rmsd < other.rmsd


class ComparisonAssembler:

    def __init__(self,
                 condition_1: str,
                 condition_2: str,
                 result_pass_path: PassPath,
                 result_info: Dict[str, Any]):
        self.condition_1: str = condition_1
        self.condition_2: str = condition_2
        self.result_pass_path: PassPath = result_pass_path
        self.info: Dict[str, Any] = result_info

        self.comparison_gene_list: List[ComparisonGene] = list()

        with open(os.path.join(result_pass_path["library_path"], "fas_data", "fas_index.json"), "r") as f:
            self.fas_index: Dict[str, str] = json.load(f)

        self.fas_scores_directory: str = os.path.join(result_pass_path["library_path"], "fas_data", "fas_scores")
        # Path to ewfd files for conditions
        self.condition_1_path: str = self.info["expression_imports"]["conditions"][condition_1]["ewfd_path"]
        self.condition_2_path: str = self.info["expression_imports"]["conditions"][condition_2]["ewfd_path"]

        self.biotype_filter: List[str] = list()
        self.tag_filter: List[str] = list()
        self.filter_hash = md5_hash("".join(self.biotype_filter+self.tag_filter), 6)

    def add_biotype_filter(self, filter_out: str):
        """ 
        Adds specified transcript biotype to attribute biotype_filter
        and recalculates filter setting ID using md5_hash.

        Parameters:
            filter_out: str: name of the transcript biotype to exclude from ewfd
                             calculation.
        
        Note:
            - Currently biotype `nonsense_mediated_decay` and `non_coding` is set as default.
            - `non_coding` not included in spice library, only `protein_coding` and `nonsense_mediated_decay`
            - Change this to allow user to choose.
        """
        self.biotype_filter.append(filter_out)
        self.filter_hash = md5_hash("".join(self.biotype_filter + self.tag_filter), 6)

    def add_tag_filter(self, filter_out: str):
        """ 
        Adds specified transcript tag to attribute tag_filter
        and recalculates filter setting ID using md5_hash.

        Parameters:
            filter_out: str: name of the transcript tag to exclude from ewfd
                             calculation.
        
        Note:
            - Currently tag `incomplete` is set as default. However, if no
              protein sequence is available and it is not NMD, this transcript is 
              not included in the Spice library.
            - Change this to allow user to choose.
        """
        self.tag_filter.append(filter_out)
        self.filter_hash = md5_hash("".join(self.biotype_filter + self.tag_filter), 6)

    def compare_genes(self):
        # Reads ewfd/conditions data for condition 1
        with open(self.condition_1_path, "r") as f_1:
            data_cond_1: Dict[str, Dict[str, List[Any]]] = json.load(f_1)["data"]
        # Reads ewfd/conditions data for condition 2
        with open(self.condition_2_path, "r") as f_2:
            data_cond_2: Dict[str, Dict[str, List[Any]]] = json.load(f_2)["data"]
        # Loops over each gene id and loads FAS score data. Creates ComparisonGene object
        for gene_id in data_cond_1.keys():
            with open(os.path.join(self.fas_scores_directory, self.fas_index[gene_id]), "r") as f:
                fas_adjacency_matrix: Dict[str, Dict[str, float]] = json.load(f)[gene_id]
            self.comparison_gene_list.append(ComparisonGene(gene_id,
                                                            data_cond_1[gene_id],
                                                            data_cond_2[gene_id],
                                                            self.biotype_filter,
                                                            self.tag_filter,
                                                            fas_adjacency_matrix))

    def sort_genes_by_score(self):
        self.comparison_gene_list.sort(key=lambda g: g.score, reverse=True)


    def delete_zero_rmsd_genes(self):
        delete_list: List[int] = list()
        for i, comparison_gene in enumerate(self.comparison_gene_list):
            if comparison_gene.rmsd == 0:
                delete_list.append(i)
        for i in sorted(delete_list, reverse=True):
            self.comparison_gene_list.pop(i)

    def __str__(self):
        headers = [
            "gene_id", "edist", "rmsd", "rmsd_rel",
            "mpd_included", "sd_included", "mpd_all", "sd_all",
            "expressed_div_ratio", "jsd", "score", "num_considered_transcripts"
        ]
        rows = ["\t".join(headers)]
        for g in self.comparison_gene_list:
            row = [
                g.gene_id, f"{g.edist:.5f}", f"{g.rmsd:.5f}",
                f"{g.rmsd_rel:.5f}", f"{g.mpd_included:.5f}", f"{g.sd_included:.5f}",
                f"{g.mpd_all:.5f}", f"{g.sd_all:.5f}", f"{g.expressed_div_ratio:.5f}",
                f"{g.jsd:.5f}", f"{g.score:.5f}", str(g.num_considered_transcripts)
            ]
            rows.append("\t".join(row))
        return "\n".join(rows)

    def save(self, out_dir):
        output_path = os.path.join(
            out_dir,
            self.condition_1 + "@" + self.condition_2 + "@" + self.filter_hash + ".tsv"
        )
        with open(output_path, "w") as f:
            f.write(str(self))


def md5_hash(data, length: int = 32):
    """
    Encodes filter settings (biotypes and tags) into a short ID used in
    mode "compare" result file names.
    """
    hasher = hashlib.md5()
    hasher.update(data.encode("utf-8"))
    hash_string = hasher.hexdigest()
    return hash_string[:length]
