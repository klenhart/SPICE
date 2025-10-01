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
import copy
import numpy as np
from scipy.spatial.distance import jensenshannon
from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.EWFDHandling.EWFDAssembler import EWFDAssembler


class ComparisonGene:

    import copy
from typing import Any, Dict, List
import numpy as np

class ComparisonGene:
    def __init__(self, gene_id: str,
                 data_dict_1: Dict[str, List[Any]],
                 data_dict_2: Dict[str, List[Any]],
                 biotype_filter: List[str], tag_filter: List[str],
                 fas_adjacency_matrix: Dict[str, Dict[str, float]]):
        self.gene_id = gene_id

        # --- metrics / outputs ---
        self.rmsd: float = 0.0
        self.rmsd_max: float = 0.0
        self.rmsd_rel: float = 0.0
        # self.cosine_diss: float = 0.0

        self.mpd_included: float = 0.0
        self.sd_included: float = 0.0
        # self.mpd_all: float = 0.0
        # self.sd_all: float = 0.0

        # self.expressed_div_ratio: float = 0.0
        self.jsd: float = 0.0
        # self.score: float = 0.0
        self.num_considered_transcripts: int = 0

        # If no expression in both conditions, keep defaults and return.
        if ComparisonGene.check_gene_not_expressed(data_dict_1, data_dict_2):
            return

        # ---------------------------------------------------------------------
        # 1) Filter to transcripts expressed in at least one condition
        #    (this function should also apply any biotype/tag filters upstream)
        #    It must return dicts that already contain ONLY the kept transcripts.
        # ---------------------------------------------------------------------
        data_cond1, data_cond2 = ComparisonGene.indices_expressed_in_one(
            data_dict_1, data_dict_2
        )  # filtered dicts

        # Aligned IDs and expression vectors AFTER filtering
        ids_used: List[str] = list(data_cond1["ids"])
        expr1_used: List[float] = list(data_cond1["expression_rel_avg"])
        expr2_used: List[float] = list(data_cond2["expression_rel_avg"])

        # Safety: both filtered lists must match and be same length
        assert ids_used == list(data_cond2["ids"]), "Filtered IDs must match across conditions."
        assert len(ids_used) == len(expr1_used) == len(expr2_used), "Lengths must align after filtering."

        self.num_considered_transcripts = len(ids_used)

        # ---------------------------------------------------------------------
        # 2) Build FAS graph for USED transcripts only (for rmsd_max & mpd_included)
        #    Keep also the full inverted graph for mpd_all.
        # ---------------------------------------------------------------------
        fas_used = ComparisonGene.update_fas_adjacency_matrix(
            ids_used, copy.deepcopy(fas_adjacency_matrix)
        )
        fas_used_inv = ComparisonGene.invert_fas_adjacency_matrix(fas_used)

        # fas_full_inv = ComparisonGene.invert_fas_adjacency_matrix(fas_adjacency_matrix)

        # Sanity: graph node set matches ids_used
        assert set(fas_used_inv.keys()) == set(ids_used), "FAS(used) must match filtered IDs."

        # ---------------------------------------------------------------------
        # 3) Recompute EWFD vectors for the SAME filtered set (order must match ids_used)
        # ---------------------------------------------------------------------
        # data_cond1 / data_cond2 are already filtered; recalculate_ewfd should respect that order
        ewfd_1: List[float] = ComparisonGene.recalculate_ewfd(data_cond1, fas_adjacency_matrix)
        ewfd_2: List[float] = ComparisonGene.recalculate_ewfd(data_cond2, fas_adjacency_matrix)

        assert len(ewfd_1) == len(ewfd_2) == len(ids_used), "EWFD and IDs must have same length."

        # ---------------------------------------------------------------------
        # 4) Raw RMSD, RMSD_max (on USED set only), RMSD_rel
        # ---------------------------------------------------------------------
        self.rmsd = ComparisonGene.calc_rmsd(ewfd_1, ewfd_2)

        # IMPORTANT: pass ONLY filtered expressions/IDs + filtered FAS
        self.calc_max_rmsd(fas_used_inv, expr1_used, expr2_used, ids_used)
        self.rmsd_rel = self.rmsd / self.rmsd_max if self.rmsd_max > 0 else 0.0

        # Cosine dissimilarity (on EWFD for the filtered set)
        # self.cosine_diss = ComparisonGene.calc_inverted_cosine_similarity(ewfd_1, ewfd_2)

        # ---------------------------------------------------------------------
        # JSD on filtered & renormalized relative expression vectors
        # (safe if sums changed slightly after trimming)
        # ---------------------------------------------------------------------
        p = np.asarray(expr1_used, dtype=float)
        q = np.asarray(expr2_used, dtype=float)
        p = p / p.sum() if p.sum() > 0 else p
        q = q / q.sum() if q.sum() > 0 else q
        self.jsd = ComparisonGene.calc_jsd(p.tolist(), q.tolist())

        # ---------------------------------------------------------------------
        # 6) Diversity (MPD) metrics
        #    - included: only used transcripts
        #    - all:      all transcripts in the library for this gene
        # ---------------------------------------------------------------------
        self.mpd_included, self.sd_included = ComparisonGene.estimate_diversity(fas_used_inv)
        # self.mpd_all, self.sd_all = ComparisonGene.estimate_diversity(fas_full_inv)

    @staticmethod
    def update_fas_adjacency_matrix(ids_to_keep, fas_adjacency_matrix):
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
    def dict_to_matrix(inverted_fas_dict: Dict[str, Dict[str, float]], transcript_order: List[str]) -> np.ndarray:
        return np.array([
            [inverted_fas_dict[row][col] for col in transcript_order]
            for row in transcript_order
        ])

    @staticmethod
    def calc_rmsd(ewfd_1: List[float], ewfd_2: List[float]):
        squared_delta_list = [(x - y)**2 for x, y in zip(ewfd_1, ewfd_2)]
        if len(ewfd_1) == 0:
            rmsd = 0.0
        else:
            rmsd = math.sqrt(sum(squared_delta_list) / len(ewfd_1))
        return rmsd

    def calc_max_rmsd(self,
                inverted_fas: Dict[str, Dict[str, float]],
                avg_rel_expr_v1: List[float], 
                avg_rel_expr_v2: List[float],
                transcript_order: List[str]):
        """
        Calculates the theoretical maximum RMSD based on the change vector v and
        the largest eigenvalue of the dissimilarity-induced metric M.

        Parameters:
            inverted_fas (Dict[str, Dict[str, float]]): Nested dict (dissimilarity matrix)
                                                        where values are 1 - FAS scores.
            avg_rel_expr_v1 (List[float]): Relative expression vector for condition 1.
            avg_rel_expr_v2 (List[float]): Relative expression vector for condition 2.
            transcript_order (List[str]): Ordered list of transcript IDs, ensures alignment
                                        between expression vectors and dissimilarity matrix.
        """
        a = np.array(avg_rel_expr_v1)
        b = np.array(avg_rel_expr_v2)
        v = a - b  # Change vector

        # Ensure matrix is aligned to transcript order
        diss_matrix = ComparisonGene.dict_to_matrix(inverted_fas, transcript_order)

        # Compute M = Dᵗ D
        M = diss_matrix.T @ diss_matrix

        # Eigenvalue decomposition (symmetric matrix)
        lambda_max = np.max(np.linalg.eigvalsh(M))

        n = len(v)
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

    @staticmethod
    def calc_jsd(vec1: List[float], vec2: List[float]) -> float:
        a = ComparisonAssembler.normalize(vec1)
        b = ComparisonAssembler.normalize(vec2)
        if a is None or b is None:
            return 1.0  # Max distance
        return jensenshannon(a, b, base=2.0)

    # def __str__(self):
    #     output: str = ",".join([self.gene_id, str(self.rmsd)])
    #     return output

    @staticmethod
    def recalculate_ewfd(data_dict: Dict[str, List[Any]], fas_adjacency_matrix: Dict[str, Dict[str, float]]):
        return EWFDAssembler.calculate_ewfd(fas_adjacency_matrix,
                                            data_dict["expression_rel_avg"],
                                            data_dict["ids"])

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
    
    @staticmethod
    def normalize(vec):
        vec = np.array(vec, dtype=np.float64)
        return vec / vec.sum() if vec.sum() > 0 else None


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
        self.comparison_gene_list.sort(key=lambda g: g.rmsd, reverse=True)


    def delete_zero_rmsd_genes(self):
        delete_list: List[int] = list()
        for i, comparison_gene in enumerate(self.comparison_gene_list):
            if comparison_gene.rmsd == 0:
                delete_list.append(i)
        for i in sorted(delete_list, reverse=True):
            self.comparison_gene_list.pop(i)

    def __str__(self):
        headers = [
            "gene_id", "rmsd", "rmsd_rel",
            "mpd", "jsd", "AT-bin"
        ]
        rows = ["\t".join(headers)]
        for g in self.comparison_gene_list:
            row = [
                g.gene_id, f"{g.rmsd:.5f}",
                f"{g.rmsd_rel:.5f}", f"{g.mpd_included:.5f}",
                f"{g.jsd:.5f}", str(g.num_considered_transcripts)
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
