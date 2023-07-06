#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  AssemblyVisualizer is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AssemblyVisualizer is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.SequenceHandling.GeneAssembler import GeneAssembler

import json
import math
from typing import Dict, List, Any, Tuple
import argparse
import numpy as np
import os
import plotly.express as px
import matplotlib.pyplot as plt

import pandas


class AssemblyVisualizer:

    incomplete_tags: List[str] = ['mRNA_end_NF', 'cds_end_NF', 'mRNA_start_NF', 'cds_start_NF']

    def __init__(self, gene_assembler: GeneAssembler, output_path: str):
        self.gene_assembler = gene_assembler
        self.output_path = output_path

    def generate_fas_diversity_among_genes_boxplot(self) -> None:
        dataframe: pandas.DataFrame = self.generate_fas_diversity_among_genes()
        fig = px.box(dataframe, x="transcript_count", y="fas_score")
        fig.show()
        # fig.write_image(os.path.join(self.output_path, "fas_diversity_among_genes.png"))

    def generate_fas_diversity_among_genes(self) -> pandas.DataFrame:
        group_list: List[int] = [0]
        for entry in [[value] * 5 for value in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]]:
            group_list += entry
        group_list += [50] * 100

        data_dict: Dict[str, List[Any]] = {"transcript_count": [],
                                           "fas_score": []}
        for gene in self.gene_assembler.get_genes():
            count: int = 0
            fas_scores_count: int = 0
            for transcript1 in gene.get_transcripts():
                if transcript1.get_biotype() != "nonsense_mediated_decay":
                    count += 1
                    for transcript2 in gene.get_transcripts():
                        if transcript1 != transcript2 and transcript2.get_biotype() != "nonsense_mediated_decay":
                            data_dict["fas_score"].append(gene.fas_dict[transcript1.get_id()][transcript2.get_id()])
                            fas_scores_count += 1
            if count == 1:
                continue
            else:
                data_dict["transcript_count"] += [group_list[count]] * fas_scores_count
        return pandas.DataFrame(data_dict)

    def generate_incomplete_fas_distribution(self) -> None:
        dataframe: pandas.DataFrame = self.generate_fas_comparison_dataframe()
        fig = px.histogram(dataframe, x="fas_score", color="complete_status", barmode="group", nbins=20)
        fig.show()

    def generate_fas_comparison_dataframe(self) -> pandas.DataFrame:
        data_dict: Dict[str, List[Any]] = {"fas_score": [],
                                           "tsl_dist": [],
                                           "complete_status": []}
        for gene in self.gene_assembler.get_genes():
            for transcript1 in gene.get_transcripts():
                if transcript1.get_biotype() == "nonsense_mediated_decay":
                    continue
                for transcript2 in gene.get_transcripts():
                    if transcript1 == transcript2 or transcript2.get_biotype() == "nonsense_mediated_decay":
                        continue
                    else:
                        data_dict["fas_score"].append(gene.fas_dict[transcript1.get_id()][transcript2.get_id()])
                        tsl1: int = transcript1.get_transcript_support_level()
                        tsl2: int = transcript2.get_transcript_support_level()
                        data_dict["tsl_dist"].append(abs(tsl1 - tsl2))
                        status1: bool = any([tag in self.incomplete_tags for tag in transcript1.get_tags()])
                        status2: bool = any([tag in self.incomplete_tags for tag in transcript2.get_tags()])
                        if status1 and status2:
                            data_dict["complete_status"].append("both incomplete")
                        elif status1 and not status2 or not status1 and status2:
                            data_dict["complete_status"].append("one incomplete")
                        else:
                            data_dict["complete_status"].append("both complete")
        return pandas.DataFrame(data_dict)

    def generate_tsl_biotype_histogram(self) -> None:
        dataframe: pandas.DataFrame = self.generate_transcript_dataframe()
        fig = px.histogram(dataframe, x='tsl', color='biotype', barmode='group')
        fig.show()
        # fig.write_image(os.path.join(self.output_path, "tsl_biotype_histogram.png"))

    def generate_tsl_complete_status_histogram(self) -> None:
        dataframe: pandas.DataFrame = self.generate_transcript_dataframe()
        fig = px.histogram(dataframe, x='tsl', color='complete_status', barmode='group')
        fig.show()
        # fig.write_image(os.path.join(self.output_path, "tsl_complete_status_histogram.png"))

    def generate_transcript_dataframe(self) -> pandas.DataFrame:
        tags: List[str] = self.gene_assembler.extract_tags()
        data_dict: Dict[str, List[Any]] = {"_id": [],
                                           "transcript_name": [],
                                           "feature": [],
                                           "gene_id": [],
                                           "taxon_id": [],
                                           "tsl": [],
                                           "biotype": [],
                                           "complete_status": []}
        for tag in tags:
            data_dict[tag] = []

        for transcript in self.gene_assembler.get_transcripts():
            data_dict["_id"].append(transcript.get_id())
            data_dict["transcript_name"].append(transcript.get_name())
            data_dict["feature"].append(transcript.get_feature())
            data_dict["gene_id"].append(transcript.get_id_gene())
            data_dict["taxon_id"].append(transcript.get_id_taxon())
            data_dict["tsl"].append(transcript.get_transcript_support_level())
            data_dict["biotype"].append(transcript.get_biotype())
            for tag in tags:
                if tag in transcript.get_tags():
                    data_dict[tag].append(True)
                else:
                    data_dict[tag].append(False)

            if any(tag in self.incomplete_tags for tag in transcript.get_tags()):
                data_dict["complete_status"].append("incomplete")
            else:
                data_dict["complete_status"].append("complete")
        return pandas.DataFrame(data_dict)


class ResultVisualizer:

    def __init__(self, library_path: str):
        self.library_path = library_path

    def simulate_transcript(self, gene_id: str,
                            transcript_1: str,
                            transcript_2: str) -> List[Tuple[float, float]]:
        with open(os.path.join(self.library_path, "fas_data", "fas_index.json"), "r") as f:
            file_name: str = json.load(f)[gene_id]
        with open(os.path.join(self.library_path, "fas_data", "fas_scores", file_name), "r") as f:
            fas_adjacency_matrix: Dict[str, Dict[str, float]] = json.load(f)[gene_id]

        transcript_list = [transcript_1, transcript_2]
        expression_lists = [[0.0] * len(transcript_list) for _ in range(11)]
        for i, value in enumerate(np.linspace(0.0, 1, 11)):
            expression_lists[i][0] = round(value, 1)
            expression_lists[i][1] = round(1.0 - value, 1)

        ewfd_lists = list()
        for expression_list in expression_lists:
            ewfd_lists.append(ResultVisualizer.calculate_ewfd(fas_adjacency_matrix, expression_list, transcript_list))
        rmsd_list = list()
        log2fc_list = list()
        for i, ewfd_vals_1 in enumerate(ewfd_lists):
            for j, ewfd_vals_2 in enumerate(ewfd_lists):
                if expression_lists[i][0] == 0 or expression_lists[j][0] == 0:
                    log2fold_change = float("inf")
                else:
                    log2fold_change: float = math.log(expression_lists[i][0] / expression_lists[j][0], 2)
                rmsd: float = ResultVisualizer.calc_rmsd(ewfd_vals_1, ewfd_vals_2)
                rmsd_list.append(rmsd)
                log2fc_list.append(round(abs(log2fold_change), 2))
                # print("RMSD:", rmsd, "| lg2fc:", log2fold_change)
        zipped_list = sorted(list(set(zip(rmsd_list, log2fc_list))), key=lambda x: x[0])
        #  for entry in zipped_list:
        #      print("RMSD:", entry[0], "| lg2fc:", entry[1])
        return zipped_list

    @staticmethod
    def calc_rmsd(ewfd_1: List[float], ewfd_2: List[float]):
        squared_delta_list: List[float] = list()
        for i, _ in enumerate(ewfd_1):
            squared_delta_list.append((ewfd_1[i] - ewfd_2[i]) ** 2)
        return math.sqrt(sum(squared_delta_list)) / len(ewfd_1)

    @staticmethod
    def calculate_ewfd(fas_adjacency_matrix: Dict[str, Dict[str, float]],
                       rel_expressions: List[float],
                       transcript_ids: List[str]) -> List[float]:
        ewfd_list: List[float] = [0.0] * len(transcript_ids)
        for s, seed_id in enumerate(transcript_ids):
            for q, query_id in enumerate(transcript_ids):
                ewfd_list[s] += rel_expressions[q] * fas_adjacency_matrix[seed_id][query_id]
        ewfd_list = [round(1 - movement_value, 4) for movement_value in ewfd_list]
        return ewfd_list

    def plot_rmsd_distribution(self, result_directory: str,
                               inclusion_count: int,
                               max_rmsd_inclusion: List[str],
                               slight_shift_inclusion: List[str],
                               included_transcript_pairs: List[List[str]],
                               inclusion_synonym: List[str],
                               inclusion_color: List[List[str]],
                               outfile: str):

        max_rmsd_inclusion_stats: List[Tuple[float, str, str]] = list()  # max_rmsd, label text, color
        for i, gene_id in enumerate(max_rmsd_inclusion):
            zipped_list: List[Tuple[float, float]] = self.simulate_transcript(gene_id,
                                                                              included_transcript_pairs[i][0],
                                                                              included_transcript_pairs[i][1])
            label: str = "{0} Max RMSD".format(inclusion_synonym[i])
            max_rmsd_inclusion_stats.append((max([rmsd for rmsd, _ in zipped_list]),
                                             label,
                                             inclusion_color[i][0]))

        slight_rmsd_inclusion_stats: List[Tuple[float, str, str]] = list()  # max_rmsd, label text, color
        for i, gene_id in enumerate(slight_shift_inclusion):
            zipped_list: List[Tuple[float, float]] = self.simulate_transcript(gene_id,
                                                                              included_transcript_pairs[i][0],
                                                                              included_transcript_pairs[i][1])
            label: str = "{0} slight expr shift RMSD".format(inclusion_synonym[i])
            slight_rmsd: float = 0.0
            slight_log2fold_change = float("inf")
            for rmsd, log2fold_change in zipped_list:
                if 1 < log2fold_change < slight_log2fold_change:
                    slight_rmsd = rmsd
                    slight_log2fold_change = log2fold_change
            slight_rmsd_inclusion_stats.append((slight_rmsd,
                                                label,
                                                inclusion_color[i][1]))

        #  Import the RMSD ranks.
        rank_entries: List[List[float]] = [[] for _ in range(inclusion_count)]
        for entry in os.listdir(result_directory):
            with open(os.path.join(result_directory, entry), "r") as f:
                for i, line in enumerate(f):
                    if i > inclusion_count:
                        break
                    elif i == 0:
                        continue
                    rank_entries[i-1].append(float(line.split(",")[1]))

        # General setup
        fig, ax = plt.subplots()
        positions = range(1, inclusion_count+1)
        bp = ax.boxplot(rank_entries, positions=positions, showfliers=True, flierprops=dict(marker='.', markersize=4))

        # General axis labels
        ax.set_xlabel('Rank')
        ax.set_ylabel('EWFD RMSD')
        ax.set_title('EWFD RMSD by rank')

        # Get less x labels
        ax.set_xlim(0, inclusion_count+1)
        x_ticks = [1] + list(range(20, inclusion_count-19, 20)) + [inclusion_count]
        ax.set_xticks(x_ticks)
        tick_labels = ['{}'.format(tick) for tick in x_ticks]  # Modify tick labels as desired
        ax.set_xticklabels(tick_labels)

        # Get more y labels.
        num_y_ticks = 10
        y_ticks = np.linspace(0.1, 1, num_y_ticks)
        y_tick_labels = ['{:.1f}'.format(tick) for tick in y_ticks]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)

        # Change median color and make them thicker
        for median in bp['medians']:
            median.set(color='red', linewidth=2)

        ax.grid(True, linestyle='-')

        # Change the colors of the boxes
        for box in bp['boxes']:
            box.set(color="gray")

        for whisker in bp['whiskers']:
            whisker.set(color="gray")

        for cap in bp['caps']:
            cap.set(color="gray")

        # Add the lines indicating the RMSD of the interesting candidates.
        for rmsd, label, color in max_rmsd_inclusion_stats:
            ax.axhline(rmsd, color=color, linestyle='--', label=label)

        for rmsd, label, color in slight_rmsd_inclusion_stats:
            ax.axhline(rmsd, color=color, linestyle='--', label=label)

        # Add a legend
        ax.legend()

        plt.savefig(outfile, format='svg')


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        action="store",
                        help="Path to directory containing result csv files.")
    parser.add_argument("-l",
                        "--library",
                        type=str,
                        action="store",
                        help="Path to library.")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        action="store",
                        help="Name of the output file.")
    argument_dict: Dict[str, str] = vars(parser.parse_args())

    result_visualizer = ResultVisualizer(argument_dict["library"])
    result_visualizer.plot_rmsd_distribution(argument_dict["input"], 200,
                                             ["ENSG00000184047"],
                                             ["ENSG00000184047"],
                                             [["ENSP00000411638", "ENSP00000320343"]],
                                             ["DIABLO"],
                                             [["green", "cyan"]],
                                             argument_dict["outfile"])


if __name__ == "__main__":
    main()
