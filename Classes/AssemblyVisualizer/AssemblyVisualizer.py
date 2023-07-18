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
from scipy.interpolate import make_interp_spline

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
                            transcript_2: str) -> float:
        with open(os.path.join(self.library_path, "fas_data", "fas_index.json"), "r") as f:
            file_name: str = json.load(f)[gene_id]
        with open(os.path.join(self.library_path, "fas_data", "fas_scores", file_name), "r") as f:
            fas_adjacency_dict: Dict[str, Dict[str, float]] = json.load(f)[gene_id]

        fas_adjacency_matrix: np.array = ResultVisualizer.distance_dict_to_matrix(fas_adjacency_dict)

        all_transcripts: List[str] = list(fas_adjacency_dict.keys())
        expression_vector_1 = np.zeros(len(all_transcripts))
        expression_vector_2 = np.zeros(len(all_transcripts))

        transcript_1_index: int = all_transcripts.index(transcript_1)
        transcript_2_index: int = all_transcripts.index(transcript_2)

        expression_vector_1[transcript_1_index] = 1.0
        expression_vector_2[transcript_2_index] = 1.0

        expressed_indices = [transcript_1_index, transcript_2_index]

        diff = np.dot(fas_adjacency_matrix, expression_vector_1) - np.dot(fas_adjacency_matrix, expression_vector_2)
        expressed_diff_list = list()
        for i in expressed_indices:
            expressed_diff_list.append(diff[i])
        expressed_diff = np.array(expressed_diff_list)
        return np.sqrt(np.sum(expressed_diff ** 2) / 2)

    @staticmethod
    def distance_dict_to_matrix(distance_dict: Dict[str, Dict[str, float]]) -> np.array:
        raw_matrix: List[List[float]] = list()
        for protein_id_outer in distance_dict.keys():
            raw_matrix.append(list())
            for protein_id_inner in distance_dict[protein_id_outer].keys():
                raw_matrix[-1].append(1 - distance_dict[protein_id_outer][protein_id_inner])
        return np.array(raw_matrix)

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

    def plot_rmsd_distribution(self,
                               result_directory: str,
                               rank_count: int,
                               max_rmsd_genes: List[str],
                               max_rmsd_path: str,
                               max_rmsd_category: str,
                               max_rmsd_gene_colors: List[str],
                               sim_switch_genes: List[str],
                               sim_switch_transcripts: List[List[str]],
                               sim_switch_synonyms: List[List[str]],
                               sim_switch_color: List[str],
                               gene_synonyms: List[str],
                               filter_tag: str,
                               outfile: str,
                               library_name: str,
                               switch_label_flag: bool):

        # Extract the max rmsd inclusion file here.
        max_rmsd_dict: Dict[str, float] = ResultVisualizer.extract_max_rmsd_file(max_rmsd_path,
                                                                                 max_rmsd_category)

        max_rmsd_stats: List[Tuple[float, str, str]] = ResultVisualizer.make_max_rmsd_stats(max_rmsd_genes,
                                                                                            max_rmsd_gene_colors,
                                                                                            gene_synonyms,
                                                                                            max_rmsd_dict)

        sim_switch_rmsds: List[float] = list()
        # Simulate a transcript
        for i, gene_id in enumerate(sim_switch_genes):
            sim_switch_rmsds.append(self.simulate_transcript(gene_id,
                                                             sim_switch_transcripts[i][0],
                                                             sim_switch_transcripts[i][1]))

        sim_switch_stats: List[Tuple[float, str, str]] = ResultVisualizer.make_sim_switch_stats(sim_switch_rmsds,
                                                                                                sim_switch_color,
                                                                                                gene_synonyms,
                                                                                                sim_switch_synonyms)

        rank_entries, rank_max = ResultVisualizer.import_rmsd_ranks(result_directory, max_rmsd_dict, rank_count)

        # General setup
        fig, ax = plt.subplots()
        positions = range(1, rank_count+1)
        bp = ax.boxplot(rank_entries, positions=positions, showfliers=False, flierprops=dict(marker='.', markersize=1),
                        zorder=2)

        # General axis labels
        ax.set_xlabel('Rank')
        ax.set_ylabel(r'$RMSD_{EWFD}$')
        if filter_tag != "" and library_name != "":
            ax.set_title(r'$RMSD_{EWFD}$ by rank ' + '\n (filter by ' + '{0}) \n ({1} library)'.format(filter_tag,
                                                                                                       library_name))
        elif filter_tag != "":
            ax.set_title(r'$RMSD_{EWFD}$ by rank ' + '\n (filter by ' + '{0})'.format(filter_tag))
        elif library_name != "":
            ax.set_title(r'$RMSD_{EWFD}$ by rank ' + '\n (no filter) \n (' + '{0} library)'.format(library_name))
        else:
            ax.set_title(r'$RMSD_{EWFD}$ by rank' + '\n (no filter)')

        # Get less x labels
        ax.set_xlim(0, rank_count+1)
        x_ticks = [1] + list(range(20, rank_count-19, 20)) + [rank_count]
        ax.set_xticks(x_ticks)
        tick_labels = ['{}'.format(tick) for tick in x_ticks]  # Modify tick labels as desired
        ax.set_xticklabels(tick_labels)

        # Get more y labels.
        num_y_ticks = 11
        y_ticks = np.linspace(0.0, 1, num_y_ticks)
        y_tick_labels = ['{:.1f}'.format(tick) for tick in y_ticks]
        ax.set_ylim(-0.1, 1.1)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)

        # Change median color and make them thicker
        for median in bp['medians']:
            median.set(color='royalblue', linewidth=2, zorder=3)

        ax.grid(True, linestyle='-', zorder=1, color="lightgrey")

        # Change the colors of the boxes
        for box in bp['boxes']:
            box.set(color="darkgray")

        for whisker in bp['whiskers']:
            whisker.set(color="darkgray")

        for cap in bp['caps']:
            cap.set(color="darkgray")

        # Add the lines indicating the RMSD of the interesting candidates.
        for rmsd, label, color in max_rmsd_stats:
            ax.axhline(rmsd, color=color, linestyle='--', zorder=8)
            if switch_label_flag:
                ax.text(rank_count / 2, rmsd + 0.03, label, ha='left', va='center', color=color, fontsize=8, zorder=5)
            else:
                ax.text(rank_count/2, rmsd - 0.03, label, ha='left', va='center', color=color, fontsize=8, zorder=5)

        for rmsd, label, color in sim_switch_stats:
            ax.axhline(rmsd, color=color, linestyle='--', zorder=7)
            ax.text(rank_count/2.5, rmsd + 0.03, label, ha='left', va='center', color=color, fontsize=8, zorder=6)

        # Smooth the line using interpolation
        x_new = np.linspace(1, rank_count, 13)
        spl = make_interp_spline(range(1, rank_count+1), rank_max, k=3)
        line_smooth = spl(x_new)
        ax.text(rank_count/2 - 20, min(line_smooth) - 0.03,
                r"Regression of mean $RMSD_{EWFD}$ Max by rank", ha='left', va='center', color="royalblue",
                fontsize=8, zorder=6)

        ax.plot(np.linspace(1, rank_count, 13), line_smooth, color='royalblue', label=r'Mean $RMSD_{EWFD}$ Max')

        plt.savefig(outfile, format='svg')
        print(outfile)

    @staticmethod
    def import_rmsd_ranks(result_directory: str, max_rmsd_dict: Dict[str, float], rank_count: int):
        max_entries: List[List[float]] = [[] for _ in range(rank_count)]
        rank_entries: List[List[float]] = [[] for _ in range(rank_count)]

        for entry in os.listdir(result_directory):
            with open(os.path.join(result_directory, entry), "r") as f:
                for i, line in enumerate(f):
                    if i > rank_count:
                        break
                    elif i == 0:
                        continue
                    if line.split(",")[0] in max_rmsd_dict.keys():
                        max_entries[i - 1].append(max_rmsd_dict[line.split(",")[0]])
                    rank_entries[i - 1].append(float(line.split(",")[1]))

        max_rmsd_avgs: List[float] = [sum(entry)/len(entry) if len(entry) > 0 else 0.0 for entry in max_entries]
        return rank_entries, max_rmsd_avgs

    @staticmethod
    def make_max_rmsd_stats(max_rmsd_genes: List[str],
                            max_rmsd_gene_colors: List[str],
                            gene_synonyms: List[str],
                            max_rmsd_dict: Dict[str, float]):
        max_rmsd_stats: List[Tuple[float, str, str]] = list()
        for i, gene_id in enumerate(max_rmsd_genes):
            max_rmsd_stats.append((max_rmsd_dict[gene_id],
                                  gene_synonyms[i] + r": Max $RMSD_{EWFD}$",
                                  max_rmsd_gene_colors[i]))
        return max_rmsd_stats

    @staticmethod
    def make_sim_switch_stats(sim_switch_rmsds,
                              sim_switch_color,
                              gene_synonyms,
                              sim_switch_synonyms):
        sim_switch_stats: List[Tuple[float, str, str]] = list()
        for i, rmsd in enumerate(sim_switch_rmsds):
            sim_switch_stats.append((rmsd,
                                     "{0} switch: {1} to {2} ({3})".format(gene_synonyms[i],
                                                                           sim_switch_synonyms[i][0],
                                                                           sim_switch_synonyms[i][1],
                                                                           str(round(rmsd, 2))),
                                     sim_switch_color[i]))
        return sim_switch_stats

    @staticmethod
    def extract_interquartile_range(max_rmsd_dict: Dict[str, float]) -> Tuple[float, float]:
        all_max_rmsds = list(max_rmsd_dict.values())
        all_max_rmsds.sort(reverse=True)
        upper_end: float = all_max_rmsds[int(len(all_max_rmsds) * 3/4)]
        lower_end: float = all_max_rmsds[int(len(all_max_rmsds) * 1/4)]
        return lower_end, upper_end

    @staticmethod
    def extract_max_rmsd_file(max_rmsd_path: str, category: str) -> Dict[str, float]:
        max_rmsd_dict: Dict[str, float] = dict()
        with open(max_rmsd_path, "r") as f:
            file = f.read()
            line_list: List[str] = file.split("\n")
            categories: List[str] = line_list[0].split("\t")
            category_index = categories.index(category)
            for line in line_list[1:]:
                entries = line.split("\t")
                max_rmsd_dict[entries[0]] = float(entries[category_index])
        return max_rmsd_dict


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
    parser.add_argument("-m",
                        "--maximum",
                        type=str,
                        action="store",
                        help="Filepath to the tsv containing the maximum RMSD for each gene.")
    parser.add_argument("-c",
                        "--category_max_rmsd",
                        type=str,
                        action="store",
                        help="""Which category of the max RMSD file shall be extracted. 
                        [all], [no_incomplete],	[no_non_coding] or [no_both].""")
    parser.add_argument("-f",
                        "--filter_tag",
                        type=str,
                        action="store",
                        help="""Filter tag that shall be added to the title""")
    parser.add_argument("-n",
                        "--name_library",
                        type=str,
                        action="store",
                        help="""What is the name of the library?""")
    parser.add_argument("-s",
                        "--switch_max",
                        action="store_true",
                        help="""Switch the max RMSD label to the top.""")
    argument_dict: Dict[str, Any] = vars(parser.parse_args())

    if argument_dict["filter_tag"] is None:
        argument_dict["filter_tag"] = ""
    if argument_dict["name_library"] is None:
        argument_dict["name_library"] = ""

    result_visualizer = ResultVisualizer(argument_dict["library"])
    result_visualizer.plot_rmsd_distribution(result_directory=argument_dict["input"],
                                             rank_count=200,
                                             max_rmsd_genes=["ENSG00000184047"],
                                             max_rmsd_path=argument_dict["maximum"],
                                             max_rmsd_category=argument_dict["category_max_rmsd"],
                                             gene_synonyms=["DIABLO"],
                                             max_rmsd_gene_colors=["red"],
                                             sim_switch_genes=["ENSG00000184047"],
                                             sim_switch_transcripts=[["ENSP00000442360", "ENSP00000320343"]],
                                             sim_switch_synonyms=[["w/o TM", "w/ TM"]],
                                             sim_switch_color=["red"],
                                             outfile=argument_dict["outfile"],
                                             filter_tag=argument_dict["filter_tag"],
                                             library_name=argument_dict["name_library"],
                                             switch_label_flag=argument_dict["switch_max"])


if __name__ == "__main__":
    main()
