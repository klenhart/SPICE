#!/bin/env python
import os

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
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

from typing import Dict, List, Any

import plotly.express as px

import pandas
import numpy as np


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
        # fig.write_image(os.path.join(self.output_path, "incomplete_fas_distrib.png"))

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


def main():
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens", "9606")
    gene_assembler.load("C:/Users/chris/Desktop/git/fade_lib_homo_sapiens_107/transcript_data/transcript_set.json")
    assembly_visualizer: AssemblyVisualizer = AssemblyVisualizer(gene_assembler,
                                                                 "C:/Users/chris/Desktop/AKE/ProgReps/12")
    assembly_visualizer.generate_tsl_biotype_histogram()
    assembly_visualizer.generate_tsl_complete_status_histogram()
    assembly_visualizer.generate_incomplete_fas_distribution()
    assembly_visualizer.generate_fas_diversity_among_genes_boxplot()


if __name__ == "__main__":
    main()
