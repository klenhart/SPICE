#!/bin/env python
import os.path
from typing import List, Type

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  FASPlot is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASPlot is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from Classes.FASHandling.FASGeneAssembler import FASGeneAssembler


class FASPlot:

    def __init__(self, fas_gene_assembler: FASGeneAssembler, output_dir: str = ""):
        self.fas_gene_assembler: FASGeneAssembler = fas_gene_assembler
        self.output_dir: str = output_dir

    def generate_fas_distribution_boxplot(self):
        fig, ax = plt.subplots()
        ax.boxplot(self.fas_gene_assembler.get_fas_grouped_by_tsl())

        # Set axis labels and title
        ax.set_xticklabels(['TSL 1 (or basic)', 'TSL 2', 'TSL 3', 'TSL 4', 'TSL 5', 'NO TSL'])
        ax.set_title('FAS distribution among different TSLs.')

        # Show the plot
        plt.show()

    def generate_transcript_count_barplot(self) -> None:
        data_4: List[List[int]] = self.fas_gene_assembler.get_protein_coding_transcript_counts(5)
        data_6: List[List[int]] = self.fas_gene_assembler.get_protein_coding_transcript_counts()

        labels_4 = data_4[0][:25]
        values_4 = data_4[1][:25]

        labels_6 = data_6[0][:25]
        values_6 = data_6[1][:25]

        fig, ax = plt.subplots()

        ax.bar(labels_4, values_4, color="red", align='center', width=0.9,
               label="with TSL <= 5")
        ax.bar(labels_6, values_6, color="blue", align='center', width=0.6,
               label="all protein coding")

        ax.legend()

        ax.set_xlabel('Transcript Count')
        ax.set_ylabel('Gene Count')
        ax.set_title('Distribution of transcript counts among given genes.')

        plt.savefig(os.path.join(self.output_dir, "transcript_count_bars.png"))
        plt.clf()

    def generate_transcript_count_nmd_barplot(self) -> None:
        data_6: List[List[int]] = self.fas_gene_assembler.get_protein_coding_transcript_counts(7, False)
        data_6_nmd: List[List[int]] = self.fas_gene_assembler.get_protein_coding_transcript_counts(7, True)

        labels_6 = data_6[0][:25]
        values_6 = data_6[1][:25]

        labels_6_nmd = data_6_nmd[0][:25]
        values_6_nmd = data_6_nmd[1][:25]

        fig, ax = plt.subplots()

        ax.bar(labels_6, values_6, color="red", align='center', width=0.9,
               label="only protein coding")
        ax.bar(labels_6_nmd, values_6_nmd, color="blue", align='center', width=0.6,
               label="protein coding and NMD")

        ax.legend()

        ax.set_xlabel('Transcript Count')
        ax.set_ylabel('Gene Count')
        ax.set_title('Distribution of transcript counts among given genes.')

        plt.savefig(os.path.join(self.output_dir, "transcript_count_nmd_bars.png"))
        plt.clf()

    def generate_fas_score_histogram(self) -> None:
        data_4 = self.fas_gene_assembler.get_all_fas_values(5)
        data_6 = self.fas_gene_assembler.get_all_fas_values()

        plt.grid(True)

        plt.hist(data_6,
                 bins=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25,
                       0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                       0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                       0.9, 0.95, 1.0, 1.05],
                 align='left',
                 rwidth=0.8,
                 color='red',
                 label="all protein coding")

        plt.hist(data_4,
                 bins=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25,
                       0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                       0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                       0.9, 0.95, 1.0, 1.05],
                 align='left',
                 rwidth=0.8,
                 color='blue',
                 label="with TSL <= 5")

        plt.legend()
        plt.xlabel('FAS score')
        plt.ylabel('Frequency')
        plt.title('FAS Score Histogram')

        plt.savefig(os.path.join(self.output_dir, "fas_score_hist.png"))
        plt.clf()

    def generate_fas_score_nmd_histogram(self) -> None:
        data_6 = self.fas_gene_assembler.get_all_fas_values(7, False)
        data_6_nmd = self.fas_gene_assembler.get_all_fas_values(7, True)

        plt.grid(True)

        plt.hist(data_6_nmd,
                 bins=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25,
                       0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                       0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                       0.9, 0.95, 1.0, 1.05],
                 align='left',
                 rwidth=0.9,
                 color='blue',
                 label="protein coding and nmd")

        plt.hist(data_6,
                 bins=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25,
                       0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                       0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                       0.9, 0.95, 1.0, 1.05],
                 align='left',
                 rwidth=0.6,
                 color='red',
                 label="only protein coding")

        plt.yticks([500000, 1000000, 1500000])

        plt.legend()

        plt.xlabel('FAS score')
        plt.ylabel('Frequency')
        plt.title('FAS Score (with nonsense mediated decay transcripts) Histogram')

        plt.savefig(os.path.join(self.output_dir, "fas_score_nmd_hist.png"))
        plt.clf()

    def generate_nonsense_mediated_decay_pieplot(self) -> None:
        sizes: List[int] = self.fas_gene_assembler.get_nonsense_mediated_decay_gene_ratio()

        fig, ax = plt.subplots()
        labels = ['with NMD', 'without NMD']
        colors = ['green', 'red']

        ax.pie(sizes, labels=labels, colors=colors, autopct='%.1f%%')
        ax.set_title('Ratio of genes with nonsense-mediated decay to genes without')

        plt.savefig(os.path.join(self.output_dir, "nmd_ratio_pie.png"))
        plt.clf()


def main():
    fas_gene_assembler = FASGeneAssembler("human", "NCBI9606")
    fas_gene_assembler.load("C:/Users/chris/Desktop/git/root/extract_without_fas.json")
    fas_plot = FASPlot(fas_gene_assembler, "C:/Users/chris/Desktop/AKE/ProgReps/11")
    fas_plot.generate_transcript_count_barplot()
    fas_plot.generate_fas_score_histogram()
    fas_plot.generate_nonsense_mediated_decay_pieplot()
    fas_plot.generate_fas_score_nmd_histogram()
    fas_plot.generate_transcript_count_nmd_barplot()


if __name__ == "__main__":
    main()
