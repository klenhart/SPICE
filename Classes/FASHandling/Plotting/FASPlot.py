#!/bin/env python
from typing import List

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

    def generate_transcript_count_barplot(self):
        data: List[List[int]] = self.fas_gene_assembler.get_protein_coding_transcript_counts()

        labels = data[0][:25]
        values = data[1][:25]

        fig, ax = plt.subplots()

        ax.set_xticks(np.linspace(0, 25, 11))
        ax.set_yticks(np.linspace(0, 6000, 20))

        ax.bar(labels, values)

        ax.set_xlabel('Transcript Count')
        ax.set_ylabel('Gene Count')
        ax.set_title('Distribution of transcript counts among given genes.')

        # Show the plot
        plt.show()

    def generate_fas_score_density_plot(self):
        data_0 = self.fas_gene_assembler.get_all_fas_values(1)
        data_1 = self.fas_gene_assembler.get_all_fas_values(2)
        data_2 = self.fas_gene_assembler.get_all_fas_values(3)
        data_3 = self.fas_gene_assembler.get_all_fas_values(4)
        data_4 = self.fas_gene_assembler.get_all_fas_values(5)
        data_5 = self.fas_gene_assembler.get_all_fas_values(6)
        data_6 = self.fas_gene_assembler.get_all_fas_values()

        kde0 = gaussian_kde(data_0)
        kde1 = gaussian_kde(data_1)
        kde2 = gaussian_kde(data_2)
        kde3 = gaussian_kde(data_3)
        kde4 = gaussian_kde(data_4)
        kde5 = gaussian_kde(data_5)
        kde6 = gaussian_kde(data_6)

        x_grid = np.linspace(0, 1, 100)
        fig, ax = plt.subplots()

        ax.plot(x_grid, kde0(x_grid), label="Only Gencode Basic")
        ax.plot(x_grid, kde1(x_grid), label="TSL >= 1")
        ax.plot(x_grid, kde2(x_grid), label="TSL >= 2")
        ax.plot(x_grid, kde3(x_grid), label="TSL >= 3")
        ax.plot(x_grid, kde4(x_grid), label="TSL >= 4")
        ax.plot(x_grid, kde5(x_grid), label="TSL >= 5")
        ax.plot(x_grid, kde6(x_grid), label="All transcripts")

        ax.legend()

        ax.set_yticks(list(range(0, 14)))
        ax.set_xlabel('FAS score')
        ax.set_ylabel('Density')
        ax.set_title('FAS Score Density Plot')

        ax.grid(True)

        # Show the plot
        plt.show()


def main():
    fas_gene_assembler = FASGeneAssembler("human", "NCBI9606")
    fas_gene_assembler.load("C:/Users/chris/Desktop/git/root/extract_with_fas_no_empty_genes.json")
    fas_plot = FASPlot(fas_gene_assembler)
    #fas_plot.generate_transcript_count_barplot()
    fas_plot.generate_fas_score_density_plot()


if __name__ == "__main__":
    main()