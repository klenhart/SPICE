#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  FASJobAssistant is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASJobAssistant is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os

from Classes.PassPath.PassPath import PassPath

RAW_SCRIPT = """#!/bin/bash
#SBATCH --partition={14}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={15}G
#SBATCH --job-name="fas_{4}{12}"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
#SBATCH --array={10}-{11}
gene=$(awk FNR==$SLURM_ARRAY_TASK_ID "{13}gene_ids{12}.txt")
{0} {1} \
-m \
-g $gene \
-c {2} \
&& \
{3} \
--seed {5} \
--query {5} \
--annotation_dir {6} \
--out_dir {7} \
--bidirectional \
--pairwise {8}$gene.tsv \
--out_name $gene \
--tsv \
--phyloprofile {9} \
--empty_as_1 \
"""

RAW_SLURM_2 = """; \
{0} {1} \
-r \
-g $gene \
-c {2}"""


class FASJobAssistant:

    def __init__(self, pass_path: PassPath, fas_dir: str):
        self.fas_anno: str = os.path.join(fas_dir, "fas.doAnno")
        self.fas_run: str = os.path.join(fas_dir, "fas.run")
