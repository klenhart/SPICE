#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  FASTools is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FASTools is distributed in the hope that it will be useful,
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
import sys
from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

RAW_SCRIPT_1 = """#!/bin/bash
#SBATCH --partition={10}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={11}G
#SBATCH --job-name="fas_{9}"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
#SBATCH --array={7}-{8}
gene=$(awk FNR==$SLURM_ARRAY_TASK_ID "{5}/gene_ids{9}.txt")
{0} {1} \\
--pairings_path {12} \\
--mode unpack \\
--gene_id $gene \\
--out_dir {5} \\
&& \\
{2} \\
--seed {3} \\
--query {3} \\
--annotation_dir {4} \\
--out_dir {5} \\
--bidirectional \\
--pairwise {5}/$gene.tsv \\
--out_name $gene \\
--tsv \\
--phyloprofile {6} \\
--empty_as_1 \\
--featuretypes {13} \\
"""

RAW_SCRIPT_2 = """; \\
{0} {1} \\
--mode concat \\
--gene_id $gene \\
--out_dir {2} \\
--anno_dir {3} \\
; \\
{0} {1} \\
--mode delete \\
--gene_id $gene \\
--out_dir {2}"""

RAW_SCRIPT_3 = """#!/bin/bash
#SBATCH --partition={2}
#SBATCH --cpus-per-task=1
#SBATCH --job-name="fas_final"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
{0} {1} \\
--mode integrate \\
--anno_dir {3}
"""

DO_ANNO_SCRIPT = """#!/bin/bash
#SBATCH --partition={0}
#SBATCH --profile=task
#SBATCH --cpus-per-task={1}
#SBATCH --job-name="fasdoAnno"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null

{2} \\
-i {3} \\
-o {4} \\
-t {5} \\
-n annotations \\
--cpus {1} \\
&& \\
{6} \\
{7} \\
-i {8} \\
-o {4} \\
&& \\
{6} \\
{9} \\
-i {8} \\
-o {4}/architectures"""


class FASJobAssistant:

    def __init__(self, pass_path: PassPath, memory: str, partitions: List[str], fas_dir: str, out_dir: str,
                 anno_partitions: List[str], anno_cpus: int, anno_tools: str):
        self.fas_anno: str = os.path.join(fas_dir, "fas.doAnno")
        self.fas_run: str = os.path.join(fas_dir, "fas.run")
        self.out_dir: str = out_dir
        self.python_path: str = sys.executable

        self.spice_path = os.path.abspath(__file__).split("/")[:-1]
        self.fas_result_handler = "/".join((self.spice_path + ["FASResultHandler.py"]))
        self.get_domain_importance = "/".join((self.spice_path + ["get_domain_importance.py"]))
        self.restructure_anno = "/".join((self.spice_path + ["restructure_anno.py"]))
        self.anno_tools: str = anno_tools

        self.lib_pass_path: PassPath = pass_path

        self.memory: str = memory
        self.partitions: str = ",".join(partitions)
        self.anno_partitions: str = ",".join(anno_partitions)
        self.anno_cpus: str = str(anno_cpus)
        self.phyloprofile_path: str = self.lib_pass_path["transcript_ids"]

    def make_fas_run_jobs(self):
        gene_count: int = FASJobAssistant.make_gene_txt(self.lib_pass_path)

        with open(os.path.join(self.lib_pass_path["transcript_data"], "genes.txt"), "r") as f:
            gene_ids: List[str] = f.read().split("\n")

        jobs_ranges = FASJobAssistant.start_stop_range(gene_count, 1000)
        for i, entry in enumerate(jobs_ranges):
            start = 1
            if (entry[1] + 1) % 1000 == 0:
                stop = 1000
            else:
                stop = ((entry[1] + 1) % 1000)
            output_ids = gene_ids[entry[0]:entry[1] + 1]
            with open(os.path.join(self.lib_pass_path["fas_temp"], "gene_ids{0}.txt".format(str(i))),
                      "w") as gene_chunk:
                gene_chunk.write("\n".join(output_ids))

            output = RAW_SCRIPT_1.format(self.python_path,  # 0
                                         self.fas_result_handler,  # 1
                                         self.fas_run,  # 2
                                         os.path.join(self.lib_pass_path["transcript_data"], "annotations.fasta"),  # 3
                                         self.lib_pass_path["fas_data"],  # 4
                                         self.lib_pass_path["fas_temp"],  # 5
                                         self.phyloprofile_path,  # 6
                                         str(start),  # 7
                                         str(stop),  # 8
                                         str(i),  # 9
                                         self.partitions,  # 10
                                         self.memory,  # 11
                                         self.lib_pass_path["transcript_pairings"],  # 12
                                         self.lib_pass_path["fas_annoTools"])  # 13

            output_2 = RAW_SCRIPT_2.format(self.python_path,  # 0
                                           self.fas_result_handler,  # 1
                                           self.lib_pass_path["fas_temp"],  # 2
                                           self.lib_pass_path["fas_data"])  # 3

            job_name = "FAS_job" + str(i) + ".job"

            output = output + output_2

            with open(os.path.join(self.out_dir, job_name), "w") as f:
                f.write(output)
        with open(os.path.join(self.out_dir, "FAS_job_final.job"), "w") as f:
            output: str = RAW_SCRIPT_3.format(self.python_path,
                                              self.fas_result_handler,
                                              self.partitions,
                                              self.lib_pass_path["fas_data"])
            f.write(output)

    def make_fas_do_anno_jobs(self):
        output = DO_ANNO_SCRIPT.format(self.anno_partitions,  # 0 partition
                                       self.anno_cpus,  # 1 cpus
                                       self.fas_anno,  # 2 fas.doAnno binary
                                       self.lib_pass_path["transcript_fasta"],  #  3 transcript_set.fasta
                                       self.lib_pass_path["fas_data"],  # 4 /fas_data/
                                       self.anno_tools,  # 5 /annoTools
                                       self.python_path,   # 6 Python path
                                       self.get_domain_importance, # 7 Get domain importance path
                                       os.path.join(self.lib_pass_path["fas_data"],
                                                    "annotations.json"),  # 8 /fas_data/annotations.json
                                       self.restructure_anno)  # 9 restructure anno path
        with open(os.path.join(self.out_dir, "doAnno.job"), "w") as f:
            f.write(output)

    def make_fas_run_output(self) -> None:
        if not os.path.exists(os.path.join(self.lib_pass_path["fas_data"], "forward.domains")):
            with open(os.path.join(self.lib_pass_path["fas_data"], "forward.domains"), "w") as f:
                f.write("")
        if not os.path.exists(os.path.join(self.lib_pass_path["fas_data"], "reverse.domains")):
            with open(os.path.join(self.lib_pass_path["fas_data"], "reverse.domains"), "w") as f:
                f.write("")
        if not os.path.exists(os.path.join(self.lib_pass_path["fas_data"], "fas.phyloprofile")):
            with open(os.path.join(self.lib_pass_path["fas_data"], "fas.phyloprofile"), "w") as f:
                f.write("geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n")

    @staticmethod
    def make_gene_txt(pass_path: PassPath) -> int:
        with open(pass_path["transcript_pairings"], "r") as f:
            genes_list: List[str] = json.load(f).keys()

        with open(os.path.join(pass_path["transcript_data"], "genes.txt"), "w") as f:
            f.write("\n".join(genes_list))

        return len(genes_list)

    @staticmethod
    def start_stop_range(length: int, chunk_size: int):
        for i in range(0, length, chunk_size):
            yield i, min(i + chunk_size - 1, length)


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--Lib_dir", "--memory", "--partitions",
                                                    "--dir_fas", "--outdir",
                                                    "--anno_partitions", "--Anno_cpus", "--tools_anno"],
                                                   [str, str, str, str, str, str, int, str],
                                                   ["store", "store", "store", "store", "store",
                                                    "store", "store", "store"],
                                                   [1, 1, "*", 1, 1, "*", 1, 1],
                                                   ["Path to a config file of a library.",
                                                    "Required memory specified in the SLURM script.",
                                                    "Set of partitions to be used for the job.",
                                                    "Path to directory containing the fas.doAnno and fas.run binary.",
                                                    "Directory to which the jobs will be saved.",
                                                    "Partitions that shall be used for the annotation job.",
                                                    "Number of CPUs that shall be used for the annotation job.",
                                                    "Path to the directory containing the annotation tools."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict['Lib_dir'] = argument_dict['Lib_dir'][0]
    argument_dict['memory'] = argument_dict['memory'][0]
    argument_dict['dir_fas'] = argument_dict['dir_fas'][0]
    argument_dict['outdir'] = argument_dict['outdir'][0]
    argument_dict['Anno_cpus'] = argument_dict['Anno_cpus'][0]
    argument_dict['tools_anno'] = argument_dict['tools_anno'][0]

    with open(os.path.join(argument_dict['Lib_dir'], "paths.json"), "r") as f:
        lib_pass_path: PassPath = PassPath(json.load(f))

    fas_job_assist: FASJobAssistant = FASJobAssistant(lib_pass_path,
                                                      argument_dict['memory'],
                                                      argument_dict['partitions'],
                                                      argument_dict['dir_fas'],
                                                      argument_dict['outdir'],
                                                      argument_dict['anno_partitions'],
                                                      argument_dict['Anno_cpus'],
                                                      argument_dict['tools_anno'])

    fas_job_assist.make_fas_run_jobs()
    fas_job_assist.make_fas_run_output()
    fas_job_assist.make_fas_do_anno_jobs()


if __name__ == "__main__":
    main()
