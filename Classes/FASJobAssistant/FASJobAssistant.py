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

import json
import os
import sys
from typing import Dict, Any, List

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

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

    def __init__(self, pass_path: PassPath, memory: str, partitions: List[str], fas_dir: str):
        self.fas_anno: str = os.path.join(fas_dir, "fas.doAnno")
        self.fas_run: str = os.path.join(fas_dir, "fas.run")
        self.python_path: str = sys.executable

        self.lib_pass_path: PassPath = pass_path

        self.memory: str = memory
        self.partitions: List[str] = partitions
        self.phyloprofile_path: str = self.lib_pass_path["transcript_ids"]

    def make_fas_run_jobs(self):
        gene_count: int = FASJobAssistant.make_gene_txt(self.lib_pass_path)

        partitions = ",".join(self.partitions)

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
            with open(fas_lib.get_config("slurm_path") + "gene_ids{0}.txt".format(str(i)), "w") as gene_chunk:
                gene_chunk.write("\n".join(output_ids))
            if lcr_flag:
                outdir = fas_lib.get_config("fas_buffer_lcr_path")
            elif tmhmm_flag:
                outdir = fas_lib.get_config("fas_buffer_tmhmm_path")
            else:
                outdir = fas_lib.get_config("fas_buffer_path")
            output = RAW_SLURM_1.format(python_path,  # 0
                                        FAS_handler_path,  # 1
                                        fas_lib.get_config("self_path"),  # 2
                                        fas_path,  # 3
                                        fas_lib.get_config("species"),  # 4
                                        fas_lib.get_config("isoforms_path"),  # 5
                                        fas_lib.get_config("annotation_path"),  # 6
                                        outdir,  # 7
                                        fas_lib.get_config("tsv_buffer_path"),  # 8
                                        fas_lib.get_config("phyloprofile_ids_path"),  # 9
                                        str(start),  # 10
                                        str(stop),  # 11
                                        str(i),  # 12
                                        fas_lib.get_config("slurm_path"),  # 13
                                        partitions,  # 14
                                        mem_per_cpu)  # 15
            output_2 = RAW_SLURM_2.format(python_path,  # 0
                                          FAS_handler_path,  # 1
                                          fas_lib.get_config("self_path"))  # 2
            if lcr_flag or tmhmm_flag:
                if lcr_flag:
                    mod_path = fas_lib.get_config("lcr_path")
                    job_name = "lcr_FAS_job{0}.job"
                elif tmhmm_flag:
                    mod_path = fas_lib.get_config("tmhmm_path")
                    job_name = "tmhmm_FAS_job{0}.job"
                output = output + "-d " + mod_path + " "
            else:
                job_name = "FAS_job{0}.job"
            output = output + output_2
            with open(fas_lib.get_config("slurm_path") + job_name.format(str(i)), "w") as f:
                f.write(output)

    def make_fas_do_anno_jobs(self):
        pass

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
                                                    "--dir_fas", "--mode_fas", "--outdir"],
                                                   [str, str, str, str, str, str],
                                                   ["store", "store", "store", "store", "store", "store"],
                                                   [1, 1, "*", 1, 1, 1],
                                                   ["Path to a config file of a library.",
                                                    "Required memory specified in the SLURM script.",
                                                    "Set of partitions to be used for the job.",
                                                    "Path to directory containing the fas.doAnno and fas.run binary.",
                                                    "FAS mode to run. Either 'run' or 'doAnno'.",
                                                    "Directory to which the jobs will be saved."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict['Lib_dir'] = argument_dict['Lib_dir'][0]
    argument_dict['memory'] = argument_dict['memory'][0]
    argument_dict['dir_fas'] = argument_dict['dir_fas'][0]
    argument_dict['mode_fas'] = argument_dict['mode_fas'][0]
    argument_dict['outdir'] = argument_dict['outdir'][0]

    with open(os.path.join(argument_dict['Lib_dir'], "paths.json"), "r") as f:
        lib_pass_path: PassPath = json.load(f)

    fas_job_assist: FASJobAssistant = FASJobAssistant(lib_pass_path,
                                                      argument_dict['memory'],
                                                      argument_dict['partitions'],
                                                      argument_dict['fas_dir'])

    if argument_dict["mode_fas"] == "run":
        fas_job_assist.make_fas_run_jobs()
    elif argument_dict["dir_fas"]:
        fas_job_assist.make_fas_do_anno_jobs()
    else:
        print("Chosen FAS mode not recognised. No jobs were generated.")


if __name__ == "__main__":
    main()
