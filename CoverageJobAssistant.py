#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  CoverageJobAssistant is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  CoverageJobAssistant is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import os
from typing import Dict, Any, List
from pathlib import Path

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

JOB_ARRAY: str = """#!/bin/bash
#SBATCH --partition={0}
#SBATCH --cpus-per-task=8
#SBATCH --job-name="coverage"
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-{1}

alignment_job=$(awk FNR==$SLURM_ARRAY_TASK_ID "{2}")

$alignment_job

"""

RAW_SCRIPT: str = "{0} -eB -G {1} -o {2}.gtf {3} -p 8"


class CoverageJobAssistant:

    def __init__(self, aligner_path: str, out_path: str, alignment_id: str, annotation_id: str):
        self.bam_file_path: str = os.path.join(out_path, alignment_id + ".bam")
        self.gtf_file_path: str = os.path.join(out_path, annotation_id + ".gtf")
        self.path_prefix: str = os.path.join(out_path, "coverage_" + alignment_id)
        self.out_file_path: str = os.path.join(self.path_prefix, alignment_id)

        self.command = RAW_SCRIPT.format(aligner_path,
                                         self.gtf_file_path,
                                         self.out_file_path,
                                         self.bam_file_path)

    def __str__(self):
        return self.command


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--input_data", "--aligner_path", "--partitions",
                                                    "--out_path"],
                                                   [str, str, str, str],
                                                   ["store", "store", "store", "store"],
                                                   [1, 1, "*", 1],
                                                   ["Path to directory containing collected data.",
                                                    "Path to the aligner to be used.",
                                                    "Set of partitions to be used for the jobs.",
                                                    "Path to directory that shall contain the 'align_list.txt'."
                                                    "Directory to which the jobs will be saved."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict['input_data'] = argument_dict['input_data'][0]
    argument_dict['aligner_path'] = argument_dict['aligner_path'][0]
    argument_dict['out_path'] = argument_dict['out_path'][0]

    job_list: List[str] = list()

    with open(os.path.join(argument_dict["input_data"], "experiment_list.txt"), "r") as f:
        experiment_id_list: List[str] = f.read().split("\n")

    coverage_list_path: str = os.path.join(argument_dict["input_data"], "coverage_list.txt")
    if not Path(coverage_list_path).exists():
        Path(coverage_list_path).touch()

    for exp_id in experiment_id_list:
        experiment_directory: str = os.path.join(argument_dict["input_data"], exp_id)

        with open(os.path.join(experiment_directory, "annotation_list.txt"), "r") as f:
            annotation_list: List[str] = f.read().split("\n")

        with open(os.path.join(experiment_directory, "alignment_list.txt"), "r") as f:
            alignment_list: List[str] = f.read().split("\n")

        for i, annotation_id in enumerate(annotation_list):
            alignment_id = alignment_list[i]
            coverage_job: CoverageJobAssistant = CoverageJobAssistant(argument_dict["aligner_path"],
                                                                      experiment_directory,
                                                                      alignment_id,
                                                                      annotation_id)
            Path(coverage_job.path_prefix).mkdir(parents=True, exist_ok=True)
            job_list.append(str(coverage_job))

    alignment_job_list_path: str = os.path.join(argument_dict['out_path'], "coverage_job_list.txt")
    with open(alignment_job_list_path, "w") as f:
        f.write("\n".join(job_list))

    alignment_job_array_path: str = os.path.join(argument_dict['out_path'], "coverage_job_array.job")
    alignment_job_array: str = JOB_ARRAY.format(",".join(argument_dict["partitions"]),
                                                str(len(job_list)),
                                                alignment_job_list_path)
    with open(alignment_job_array_path, "w") as f:
        f.write(alignment_job_array)


if __name__ == "__main__":
    main()
