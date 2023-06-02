#!/bin/env python
import gzip
import os
import shutil
from pathlib import Path
#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  EncodeScrapy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EncodeScrapy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from typing import Dict, Any, List

import requests
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

URL: str = "https://www.encodeproject.org/batch_download/?type=Experiment&@id=/experiments/{0}/&files.processed=true"


def extract_id(link: str) -> str:
    split_link: List[str] = link.split("/")
    if split_link[-1] == "":
        return split_link[-2]
    else:
        return split_link[-1]


def extract_experiment_ids_from_links(experiment_file_path: str) -> List[str]:
    return [extract_id(line) for line in open_url_file(experiment_file_path)]


def open_url_file(in_path: str) -> List[str]:
    with open(in_path, "r") as f:
        file = f.read()
    return file.split("\n")


def download_gtf_gz(url: str, file_id: str, exp_id: str, out_path: str):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        Path(os.path.join(out_path, exp_id)).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(out_path, exp_id, file_id + ".gtf.gz"), 'wb') as file:
            response.raw.decode_content = True
            shutil.copyfileobj(response.raw, file)


def download_annotation(exp_id: str, out_path: str):
    url: str = URL.format(exp_id)
    response_files = requests.get(url)
    if response_files.status_code == 200:
        download_link_list: List[str] = response_files.content.decode('utf-8').split("\n")
        download_link_list = [link for link in download_link_list if link.endswith("gtf.gz")]
        for download_link in download_link_list:
            file_id: str = download_link.split("@@download/")[1].split(".")[0]
            download_gtf_gz(download_link, file_id, exp_id, out_path)
    else:
        print(f"Error: {response_files.status_code} - {response_files.text}")
    #
    ## output_type == "alignments"
    #url = "https://www.encodeproject.org/files/ENCFF203MFP/?format=json"
#
    #response = requests.get(url)
#
    #if response.status_code == 200:
    #    text = response.json()
    #    for entry in text:
    #        print(entry, ":", text[entry])


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--mode", "--input", "--outdir"],
                                                   [str, str, str],
                                                   ["store", "store", "store"],
                                                   [1, 1, 1],
                                                   ["""Either [annotation] or [alignment]
                                                    depending on what shall get scraped.""",
                                                    "Path to the file containing the encode links.",
                                                    "Path to the directory that shall contain the downloaded files."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()
    argument_dict["mode"] = argument_dict["mode"][0]
    argument_dict["input"] = argument_dict["input"][0]
    argument_dict["outdir"] = argument_dict["outdir"][0]

    if argument_dict["mode"] == "alignment":
        pass
    elif argument_dict["mode"] == "annotation":
        experiment_ids = extract_experiment_ids_from_links(argument_dict["input"])
        # experiment_urls: List[str] = open_url_file(argument_dict["input"])
        for exp_id in experiment_ids:
            download_annotation(exp_id, argument_dict["outdir"])
            break

    else:
        print("Mode", argument_dict["mode"], "not recognised. Shutting down.")


if __name__ == "__main__":
    main()
