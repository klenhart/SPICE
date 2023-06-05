#!/bin/env python
import sys
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

import json
import requests
import os
import shutil
import datetime
from pathlib import Path
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

URL: str = "https://www.encodeproject.org/batch_download/?type=Experiment&@id=/experiments/{0}/&files.processed=true"
EXP_URL: str = "https://www.encodeproject.org/experiments/{0}/"


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


def download_filtered_alignment(exp_id, out_path: str):
    pass


def download_annotation(exp_id: str, out_path: str):
    url: str = URL.format(exp_id)
    response_files = requests.get(url)
    if response_files.status_code == 200:
        download_link_list: List[str] = response_files.content.decode('utf-8').split("\n")
        download_link_list = [link for link in download_link_list if check_link_meta(link, "annotation")]

        for download_link in download_link_list:
            file_id: str = download_link.split("@@download/")[1].split(".")[0]
            download_gtf_gz(download_link, file_id, exp_id, out_path)

        with open(os.path.join(out_path, exp_id, "info.json"), "r") as f:
            info_dict: Dict[str, Any] = json.load(f)
        info_dict["annotation_urls"] = download_link_list
        with open(os.path.join(out_path, exp_id, "info.json"), "w") as f:
            json.dump(info_dict, f, indent=4)
    else:
        print(f"Error: {response_files.status_code} - {response_files.text}")


def check_link_meta(url: str, output_category: str) -> bool:
    if "processed=true" in url or len(url) == 0:
        return False
    else:
        meta_url: str = url.split("@@download")[0] + "?format=json"
        response = requests.get(meta_url)
        if response.status_code == 200:
            meta_json = response.json()
            if meta_json["output_category"] == output_category:
                return True
            else:
                return False
        else:
            print(f"Error: {response.status_code} - {response.text}")

    # output_type == "alignments"
    # url = "https://www.encodeproject.org/files/ID/?format=json"



def make_log_file(exp_id: str, out_path: str):
    directory: str = os.path.join(out_path, exp_id)
    info_file: str = os.path.join(directory, "info.json")
    Path(directory).mkdir(parents=True, exist_ok=True)
    if not Path(info_file).exists():
        info_dict: Dict[str, Any] = {
            "experiment_id": exp_id,
            "replicate_count": 0,
            "replicate_names": list(),
            "experiment_url": EXP_URL.format(exp_id),
            "annotation_urls": list(),
            "replicate_urls": list(),
            "experiment_info": ""
        }
        with open(info_file, "w") as f:
            json.dump(info_dict, f, indent=4)


def make_experiment_list_file(experiment_ids: List[str], out_path: str):
    list_path: str = os.path.join(out_path, "experiment_list.txt")
    if not Path(list_path).exists():
        with open(list_path, "w") as f:
            f.write("\n".join(experiment_ids))
    else:
        with open(list_path, "a") as f:
            f.write("\n".join(experiment_ids))


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
        experiment_ids: List[str] = extract_experiment_ids_from_links(argument_dict["input"])

        make_experiment_list_file(experiment_ids, argument_dict["outdir"])

        id_count: int = len(experiment_ids)
        for i, exp_id in enumerate(experiment_ids):

            make_log_file(exp_id, argument_dict["outdir"])

            time = datetime.datetime.fromtimestamp(datetime.datetime.now().timestamp())
            timestamp = time.strftime('%H:%M:%S')

            print(str(i+1) + "/" + str(id_count) + ":", "Downloading annotation for", exp_id, "|", timestamp)
            download_annotation(exp_id, argument_dict["outdir"])
        print("All done!")
    else:
        print("Mode", argument_dict["mode"], "not recognised. Shutting down.")


if __name__ == "__main__":
    main()
