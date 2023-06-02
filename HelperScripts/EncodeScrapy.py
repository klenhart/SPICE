#!/bin/env python

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

import requests

from typing import Dict, Any, List

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse

URL: str = "https://www.encodeproject.org/batch_download/?type=Experiment&@id=/experiments/{0}/&files.processed=true"

# response = requests.get(url)

# if response.status_code == 200:
#     data = response.json()
#     file_links = data['@graph'][0]['files']
#     for file_link in file_links:
#         download_url = file_link['href']
#         print(download_url)
# else:
#     print(f"Error: {response.status_code} - {response.text}")


def extract_id(link: str) -> str:
    split_link: List[str] = link.split("/")
    if split_link[-1] == "":
        return split_link[-2]
    else:
        return split_link[-1]


def extract_experiment_ids_from_links(experiment_file_path: str) -> List[str]:
    with open(experiment_file_path, "r") as f:
        file = f.read()
    lines: List[str] = file.split("\n")
    return [extract_id(line) for line in lines]


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

    if argument_dict["mode"] == "alignment":
        pass
    elif argument_dict["mode"] == "annotation":
        experiment_ids: List[str] = extract_experiment_ids_from_links(argument_dict["input"])

    else:
        print("Mode", argument_dict["mode"], "not recognised. Shutting down.")


if __name__ == "__main__":
    main()
