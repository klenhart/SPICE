#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  RemoteEnsembl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RemoteEnsembl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


from typing import Type, List, Dict, Any
from time import sleep

import requests
import sys

from Classes.API.ensembl_mod.EnsemblUtils import chunks, make_request_data
from Classes.SequenceHandling.Protein import Protein


class RemoteEnsembl:

    @staticmethod
    def collect_sequences(proteins: List[Protein]) -> List[Dict[str, str]]:
        error_flag: bool = False

        request_chunks: List[List[str]] = list(chunks([protein.get_id() for protein in proteins], 50))
        ensembl_requests: List[str] = list()

        for chunk in request_chunks:
            ensembl_requests.append(make_request_data(chunk))

        server: str = "https://rest.ensembl.org/sequence/id"
        headers: Dict[str, str] = {"Content-Type": "application/json", "Accept": "application/json"}

        output_list: List[Dict[str, str]] = list()

        for i, request in enumerate(ensembl_requests):
            for x in range(4):
                r = requests.post(server, headers=headers, data=request)
                decoded: List[Dict[str, str]] = r.json()
                if r.ok:
                    break
                elif x >= 3:
                    decoded: Dict[str, str] = r.json()
                    if "error" in decoded.keys():
                        if decoded["error"] == "No results found":
                            error_flag = True
                            decoded = [decoded]
                    else:
                        r.raise_for_status()
                        sys.exit()
                sleep(0.1)
            for result in decoded:
                if "error" not in result.keys():
                    output_list.append(result)
            if len(output_list) == 0 and error_flag:
                output_list.append({"error": "No results found"})
        return output_list


def main():
    pass


if __name__ == "__main__":
    main()
