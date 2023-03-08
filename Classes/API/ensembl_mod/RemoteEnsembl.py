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
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
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
    def collect_sequences(self, proteins: List[Protein]) -> Dict[str, str]:
        request_chunks: List[List[str]] = list(chunks([protein.get_id() for protein in proteins], 50))
        ensembl_requests: List[str] = list()

        for chunk in request_chunks:
            ensembl_requests.append(make_request_data(chunk))

        server: str = "https://rest.ensembl.org/sequence/id"
        headers: Dict[str, str] = {"Content-Type": "application/json", "Accept": "application/json"}

        for i, request in enumerate(ensembl_requests):
            for x in range(10):
                r = requests.post(server, headers=headers, data=request)
                if r.ok:
                    break
                elif x >= 9:
                    r.raise_for_status()
                    sys.exit()
                sleep(0.05)
            decoded: Any = r.json()
            print(decoded)


def main():
    pass


if __name__ == "__main__":
    main()
