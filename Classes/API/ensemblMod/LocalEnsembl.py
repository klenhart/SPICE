#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  LocalEnsembl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  LocalEnsembl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import gzip
import os
import shutil
from contextlib import closing
from typing import Type
from typing import List
from urllib import request

from Classes.API.ensemblMod.ensemblUtils import ping_ensembl, get_current_release, get_id_taxon, get_species_info


class LocalEnsembl:
    ftp_template: str = "http://ftp.ensembl.org/pub/release-{0}/gtf/{1}/{2}.{3}.{4}.gtf.gz"
    filename_template: str = "{0}.{1}.{2}.gtf.gz"
    # 0=release_num
    # 1=species_name
    # 2=url_species_name
    # 3=assembly_default_species_name
    # 4=release_num

    def __init__(self, raw_species: str, goal_directory: str, release_num: str = "default") -> None:
        self.goal_directory: str = goal_directory

        self.raw_species_name: str = raw_species

        if release_num == "default":
            self.release_num: str = get_current_release()
        else:
            self.release_num: str = release_num

        self.species_name: str
        self.assembly_default_species_name: str
        self.url_species_name: str
        self.species_name, self.url_species_name, self.assembly_default_species_name = get_species_info(raw_species)
        self.taxon_id: str = get_id_taxon(self.species_name)

        self.local_zipname: str = self.filename_template.format(self.url_species_name,
                                                                self.assembly_default_species_name,
                                                                self.release_num)

        self.local_filename: str = self.local_zipname[:-3]

        self.ftp_address: str = self.ftp_template.format(self.release_num,
                                                         self.species_name,
                                                         self.url_species_name,
                                                         self.assembly_default_species_name,
                                                         self.release_num)

    def download(self) -> str:
        with closing(request.urlopen(self.ftp_address)) as r:
            with open(self.goal_directory + self.local_zipname, 'wb') as f:
                shutil.copyfileobj(r, f)

        # Unpacking file
        with gzip.open(self.goal_directory + self.local_zipname, 'rb') as f_in:
            with open(self.goal_directory + self.local_filename, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(self.goal_directory + self.local_zipname)
        return self.goal_directory + self.local_filename

    @property
    def ping(self) -> bool:
        return ping_ensembl()


def main():
    local_ensembl1 = LocalEnsembl("human", "C:/Users/chris/Desktop/git/root/", "107")
    print(local_ensembl1.ping)
    print(local_ensembl1.release_num)
    print(local_ensembl1.taxon_id)
    print(local_ensembl1.download())

    local_ensembl2 = LocalEnsembl("human", "C:/Users/chris/Desktop/git/root/")
    print(local_ensembl2.ping)
    print(local_ensembl2.release_num)
    print(local_ensembl2.taxon_id)


if __name__ == "__main__":
    main()
