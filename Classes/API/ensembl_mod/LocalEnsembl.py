#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
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
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import gzip
import os
import shutil
from contextlib import closing
from urllib import request

from Classes.API.ensembl_mod.EnsemblUtils import ping_ensembl, get_current_release, get_id_taxon, get_species_info


class LocalEnsembl:
    ftp_template: str = "http://ftp.ensembl.org/pub/release-{0}/gtf/{1}/{2}.{3}.{4}.gtf.gz"
    ftp_pep_template: str = "http://ftp.ensembl.org/pub/release-{0}/fasta/{1}/pep/{2}.{3}.pep.all.fa.gz"
    # 0=release_num
    # 1=species_name
    # 2=url_species_name
    # 3=assembly_default_name
    # 4=release_num

    filename_template: str = "{0}.{1}.{2}.gtf.gz" # url_species_name, assembly_default_name, release_num
    filename_pep_template: str = "{0}.{1}.pep.all.fa.gz" # url_species_name, assembly_default_name

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
        self.local_pep_zipname: str = self.filename_pep_template.format(self.url_species_name,
                                                                        self.assembly_default_species_name)

        self.local_filename: str = self.local_zipname[:-3]
        self.local_pep_filename: str = self.local_pep_zipname[:-3]

        self.ftp_address: str = self.ftp_template.format(self.release_num,
                                                         self.species_name,
                                                         self.url_species_name,
                                                         self.assembly_default_species_name,
                                                         self.release_num)

        self.ftp_pep_address: str = self.ftp_pep_template.format(self.release_num,
                                                             self.species_name,
                                                             self.url_species_name,
                                                             self.assembly_default_species_name)

    def get_species_name(self) -> str:
        return self.species_name

    def get_taxon_id(self) -> str:
        return self.taxon_id

    def download(self) -> str:
        if not self.is_downloaded():
            print("\tDownloading " + self.ftp_address + ".")
            with closing(request.urlopen(self.ftp_address)) as r:
                with open(os.path.join(self.goal_directory, self.local_zipname), 'wb') as f:
                    shutil.copyfileobj(r, f)

            # Unpacking file
            with gzip.open(os.path.join(self.goal_directory, self.local_zipname), 'rb') as f_in:
                with open(os.path.join(self.goal_directory, self.local_filename), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(os.path.join(self.goal_directory, self.local_zipname))
        else:
            print("Ensembl file already downloaded. (" + os.path.join(self.goal_directory, self.local_filename) + ")")
        return os.path.join(self.goal_directory, self.local_filename)

    def download_pep(self) -> str:
        if not self.is_pep_downloaded():
            print("\tDownloading " + self.ftp_pep_address + ".")
            with closing(request.urlopen(self.ftp_pep_address)) as r:
                with open(os.path.join(self.goal_directory, self.local_pep_zipname), 'wb') as f:
                    shutil.copyfileobj(r, f)

            # Unpacking file
            with gzip.open(os.path.join(self.goal_directory, self.local_pep_zipname), 'rb') as f_in:
                with open(os.path.join(self.goal_directory, self.local_pep_filename), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(os.path.join(self.goal_directory, self.local_pep_zipname))
        else:
            print("Coding sequences already downloaded. (" + os.path.join(self.goal_directory,
                                                                          self.local_pep_filename) + ")")
        return os.path.join(self.goal_directory, self.local_pep_filename)

    def get_release_num(self) -> str:
        return self.release_num

    def remove(self) -> None:
        if self.is_downloaded():
            os.remove(os.path.join(self.goal_directory, self.local_filename))

    def remove_pep(self) -> None:
        if self.is_pep_downloaded():
            os.remove(os.path.join(self.goal_directory, self.local_pep_filename))

    @property
    def ping(self) -> bool:
        return ping_ensembl()

    def is_downloaded(self) -> bool:
        if os.path.isfile(os.path.join(self.goal_directory, self.local_filename)):
            return True
        else:
            return False

    def is_pep_downloaded(self) -> bool:
        if os.path.isfile(os.path.join(self.goal_directory, self.local_pep_filename)):
            return True
        else:
            return False


def main():
    pass


if __name__ == "__main__":
    main()
