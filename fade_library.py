#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  fade_library is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  fade_library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.API.ensembl_mod.LocalEnsembl import LocalEnsembl
from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.TreeGrow.TreeGrow import TreeGrow

from typing import Dict, Any
import os.path


def main():
    # Set up the args parser.
    argument_parser: ReduxArgParse = ReduxArgParse(["--outdir", "--species", "--release"],
                                                   [str, str, str],
                                                   ["store", "store", "store"],
                                                   [None, None, None],
                                                   ["Directory the library will be generated in.",
                                                    "Species of the library.",
                                                    "Ensembl release of the library."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()

    ####################################################################

    # Acquire the local ensembl file.
    local_ensembl: LocalEnsembl = LocalEnsembl(argument_dict["species"],
                                               argument_dict["outdir"],
                                               argument_dict["release"])

    ####################################################################

    # Build the library directory system
    library_name: str = "fade_lib_" + local_ensembl.get_species_name() + "_" + argument_dict["release"]
    path_dict: Dict[str, str] = {"root": os.path.join(argument_dict["outdir"],
                                                      library_name),
                                 "info": os.path.join(argument_dict["outdir"],
                                                      library_name,
                                                      "info.tsv"),
                                 "fas_data": os.path.join(argument_dict["outdir"],
                                                          library_name,
                                                          "fas_data"),
                                 "fas_temp": os.path.join(argument_dict["outdir"],
                                                          library_name,
                                                          "fas_data",
                                                          "temp"),
                                 "transcript_data": os.path.join(argument_dict["outdir"],
                                                                 library_name,
                                                                 "transcript_data"),
                                 "transcript_json": os.path.join(argument_dict["outdir"],
                                                                 library_name,
                                                                 "transcript_data",
                                                                 "transcript_set.json")
                                 }

    tree_grow: TreeGrow = TreeGrow(path_dict)
    tree_grow.create_folders()
    tree_grow.put_path_json()

    ####################################################################

    # Download the local ensembl file.
    gtf_path: str = local_ensembl.download()

    # Extract the file
    gene_assembler: GeneAssembler = GeneAssembler(local_ensembl.get_species_name(),
                                                  local_ensembl.get_taxon_id())
    gene_assembler.update_inclusion_filter("gene_biotype", ["protein_coding"])
    gene_assembler.update_inclusion_filter("transcript_biotype", ["protein_coding", "nonsense_mediated_decay"])
    gene_assembler.extract(gtf_path)
    gene_assembler.save(path_dict["transcript_json"])

    # Delete the file after successful extraction.
    local_ensembl.remove()

    ####################################################################


if __name__ == "__main__":
    main()
