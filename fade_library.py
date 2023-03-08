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
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.TreeGrow.TreeGrow import TreeGrow

from typing import Dict, Any
import os.path
import shutil
import yaml
import json
from datetime import date


def main():
    # Set up the args parser.
    argument_parser: ReduxArgParse = ReduxArgParse(["--outdir", "--species", "--release", "--force", "--keepgtf"],
                                                   [str, str, str, bool],
                                                   ["store", "store", "store", "store_true", "store_true"],
                                                   [None, None, None, "?", "?"],
                                                   ["Directory the library will be generated in.",
                                                    "Species of the library.",
                                                    "Ensembl release of the library.",
                                                    "If the specified library already exists, it will be overwritten.",
                                                    "Keeps the ensembl GTF on the system after library setup."])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()

    if argument_dict["force"] is None:
        argument_dict["force"] = False

    if argument_dict["keepgtf"] is None:
        argument_dict["keepgtf"] = False

    ####################################################################

    # Acquire the local ensembl file.
    local_ensembl: LocalEnsembl = LocalEnsembl(argument_dict["species"],
                                               argument_dict["outdir"],
                                               argument_dict["release"])

    library_name: str = "fade_lib_" + local_ensembl.get_species_name() + "_" + argument_dict["release"]

    ####################################################################

    gene_assembler: GeneAssembler = GeneAssembler(local_ensembl.get_species_name(),
                                                  local_ensembl.get_taxon_id())

    # Check if already exists:
    if os.path.exists(os.path.join(argument_dict["outdir"], library_name)) and not argument_dict["force"]:
        print("Library \"" + os.path.join(argument_dict["outdir"], library_name) + "\" already exists.")
        print("Loading existing library.")

        with open(os.path.join(argument_dict["outdir"], library_name, "paths.json") , "r") as f:
            path_dict: Dict[str, str] = json.load(f)
        gene_assembler.load(path_dict["transcript_json"])

        print("Skipping collection of transcript information.")
    else:
        if os.path.exists(os.path.join(argument_dict["outdir"], library_name)) and argument_dict["force"]:
            shutil.rmtree(os.path.join(argument_dict["outdir"], library_name))

        print("Collecting transcript information.")

        ####################################################################

        # Build the library directory system
        path_dict: Dict[str, str] = {"root": os.path.join(argument_dict["outdir"],
                                                          library_name),
                                     "info": os.path.join(argument_dict["outdir"],
                                                          library_name,
                                                          "info.yaml"),
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
        gene_assembler.update_inclusion_filter("gene_biotype", ["protein_coding"])
        gene_assembler.update_inclusion_filter("transcript_biotype", ["protein_coding", "nonsense_mediated_decay"])
        gene_assembler.extract(gtf_path)
        gene_assembler.clear_empty_genes()
        gene_assembler.save(path_dict["transcript_json"])

        if not argument_dict["keepgtf"]:
            # Delete the file after successful extraction.
            local_ensembl.remove()

    print("Saving info.yaml at " + path_dict["info"])

    # Save base info
    info_dict: Dict[str, Any] = dict()
    info_dict["fade_version"] = "0.1"
    info_dict["date"] = str(date.today())
    info_dict["commandline_args"] = argument_dict
    info_dict["library_info"] = {"species": local_ensembl.get_species_name(),
                                 "taxon_id": local_ensembl.get_taxon_id(),
                                 "release": local_ensembl.get_release_num(),
                                 "gene_count": gene_assembler.get_gene_count(),
                                 "transcript_count": gene_assembler.get_transcript_count(),
                                 "protein_count": gene_assembler.get_protein_count(),
                                 "collected_sequences_count": gene_assembler.get_collected_sequences_count()
                                 }
    with open(path_dict["info"], "w") as f:
        yaml.dump(info_dict, f)

    ####################################################################


if __name__ == "__main__":
    main()
