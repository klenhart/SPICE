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
from Classes.API.ensembl_mod.RemoteEnsembl import RemoteEnsembl
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo
from Classes.SequenceHandling.Protein import Protein
from Classes.TreeGrow.TreeGrow import TreeGrow
from Classes.WriteGuard.WriteGuard import WriteGuard

from typing import Dict, Any, List
from tqdm import tqdm
from datetime import date

import os.path
import shutil
import json


def check_library_status(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]) -> None:
    library_info["info"]["gene_count"] = gene_assembler.get_gene_count()
    library_info["info"]["transcript_count"] = gene_assembler.get_transcript_count()
    library_info["info"]["protein_count"] = gene_assembler.get_transcript_count()
    library_info["info"]["collected_sequences_count"] = gene_assembler.get_collected_sequences_count()
    library_info["info"]["fas_scored_sequences_count"] = gene_assembler.get_fas_scored_count()

    # Check sequence collection
    if gene_assembler.get_protein_count(True) > 0:
        if library_info["status"]["02_sequence_collection"]:
            print("Proteins without sequence found.")
            print("Will attempt to acquire sequence.")
        library_info["status"]["02_sequence_collection"] = False
    else:
        library_info["status"]["02_sequence_collection"] = True

    # Check small protein removal
    flag: bool = all([len(transcript) > 10 for transcript in gene_assembler.get_transcripts() if isinstance(transcript,
                                                                                                            Protein)])
    if not flag:
        if library_info["status"]["03_small_protein_removing"]:
            print("Protein sequences below length threshold found.")
            print("Will remove them.")
        library_info["status"]["03_small_protein_removing"] = False
    else:
        library_info["status"]["03_small_protein_removing"] = True

    # Check implicit FAS scoring
    fas_dict_list: List[Dict[str, Dict[str, float]]] = [gene.fas_dict for gene in gene_assembler.get_genes(False, True)]
    flag = True
    for fas_dict in fas_dict_list:
        for key in fas_dict.keys():
            if fas_dict[key][key] == -1.0:
                flag = False
    if not flag:
        if library_info["status"]["04_implicit_fas_scoring"]:
            print("Not yet computed implicit FAS scores found.")
            print("Will compute them.")
        library_info["status"]["04_implicit_fas_scoring"] = False
    else:
        library_info["status"]["04_implicit_fas_scoring"] = True

    # Check fasta generation
    with open(path_dict["transcript_fasta"], "r") as f:
        fasta_length = len(f.read().split("\n"))
    gene_list: List[Gene] = gene_assembler.get_genes(False, True)
    output_list: List[str] = list()
    for gene in gene_list:
        output_list.append(gene.fasta)
    new_fasta_length = len("\n".join(output_list).split("\n"))
    if fasta_length != new_fasta_length:
        if library_info["status"]["05_fasta_generation"]:
            print("Fasta file differs in size from what was expected.")
            print("Will regenerate it.")
        library_info["status"]["05_fasta_generation"] = False
    else:
        library_info["status"]["05_fasta_generation"] = True

    # Check pairing generation
    pairings_dict: Dict[str, str] = dict()
    for gene in gene_list:
        pairings_dict[gene.get_id()] = gene.make_pairings()
    with open(path_dict["transcript_pairings"], "r") as f:
        old_pairing_length = len(str(json.load(f)))
    new_pairing_length = len(str(pairings_dict))
    if old_pairing_length != new_pairing_length:
        if library_info["status"]["06_pairing_generation"]:
            print("Pairing file differs in size from what was expected.")
            print("Will regenerate it.")
        library_info["status"]["06_pairing_generation"] = False
    else:
        library_info["status"]["06_pairing_generation"] = True

    # Check ids tsv generation
    with open(path_dict["transcript_ids"], "r") as f:
        old_length = len(f.read())
    id_list: List[str] = list()
    for gene in gene_list:
        protein_list: List[Protein] = gene.get_proteins(False, True)
        for protein in protein_list:
            id_list.append(protein.make_header() + "\tncbi" + str(protein.get_id_taxon()))
    new_length = len(id_list)
    if old_length != new_length:
        if library_info["status"]["07_id_tsv_generation"]:
            print("ID tsv differs in size from what was expected.")
            print("Will regenerate it.")
        library_info["status"]["07_id_tsv_generation"] = False
    else:
        library_info["status"]["07_id_tsv_generation"] = True

    # Check FAS calculation
    if library_info["info"]["fas_scored_sequences_count"] < library_info["info"]["protein_count"]:
        library_info["status"]["09_fas_scoring"] = False
    else:
        library_info["status"]["09_fas_scoring"] = True
    library_info.save()


def collect_sequences(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]) -> None:
    incomplete_gene_list: List[Gene] = gene_assembler.get_genes(True)
    save_marker: int = 0
    for gene in tqdm(incomplete_gene_list, ncols=100,
                     total=len(incomplete_gene_list), desc="Sequence collection progress:"):

        save_marker += 1
        incomplete_proteins_list: List[Protein] = gene.get_proteins(True)
        results: List[Dict[str, str]] = RemoteEnsembl.collect_sequences(incomplete_proteins_list)
        for result in results:
            if "error" in result.keys():
                for protein in incomplete_proteins_list:
                    print("\n\t", protein.get_id(), " is depreciated. Removing from library.")
                    gene.delete_transcript(protein.get_id())
                    library_info["info"]["transcript_count"] = gene_assembler.get_transcript_count()
                    library_info["info"]["protein_count"] = gene_assembler.get_protein_count()

                    library_info.save()
                break
            else:
                gene.set_sequence_of_transcript(result["query"], result["seq"])
        if save_marker % 250 == 0:
            gene_assembler.save(path_dict["transcript_json"])
            library_info["info"]["collected_sequences_count"] = gene_assembler.get_collected_sequences_count()
            library_info.save()
    gene_assembler.clear_empty_genes()
    gene_assembler.save(path_dict["transcript_json"])
    info: Dict[str, Any] = library_info["info"]
    library_info["info"]["collected_sequences_count"] = gene_assembler.get_collected_sequences_count()
    sequence_collection_flag: bool = info["protein_count"] == info["collected_sequences_count"]
    library_info["status"]["02_sequence_collection"] = sequence_collection_flag
    library_info["info"]["gene_count"] = gene_assembler.get_gene_count()
    library_info.save()


def remove_small_proteins(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]):
    gene_list: List[Gene] = gene_assembler.get_genes()
    for gene in tqdm(gene_list, ncols=100, total=len(gene_list), desc="Small protein removal progress"):
        protein_list: List[Protein] = gene.get_proteins()
        for protein in protein_list:
            if len(protein) < 11:
                gene.delete_transcript(protein.get_id())
    gene_assembler.clear_empty_genes()
    gene_assembler.save(path_dict["transcript_json"])
    library_info["info"]["gene_count"] = gene_assembler.get_gene_count()
    library_info["info"]["transcript_count"] = gene_assembler.get_transcript_count()
    library_info["info"]["protein_count"] = gene_assembler.get_protein_count()
    library_info["info"]["collected_sequences_count"] = gene_assembler.get_collected_sequences_count()
    library_info["info"]["fas_scored_sequences_count"] = gene_assembler.get_fas_scored_count()
    library_info["status"]["03_small_protein_removing"] = True
    library_info.save()


def calculate_implicit_fas_scores(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]):
    gene_list: List[Gene] = gene_assembler.get_genes()
    for gene in tqdm(gene_list, ncols=100, total=len(gene_list), desc="Implicit FAS score collection progress"):
        gene.calculate_implicit_fas_scores()
    gene_assembler.save(path_dict["transcript_json"])
    library_info["status"]["04_implicit_fas_scoring"] = True
    library_info.save()


def generate_fasta_file(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]):
    gene_list: List[Gene] = gene_assembler.get_genes(False, True)
    output_list: List[str] = list()
    for gene in tqdm(gene_list, ncols=100, total=len(gene_list), desc="Fasta generation process"):
        output_list.append(gene.fasta)
    with open(path_dict["transcript_fasta"], "w") as f:
        f.write("\n".join(output_list))
    library_info["status"]["05_fasta_generation"] = True
    library_info.save()


def generate_pairings(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]):
    gene_list: List[Gene] = gene_assembler.get_genes(False, True)
    pairings_dict: Dict[str, str] = dict()
    for gene in tqdm(gene_list, ncols=100, total=len(gene_list), desc="Pairing generation process"):
        pairings_dict[gene.get_id()] = gene.make_pairings()

    with open(path_dict["transcript_pairings"], "w") as f:
        json.dump(pairings_dict, f, indent=4)

    library_info["status"]["06_pairing_generation"] = True
    library_info.save()


def generate_ids_tsv(gene_assembler: GeneAssembler, library_info: LibraryInfo, path_dict: Dict[str, str]):
    gene_list: List[Gene] = gene_assembler.get_genes(False, True)
    output_list: List[str] = list()
    for gene in tqdm(gene_list, ncols=100, total=len(gene_list), desc="Generating phyloprofile ids"):
        protein_list: List[Protein] = gene.get_proteins(False, True)
        for protein in protein_list:
            output_list.append(protein.make_header() + "\tncbi" + str(protein.get_id_taxon()))

    with open(path_dict["transcript_ids"], "w") as f:
        f.write("\n".join(output_list))

    library_info["status"]["07_id_tsv_generation"] = True
    library_info.save()


def main():
    ####################################################################
    # SETUP ARGS PARSER

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
    # LocalEnsembl SETUP.

    # Acquire the local ensembl file.
    local_ensembl: LocalEnsembl = LocalEnsembl(argument_dict["species"],
                                               argument_dict["outdir"],
                                               argument_dict["release"])

    library_name: str = "fade_lib_" + local_ensembl.get_species_name() + "_" + argument_dict["release"]

    ####################################################################
    # DATASTRUCTURE AND PATH SETUP

    gene_assembler: GeneAssembler = GeneAssembler(local_ensembl.get_species_name(),
                                                  local_ensembl.get_taxon_id())

    print("#00 Building library.")
    # Check if already exists:
    if os.path.exists(os.path.join(argument_dict["outdir"], library_name)) and not argument_dict["force"]:
        print("\tLibrary \"" + os.path.join(argument_dict["outdir"], library_name) + "\" already exists.")
        print("\tLoading existing library.")

        with open(os.path.join(argument_dict["outdir"], library_name, "paths.json"), "r") as f:
            path_dict: Dict[str, str] = json.load(f)
        gene_assembler.load(path_dict["transcript_json"])

        print("#01 Collecting transcripts information.")
        print("\tTranscript information already collected.")

        library_info: LibraryInfo = LibraryInfo(path_dict["info"])
        check_library_status(gene_assembler, library_info, path_dict)

    else:
        if os.path.exists(os.path.join(argument_dict["outdir"], library_name)) and argument_dict["force"]:
            shutil.rmtree(os.path.join(argument_dict["outdir"], library_name))

        print("#01 Collecting transcripts information.")

        ####################################################################
        # DIRECTORY SYSTEM BUILDING

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
                                                                     "transcript_set.json"),
                                     "transcript_fasta": os.path.join(argument_dict["outdir"],
                                                                      library_name,
                                                                      "transcript_data",
                                                                      "transcript_set.fasta"),
                                     "transcript_pairings": os.path.join(argument_dict["outdir"],
                                                                         library_name,
                                                                         "transcript_data",
                                                                         "transcript_pairings.json"),
                                     "transcript_ids": os.path.join(argument_dict["outdir"],
                                                                    library_name,
                                                                    "transcript_data",
                                                                    "phyloprofile_ids.tsv")
                                     }

        tree_grow: TreeGrow = TreeGrow(path_dict)
        tree_grow.create_folders()
        tree_grow.put_path_json()

        ####################################################################
        # LOCAL ENSEMBL DOWNLOAD

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

        print("\tSaving info.yaml at " + path_dict["info"])

        library_info: LibraryInfo = LibraryInfo(path_dict["info"])

        # Save base info
        library_info["fade_version"] = "0.1"
        library_info["init_date"] = str(date.today())
        library_info["last_edit"] = str(date.today())
        library_info["commandline_args"] = argument_dict
        library_info["info"] = {"species": local_ensembl.get_species_name(),
                                "taxon_id": local_ensembl.get_taxon_id(),
                                "release": local_ensembl.get_release_num(),
                                "gene_count": gene_assembler.get_gene_count(),
                                "transcript_count": gene_assembler.get_transcript_count(),
                                "protein_count": gene_assembler.get_protein_count(),
                                "collected_sequences_count": gene_assembler.get_collected_sequences_count(),
                                "fas_scored_sequences_count": gene_assembler.get_fas_scored_count()
                                }
        library_info["status"] = {"01_id_collection": True,
                                  "02_sequence_collection": False,
                                  "03_small_protein_removing": False,
                                  "04_implicit_fas_scoring": False,
                                  "05_fasta_generation": False,
                                  "06_pairing_generation": False,
                                  "07_id_tsv_generation": False,
                                  "08_sequence_annotation": False,
                                  "09_fas_scoring": False
                                  }

    ####################################################################
    # COLLECT SEQUENCES

    print("#02 Collecting sequences.")
    if not library_info["status"]["02_sequence_collection"]:
        # Collect sequences.
        with WriteGuard(path_dict["transcript_json"], path_dict["transcript_data"]):
            # Collect the sequences for each incomplete gene.
            print("\tSequence Collection Run 0/2")
            collect_sequences(gene_assembler, library_info, path_dict)
            print("\tSequence Collection Run 1/2")
            collect_sequences(gene_assembler, library_info, path_dict)
            print("\tSequence Collection Run 2/2")
    else:
        print("\tSequences already collected.")

    ####################################################################
    # CLEAR SMALL PROTEINS

    print("#03 Removing small proteins.")
    if not library_info["status"]["03_small_protein_removing"]:
        remove_small_proteins(gene_assembler, library_info, path_dict)
    else:
        print("\tSmall proteins already removed.")

    ####################################################################
    # IMPLICIT FAS SCORING

    print("#04 Calculating implicit FAS scores.")
    if not library_info["status"]["04_implicit_fas_scoring"]:
        calculate_implicit_fas_scores(gene_assembler, library_info, path_dict)
    else:
        print("\tImplicit FAS scores already calculated.")

    ####################################################################
    # CREATE FASTA FILE

    print("#05 Generating FASTA file for all sequences missing FAS scores.")
    if not library_info["status"]["05_fasta_generation"]:
        generate_fasta_file(gene_assembler, library_info, path_dict)
    else:
        print("\tFasta file already generated.")

    ####################################################################
    # CREATE PAIRINGS FOR ALL GENES

    print("#06 Creating protein pairings for all genes.")
    if not library_info["status"]["06_pairing_generation"]:
        generate_pairings(gene_assembler, library_info, path_dict)
    else:
        print("\tPairings already generated.")

    ####################################################################
    # CREATE OF ALL IDS.

    print("#07 Generating phyloprofile IDs for all proteins missing FAS scores.")
    if not library_info["status"]["07_id_tsv_generation"]:
        generate_ids_tsv(gene_assembler, library_info, path_dict)
    else:
        print("\tIDs already generated.")

    ####################################################################

    # TODO Calculate annotation (requires FAS location)
    # TODO Filter proteins with less than 10 aminoacids.


if __name__ == "__main__":
    main()
