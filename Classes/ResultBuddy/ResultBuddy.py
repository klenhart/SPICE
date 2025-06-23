#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  ResultBuddy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ResultBuddy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
import json

from typing import Dict, Any, List

from tqdm import tqdm

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.PassPath.PassPath import PassPath
from Classes.ResultBuddy.ComparisonHandling.ComparisonAssembler import ComparisonAssembler
from Classes.ResultBuddy.ExpressionHandling.ConditionAssembler import ConditionAssembler
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.ResultBuddy.EWFDHandling.EWFDAssembler import EWFDAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.LibraryInfo import LibraryInfo
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript
from Classes.TreeGrow.TreeGrow import TreeGrow
from Classes.WriteGuard.WriteGuard import WriteGuard


class ResultBuddy:
    """
    Interface between spice_result and ExpressionAssembler, ConditionAssembler and ComparisonAssembler.
    """

    def __init__(self, library_path: str, output_path: str, initial_flag: bool = False, suffix: str = ""):
        self.library_path: str = library_path
        with open(os.path.join(library_path, "paths.json"), "r") as f:
            self.library_pass_path: PassPath = PassPath(json.load(f))
        library_info: LibraryInfo = LibraryInfo(self.library_pass_path["info"])
        species: str = library_info["info"]["species"]
        release: str = library_info["info"]["release"]
        fas_hex: str = library_info["info"]["fas_mode"]
        if len(suffix) > 0:
            filename: str = "spice_result_" + species + "_" + release + "_" + fas_hex + "_" + suffix
        else:
            filename: str = "spice_result_" + species + "_" + release + "_" + fas_hex
        self.result_path: str = output_path + "/" + filename

        if initial_flag:
            self.result_info: Dict[str, Any] = dict()

            self.result_info["species"] = species
            self.result_info["taxon_id"] = library_info["info"]["taxon_id"]
            self.result_info["release"] = release
            self.result_info["library_integrity_flag"] = all(library_info["status"].values())

            self.result_info["expression_imports"]: Dict[str, Dict[str, Dict[str, str]]] = dict()
            self.result_info["expression_imports"]["conditions"]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["replicates"]: Dict[str, Dict[str, Dict[str, str]]] = dict()

            result_paths: Dict[str, Any] = dict()
            result_paths["library_path"] = self.library_path
            result_paths["root"] = self.result_path
            result_paths["result_info"] = "info.json"
            result_paths["expression"] = "expression"
            result_paths["expression_replicates"] = "expression/replicates"
            result_paths["expression_conditions"] = "expression/conditions"
            result_paths["ewfd"] = "ewfd"
            result_paths["ewfd_replicates"] = "ewfd/replicates"
            result_paths["ewfd_conditions"] = "ewfd/conditions"
            result_paths["comparison"] = "comparison"

            self.result_pass_path: PassPath = PassPath(result_paths)
            tree_grow: TreeGrow = TreeGrow(result_paths)
            tree_grow.create_folders()
            tree_grow.put_path_json()

            self.__save_info()

        else:
            self.result_info: Dict[str, Any] = self.__load_info()
            self.result_pass_path: PassPath = PassPath(self.__load_paths())

    def __load_info(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "info.json"), "r") as f:
            result_info: Dict[str, Any] = json.load(f)
        return result_info

    def __save_info(self) -> None:
        with open(os.path.join(self.result_path, "info.json"), "w") as f:
            json.dump(self.result_info, f, indent=4)

    def __load_paths(self) -> Dict[str, Any]:
        with open(os.path.join(self.result_path, "paths.json"), "r") as f:
            result_paths: Dict[str, Any] = json.load(f)
        return result_paths

    def generate_ewfd_file(self, name: str, condition_flag: bool = False):
        """
        Initiates calculation of expression-weighted-functional-disturbance (EWFD)
        vectors and saves results.

        Parameters:
            name: str: Name of the currenty sample.
            condition_flag: bool: Is True if in mode "condition", Default: False (mode "expression")

        Loads:
            - Expreriment info saved in info.json

        Behavior:
            - Initiates an EWFDAssembler object which takes information about FAS scores and relative
              expression of transcritpts to calculate EWFD vectors
            - After calculation, info.json is updated and results are saved in either ewfd/replicates
              or ewfd/conditions

        Use Case:
            spice_result.py; mode "expression" or "condition"
        """
        if condition_flag:
            ewfd_type: str = "conditions"
        else:
            ewfd_type: str = "replicates"
        self.result_info = self.__load_info()
        species: str = self.result_info["species"]
        taxon_id: int = self.result_info["taxon_id"]
        expr_path: str = self.result_info["expression_imports"][ewfd_type][name]["expression_path"]
        # Initialization of EWFDAssembler object -> calculation of ewfd vectors
        ewfd: EWFDAssembler = EWFDAssembler(species, taxon_id, self.library_pass_path, expr_path, True, condition_flag)
        # Save results
        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, name):
            # Update result info (info.json)
            self.result_info = self.__load_info()
            filename: str = "/ewfd_" + name + ".json"
            ewfd_path: str = self.result_pass_path["ewfd_" + ewfd_type] + filename
            self.result_info["expression_imports"][ewfd_type][name]["ewfd_path"] = ewfd_path
            self.__save_info()
            # Saves ewfd instance as JSON file in /ewfd/replicates or /ewfd/conditions
            ewfd.save(self.result_info["expression_imports"][ewfd_type][name]["ewfd_path"])

    def build_condition(self, condition_name: str, replicate_names: List[str]):
        """
        Creates a ConditionAssembler object from given replicates to build a
        condition-level datastructure.

        Loads:
            - Experiment info saved in `info.json`.
            - Previously generated ExpressionAssembler datastructure for each replicate,
              calculated and saved in mode "expression".

        Behavior:
            - Initiates a ConditionAssembler object
            - Loops over all replicates, loads their pre-calculated ExpressionAssembler data and
              populates the ConditionAssembler datastructure for each gene
                - Adds absolute and relative expression values
                - Calculates average relative expression and adds values
            - Saves results in expression/conditions/expression_[cond_name].json 
        
        Note:
            - This will be called by spice_result.py in "condition" mode
        """
        # Load experiment info.json
        self.result_info = self.__load_info()
        # Initiate ConditionAssembler instance
        condition: ConditionAssembler = ConditionAssembler(self.library_pass_path,
                                                           condition_name,
                                                           True)
        # Loop over replicates of condition
        for name in replicate_names:
            # Creates an ExpressionAssembler instance but with empty ExpressionAssembler.expression_assembly["data"]
            expression: ExpressionAssembler = ExpressionAssembler(self.library_pass_path,
                                                                  name)
            # Load expression_assembly of replicate (calculated in "expression" mode)
            expression.load(self.result_info["expression_imports"]["replicates"][name]["expression_path"])
            # Calculate averaged relative expression
            condition.insert_expression(expression)

        # Keep all genes and transcripts, doesn't matter if they have expression or not.
        # condition.cleanse_assembly()

        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, condition_name):
            self.result_info = self.__load_info()
            filename: str = "/expression_" + condition_name + ".json"
            condition_path: str = self.result_pass_path["expression_conditions"] + filename
            new_condition_dict: Dict[str, Any] = {"replicates": replicate_names,
                                                  "expression_path": condition_path,
                                                  "ewfd_path": ""}
            self.result_info["expression_imports"]["conditions"][condition_name]: Dict[str, Dict[str, Any]] = dict()
            self.result_info["expression_imports"]["conditions"][condition_name] = new_condition_dict
            self.__save_info()
            condition.save(self.result_info["expression_imports"]["conditions"][condition_name]["expression_path"])

    def compare(self, condition_pair: List[str]):
        condition_pair.sort()
        comparison: ComparisonAssembler = ComparisonAssembler(condition_pair[0],
                                                              condition_pair[1],
                                                              self.result_pass_path,
                                                              self.result_info)
        return comparison

    def import_expression_gtf(self,
                              expression_path: str,
                              expression_name: str,
                              normalization: str,
                              expression_threshold: float = 1.0) -> None:
        """
        Initiaites and populates the ExpressionAssembler datastructure for spice_result
        mode "expression".

        Parameters:
            expression_path: Path to expression gtf file.
            expression_name: Name of the sample.
            normalization: Used expression normalization (e.g. "FPKM", "TPM"). Requires that
                           expression gtf file contains normalized expression value in attributes
                           column and is labeled as such (e.g. "FPKM": "5")
            expression_threshold: Expression of transcripts below this threshold will be set to
                                  zero. DEPRECATED: Invalid. Use expression_threshold = 0 for now
                                  until changed.

        Behavior:
            - Loads the expression gtf file
            - Initiates ExpressionAssembler datastructure which contains genes/transcripts
              present in Spice library
            - Adds expression to datastructure
            - Calculates relative expression of in-library transcripts for each gene
            - Saves results /expression/replicates and updates info.json
        
        Use Case:
            spice_result.py; mode "expression"
        
        Note:
            Allows for transcript and gene ids in the expression gtf file to differ
            from the reference gtf file (used in library generation). However, this requires
            that synonyms have been previously saved in the GeneAssembler datastructure
            when generating spice library.
        """
        transcript_to_protein_dict: Dict[str, str] = self.transcript_to_protein_map()
        transcript_to_gene_dict: Dict[str, str] = self.transcript_to_gene_map()
        synonym_to_transcript_dict: Dict[str, str] = self.synonym_to_transcript_map()
        expression_gtf: GTFBoy = GTFBoy(expression_path)
        expression_assembler: ExpressionAssembler = ExpressionAssembler(self.library_pass_path,
                                                                        expression_name,
                                                                        expression_path,
                                                                        normalization,
                                                                        True,
                                                                        expression_threshold)
        for line in tqdm(expression_gtf,
                         ncols=100,
                         total=expression_gtf.total_lines,
                         desc=expression_name + " GTF extraction progress"):
            if line.startswith("#"):
                continue
            else:
                split_line: List[str] = line.split("\t")
                line_dict: Dict[str, str] = GTFBoy.build_dict(split_line)
                # transcript_flag = True if current line is feature transcript (or novel_transcript?)
                transcript_flag: bool = line_dict["feature"] in ["transcript", "novel_transcript"]
                if "transcript_id" in line_dict.keys() and transcript_flag:
                    line_dict["transcript_id"] = line_dict["transcript_id"].split(".")[0].split(":")[-1]
                    # Checks if transcript in line currently looked at is in spice library
                    # Using synonym_to_transcript_dict allows for matching of alternative transcript names
                    # if they differ between reference GTF and expression GTF (however, this should not be the case)
                    transcript_in_lib_flag: bool = line_dict["transcript_id"] in synonym_to_transcript_dict.keys()
                    if transcript_in_lib_flag:
                        # if alternative transcript keys are used in expression gtf we revert them to the 
                        # transcript ids used in the reference/spice library
                        line_dict["transcript_id"] = synonym_to_transcript_dict[line_dict["transcript_id"]]
                        # This should always be True if the transcript is in the spice library
                        gene_in_lib_flag: bool = line_dict["transcript_id"] in transcript_to_gene_dict.keys()
                    else:
                        gene_in_lib_flag: bool = False
                    # This implicitly checks if the transcript is PROTEIN CODING or NMD bio-typed (because only
                    # these genes/transcripts are in the spice library)
                    if gene_in_lib_flag:
                        line_dict["transcript_id"] = synonym_to_transcript_dict[line_dict["transcript_id"]]
                        line_dict["gene_id"] = transcript_to_gene_dict[line_dict["transcript_id"]]
                        # Also add protein_id to line_dict
                        line_dict["protein_id"] = transcript_to_protein_dict[line_dict["transcript_id"]]
                        expression_assembler.insert_expression_dict(line_dict)
        # Keep all genes and transcripts, doesn't matter if they have expression or not.
        # expression_assembler.cleanse_assembly()
        # Calculate relative expression and further populate expression_assembler datastructure
        expression_assembler.calc_relative_expression()
        filename: str = "/expression_" + expression_name + ".json"
        expression_json_path: str = self.result_pass_path["expression_replicates"] + filename
        # Update info.json
        with WriteGuard(os.path.join(self.result_path, "info.json"), self.result_path, expression_name):
            self.result_info = self.__load_info()
            new_expression_dict: Dict[str, str] = {"origin": expression_path,
                                                   "expression_path": expression_json_path,
                                                   "ewfd_path": ""}
            self.result_info["expression_imports"]["replicates"][expression_name]: Dict[str, Dict[str, str]] = dict()
            self.result_info["expression_imports"]["replicates"][expression_name] = new_expression_dict
            self.__save_info()
        # Save Results
        expression_assembler.save(expression_json_path)

    def transcript_to_biotype_map(self) -> Dict[str, str]:
        """
        Generates a mapping from transcript IDs to their corresponding biotypes.

        Returns:
            Dict[str, str]: A dictionary where:
                - Keys are transcript IDs (ENST).
                - Values are biotype labels (e.g., "protein_coding", "nonsense_mediated_decay").

        Behavior:
            - Loads the GeneAssembler from the stored transcript and FAS data.
            - Iterates over all transcripts (Protein and Transcript objects).
            - Assigns the biotype based on object type:
                - Protein → "protein_coding"
                - Transcript → "nonsense_mediated_decay" (assumed default)

        Note:
            This function assumes that any Transcript object not subclassed as Protein is NMD,
            which is consistent with the design of the library generation pipeline as
            inclusion_filter_dict is specified accordingly in spice_library.py. If in the future
            also other biotypes should be included, both GeneAssembler's inclusion_filter_dict and
            this function need to be adjusted acccordingly.
        """
        transcript_to_biotype_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_biotype_map[transcript.id_transcript] = "protein_coding"
            elif isinstance(transcript, Transcript):
                transcript_to_biotype_map[transcript.id_transcript] = "nonsense_mediated_decay"
        return transcript_to_biotype_map

    def synonym_to_transcript_map(self) -> Dict[str, str]:
        """
        Creates a mapping from synonym names (or aliases) to canonical transcript IDs.

        Returns:
            Dict[str, str]: A dictionary where:
                - Keys are synonyms or canonical transcript IDs.
                - Values are the canonical transcript IDs (ENST).

        Behavior:
            - Loads the GeneAssembler instance.
            - For each transcript (Protein or Transcript object), collects:
                - Its canonical transcript ID.
                - All synonyms associated with it.
            - Maps all of them to the canonical transcript ID.

        Use Case:
            Enables flexible transcript lookup using either standard Ensembl IDs or known synonyms.
            This is useful in datasets or interfaces that use alternative naming schemes.
        Example:
                {
                    "ENST00000312345": "ENST00000312345",
                    "transcript_alias_1": "ENST00000312345",
                    "transcript_alias_2": "ENST00000312345"
                }
        """
        synonym_to_transcript_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                synonym_to_transcript_map[transcript.get_id_transcript()] = transcript.get_id_transcript()
                for synonym in transcript.get_synonyms():
                    synonym_to_transcript_map[synonym] = transcript.get_id_transcript()
            elif isinstance(transcript, Transcript):
                synonym_to_transcript_map[transcript.get_id()] = transcript.get_id()
                for synonym in transcript.get_synonyms():
                    synonym_to_transcript_map[synonym] = transcript.get_id()
        return synonym_to_transcript_map

    def transcript_to_protein_map(self) -> Dict[str, str]:
        """
        Creates a mapping from transcript IDs to corresponding protein IDs, if available.

        Returns:
            Dict[str, str]: A dictionary where:
                - Keys are transcript IDs (ENST).
                - Values are matching protein IDs (ENSP) for protein-coding transcripts.
                If no protein exists (e.g., for NMD transcripts), the value is an empty string.

        Behavior:
            - Loads the GeneAssembler instance using species/taxon info and the stored library path.
            - Iterates over all transcripts (both Protein and Transcript objects).
            - Distinguishes between protein-coding and non-coding/NMD transcripts using `isinstance`.

        Use Case:
            Allows for quick lookup of transcript id -> protein id association.

        """
        transcript_to_protein_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_protein_map[transcript.get_id_transcript()] = transcript.get_id()
            elif isinstance(transcript, Transcript):
                transcript_to_protein_map[transcript.get_id()] = ""
        return transcript_to_protein_map

    def transcript_to_gene_map(self) -> Dict[str, str]:
        """
        Creates a mapping from transcript IDs to their corresponding gene IDs.

        Returns:
            Dict[str, str]: A dictionary where:
                - Keys are transcript IDs (ENST).
                - Values are gene IDs (ENSG) that the transcripts belong to.

        Behavior:
            - Loads the GeneAssembler instance using species/taxon info and the stored library path.
            - Iterates over all Transcript and Protein objects, extracting transcript → gene relationships.

        Use Case:
            Allows for quick lookup of transcript -> gene association.
        """
        transcript_to_gene_map: Dict[str, str] = dict()
        gene_assembler: GeneAssembler = GeneAssembler(self.result_info["species"], self.result_info["taxon_id"])
        gene_assembler.load(self.library_pass_path)
        for transcript in gene_assembler.get_transcripts():
            if isinstance(transcript, Protein):
                transcript_to_gene_map[transcript.get_id_transcript()] = transcript.get_id_gene()
            elif isinstance(transcript, Transcript):
                transcript_to_gene_map[transcript.get_id()] = transcript.get_id_gene()
        return transcript_to_gene_map


def main():
    pass


if __name__ == "__main__":
    main()
