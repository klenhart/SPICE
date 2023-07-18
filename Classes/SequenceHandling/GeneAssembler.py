#!/bin/env python


#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  GeneAssembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GeneAssembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.PassPath.PassPath import PassPath
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.Transcript import Transcript
from Classes.SequenceHandling.Protein import Protein
import os
from tqdm import tqdm
import json
from typing import List, Dict, Any, Set, Iterator, Tuple


class GeneAssembler:
    """
    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
     Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl
      identifier such as a scaffold ID, without any additional content such as species or assembly.
       See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position* of the feature, with sequence numbering starting at 1.
    end - End position* of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1'
     that the second base is the first base of a codon, and so on.
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

    Source: https://www.ensembl.org/info/website/upload/gff.html
    """

    GTF_MASK: List[str] = ["seqname", "source", "feature",
                           "start", "end", "score",
                           "strand", "frame", "attribute"]

    def __init__(self, species: str, taxon_id: str):
        self.gene_assembly: Dict[str, Gene] = dict()
        self.species: str = species
        self.taxon_id: str = taxon_id
        self.inclusion_filter_dict: Dict[str, List[str]] = dict()

    def __getitem__(self, gene_id) -> Gene:
        return self.gene_assembly[gene_id]

    def __contains__(self, gene_id) -> bool:
        return gene_id in self.gene_assembly.keys()

    def update_inclusion_filter(self, key: str, possible_values: List[str]) -> None:
        self.inclusion_filter_dict.update({key: possible_values})

    def save_fas(self, pass_path: PassPath) -> None:
        index_dict: Dict[str, str] = dict()
        for index, index_sub_dict, entry_dict in GeneAssembler.fas_to_dict_iter(self.gene_assembly):
            index_dict.update(index_sub_dict)
            with open(os.path.join(pass_path["fas_scores"], index), "w") as f:
                json.dump(entry_dict, f, indent=4)
        with open(pass_path["fas_index"], "w") as f:
            json.dump(index_dict, f, indent=4)

    def save_info(self, pass_path: PassPath) -> None:
        json_dict: Dict[str, Dict[str, Any]] = GeneAssembler.to_dict(self.gene_assembly, "info")
        with open(pass_path["transcript_info"], "w") as f:
            json.dump(json_dict, f, indent=4)

    def save_seq(self, pass_path: PassPath) -> None:
        json_dict: Dict[str, Dict[str, Any]] = GeneAssembler.to_dict(self.gene_assembly, "seq")
        with open(pass_path["transcript_seq"], "w") as f:
            json.dump(json_dict, f, indent=4)

    def load(self, pass_path: PassPath) -> None:
        with open(pass_path["transcript_info"], "r") as f:
            info_dict: Dict[str, Dict[str, Any]] = json.load(f)
        with open(pass_path["transcript_seq"], "r") as f:
            seq_dict: Dict[str, Dict[str, Any]] = json.load(f)
        fas_dict: Dict[str, Dict[str, Any]] = dict()
        with open(pass_path["fas_index"], "r") as f1:
            fas_index: Dict[str, str] = json.load(f1)
            for path in set(fas_index.values()):
                with open(os.path.join(pass_path["fas_scores"], path), "r") as f2:
                    fas_sub_dict: Dict[str, Dict[str, Any]] = json.load(f2)
                    fas_dict.update(fas_sub_dict)
        self.gene_assembly = GeneAssembler.from_dict(info_dict, seq_dict, fas_dict)

    def integrate_fas_json(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            distance_dict: Dict[str, Dict[str, Dict[str, float]]] = json.load(f)
        count: int = 0
        for key_gene_id in tqdm(distance_dict.keys(),
                                ncols=100,
                                total=len(list(distance_dict.keys())),
                                desc="Integrating FAS process"):
            if key_gene_id in self:
                gene: Gene = self.gene_assembly[key_gene_id]
                for key_prot_id1 in distance_dict[key_gene_id].keys():
                    if key_prot_id1 in gene.get_fas_dict().keys():
                        for key_prot_id2 in distance_dict[key_gene_id][key_prot_id1].keys():
                            if key_prot_id2 in gene.get_fas_dict()[key_prot_id1].keys():
                                count += 1
                                value: float = distance_dict[key_gene_id][key_prot_id1][key_prot_id2]
                                gene.get_fas_dict()[key_prot_id1][key_prot_id2] = value
        print("Integrated ", count, " FAS scores.")

    def extract_tags(self) -> List[str]:
        tag_set: Set[str] = set()
        for transcript in self.get_transcripts():
            for tag in transcript.get_tags():
                tag_set.add(tag)
        return list(tag_set)

    def extract(self, gtf_path: str) -> None:
        """
        Method that extracts all genes, proteins, transcripts and exons from a gtf file and sorts them into the
        GeneAssembler.

        :param gtf_path: The absolute path to a gtf file.
        """

        # GTF Boy takes the path and can now be used to stream the gtf line by line.
        gtf_boy: GTFBoy = GTFBoy(gtf_path)

        # Iterate over the GTFBoys GTF file.
        for line in tqdm(gtf_boy, ncols=100, total=gtf_boy.total_lines, desc="Extract GTF Progress"):
            # Skip the header.
            if line.startswith("#!"):
                continue
            else:
                # Split each line.
                split_line: List[str] = line.split("\t")
                # Check the feature CDS for protein. (Coding Sequence?)
                feature: str = split_line[2]
                if not GTFBoy.has_values(self.inclusion_filter_dict, split_line):
                    continue
                elif feature == "gene":
                    # Make gene
                    gene: Gene = Gene()
                    # Use Genes method to construct itself from a split gtf line.
                    gene.from_gtf_line(split_line)
                    # Externally set the taxon id and species.
                    gene.set_id_taxon(self.taxon_id)
                    gene.set_species(self.species)
                    # Insert the gene into the SearchTree instance
                    self.gene_assembly[gene.get_id()] = gene
                elif feature == "transcript":
                    # Do not extract transcripts that will have a protein counterpart in the GTF due to their biotype.
                    if GTFBoy.has_attribute_value("transcript_biotype", "protein_coding", split_line[8]):
                        continue
                    else:
                        # Make transcript
                        transcript: Transcript = Transcript()
                        # Use Transcripts method to construct itself from a split gtf line.
                        transcript.from_gtf_line(split_line)
                        # Externally set the taxon id.
                        transcript.set_id_taxon(int(self.taxon_id))
                    self.gene_assembly[transcript.get_id_gene()].add_transcript(transcript, True)
                elif feature == "CDS":
                    # Do not extract transcripts that will have a protein counterpart in the GTF due to their biotype.
                    if GTFBoy.has_attribute_value("transcript_biotype", "nonsense_mediated_decay", split_line[8]):
                        continue
                    else:
                        # Make protein
                        protein: Protein = Protein()
                        # Use Proteins method to construct itself from a split gtf line.
                        protein.from_gtf_line(split_line)
                        # Externally set the taxon id.
                        protein.set_id_taxon(int(self.taxon_id))
                        self.gene_assembly[protein.get_id_gene()].add_transcript(protein, True)

    def get_genes(self, no_sequence_flag: bool = False, no_fas_flag: bool = False) -> List[Gene]:
        output_list: List[Gene] = list()
        if no_sequence_flag:
            for key in list(self.gene_assembly.keys()):
                gene: Gene = self.gene_assembly[key]
                if not gene.is_sequence_complete():
                    output_list.append(gene)
        elif no_fas_flag:
            for key in list(self.gene_assembly.keys()):
                gene: Gene = self.gene_assembly[key]
                if len(gene.get_proteins(False, True)) > 0:
                    output_list.append(gene)
        else:
            for key in list(self.gene_assembly.keys()):
                gene: Gene = self.gene_assembly[key]
                output_list.append(gene)
        return output_list

    def get_transcripts(self) -> List[Transcript]:
        output_list: List[Transcript] = list()
        for key in self.gene_assembly.keys():
            gene: Gene = self.gene_assembly[key]
            output_list += gene.get_transcripts()
        return output_list

    def clear_empty_genes(self) -> None:
        gene_list: List[Gene] = self.get_genes()
        for gene in gene_list:
            if len(gene.get_proteins()) == 0:
                del self.gene_assembly[gene.get_id()]

    def get_gene_count(self) -> int:
        return len(self.gene_assembly.keys())

    def get_transcript_count(self, no_sequence_flag: bool = False) -> int:
        return sum([gene.get_transcript_count(no_sequence_flag) for gene in self.get_genes()])

    def get_protein_count(self,
                          no_sequence_flag: bool = False,
                          no_fas_flag: bool = False) -> int:
        return sum([gene.get_protein_count(no_sequence_flag,
                                           no_fas_flag) for gene in self.get_genes()])

    def get_collected_sequences_count(self) -> int:
        return self.get_protein_count() - self.get_protein_count(True)

    def get_fas_scored_count(self) -> int:
        return self.get_protein_count() - self.get_protein_count(False, True)

    def get_fas_dist_matrix(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        dist_matrix: Dict[str, Dict[str, Dict[str, float]]] = dict()
        for gene in self.get_genes():
            dist_matrix[gene.get_id()] = gene.get_fas_dict()
        return dist_matrix

    def reset_fas(self) -> None:
        for gene in self.get_genes():
            gene.reset_fas()

    @staticmethod
    def fas_to_dict_iter(gene_assembly: Dict[str, Gene]) -> Iterator[List[Tuple[str, Dict[str, str], Dict[str, Dict[str, Any]]]]]:
        json_dict: Dict[str, Dict[str, Any]] = dict()
        index_dict: Dict[str, str] = dict()
        count: int = 0
        reset_count: int = 0

        for key in gene_assembly.keys():
            json_dict[key] = gene_assembly[key].to_dict("fas")
            index_dict[key] = "0" * (9 - len(str(count))) + str(count) + ".json"
            reset_count += 1
            if reset_count == 100:
                yield "0" * (9 - len(str(count))) + str(count) + ".json", index_dict, json_dict
                reset_count = 0
                json_dict: Dict[str, Dict[str, Any]] = dict()
                index_dict: Dict[str, str] = dict()
                count += 1
        if len(json_dict) > 0:
            yield "0" * (9 - len(str(count))) + str(count) + ".json", index_dict, json_dict

    @staticmethod
    def to_dict(gene_assembly: Dict[str, Gene], mode: str) -> Dict[str, Dict[str, Any]]:
        json_dict: Dict[str, Dict[str, Any]] = dict()
        for key in gene_assembly.keys():
            json_dict[key] = gene_assembly[key].to_dict(mode)
        return json_dict

    @staticmethod
    def from_dict(info_dict: Dict[str, Dict[str, Any]],
                  seq_dict: Dict[str, Dict[str, Any]],
                  fas_dict: Dict[str, Dict[str, Any]]) -> Dict[str, Gene]:
        output_dict: Dict[str, Gene] = dict()
        for key in info_dict.keys():
            new_gene: Gene = Gene()
            new_gene.from_dict(info_dict[key], seq_dict[key], fas_dict[key])
            output_dict[new_gene.get_id()] = new_gene
        return output_dict


def main() -> None:
    pass


if __name__ == "__main__":
    main()
