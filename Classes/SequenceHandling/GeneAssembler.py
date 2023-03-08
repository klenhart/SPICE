#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
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
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.Transcript import Transcript
from Classes.SequenceHandling.Protein import Protein

from tqdm import tqdm
import json
from typing import List, Dict, Any


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
     that the second base is the first base of a codon, and so on..
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

    def update_inclusion_filter(self, key: str, possible_values: List[str]) -> None:
        self.inclusion_filter_dict.update({key: possible_values})

    def save(self, output_path: str) -> None:
        json_dict: Dict[str, Dict[str, Any]] = GeneAssembler.to_dict(self.gene_assembly)
        with open(output_path, "w") as f:
            json.dump(json_dict, f, indent=4)

    def load(self, input_path: str) -> None:
        with open(input_path, "r") as f:
            self.gene_assembly = GeneAssembler.from_dict(json.load(f))

    def extract(self, gtf_path: str) -> None:
        """
        Method that extracts all genes, proteins, transcripts and exons from a gtf file and sorts them into
        a SearchTree class object.

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
                        transcript.set_id_taxon(self.taxon_id)
                    self.gene_assembly[transcript.get_id_gene()].add_transcript(transcript)
                elif feature == "CDS":
                    # Make protein
                    protein: Protein = Protein()
                    # Use Proteins method to construct itself from a split gtf line.
                    protein.from_gtf_line(split_line)
                    # Externally set the taxon id.
                    protein.set_id_taxon(self.taxon_id)
                    self.gene_assembly[protein.get_id_gene()].add_transcript(protein)

    def get_genes(self, incomplete_flag: bool = False) -> List[Gene]:
        output_list: List[Gene] = list()
        if incomplete_flag:
            for key in list(self.gene_assembly.keys()):
                gene: Gene = self.gene_assembly[key]
                if not gene.is_sequence_complete():
                    output_list.append(gene)
        else:
            for key in list(self.gene_assembly.keys()):
                gene: Gene = self.gene_assembly[key]
                output_list.append(gene)
        return output_list

    def clear_empty_genes(self) -> None:
        gene_list: List[Gene] = self.get_genes()
        for gene in gene_list:
            if len(gene.get_transcripts()) == 0:
                del self.gene_assembly[gene.get_id()]

    def get_gene_count(self) -> int:
        return len(self.gene_assembly.keys())

    def get_transcript_count(self, incomplete_flag: bool = False) -> int:
        return sum([gene.get_transcript_count(incomplete_flag) for gene in self.get_genes()])

    def get_protein_count(self, incomplete_flag: bool = False) -> int:
        return sum([gene.get_protein_count(incomplete_flag) for gene in self.get_genes()])

    def get_collected_sequences_count(self) -> int:
        return self.get_protein_count(True)

    @staticmethod
    def to_dict(gene_assembly: Dict[str, Gene]) -> Dict[str, Dict[str, Any]]:
        json_dict: Dict[str, Dict[str, Any]] = dict()
        for key in gene_assembly.keys():
            json_dict[key] = gene_assembly[key].to_dict()
        return json_dict

    @staticmethod
    def from_dict(json_dict: Dict[str, Dict[str, Any]]) -> Dict[str, Gene]:
        output_dict: Dict[str, Gene] = dict()
        for key in json_dict.keys():
            new_gene: Gene = Gene()
            new_gene.from_dict(json_dict[key])
            output_dict[new_gene.get_id()] = new_gene
        return output_dict


def main() -> None:
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens", "9606")
    gene_assembler.update_inclusion_filter("gene_biotype", ["protein_coding"])
    gene_assembler.extract("C:/Users/chris/Desktop/git/root/Homo_sapiens.GRCh38.107.gtf")
    gene_assembler.save("C:/Users/chris/Desktop/git/root/extract.json")


if __name__ == "__main__":
    main()
