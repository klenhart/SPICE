#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of main.
#
#  Gene_Assembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Gene_Assembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from Classes.SearchTree.SearchTree import SearchTree
from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.SequenceHandling.Gene import Gene
from Classes.SequenceHandling.Transcript import Transcript
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Exon import Exon

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
        self.gene_assembly: SearchTree = SearchTree(species + "_" + taxon_id)
        self.species: str = species
        self.taxon_id: str = taxon_id
        self.inclusion_filter_dict: Dict[str, List[str]] = dict()

    def update_inclusion_filter(self, key: str, possible_values: List[str]) -> None:
        self.inclusion_filter_dict.update({key: possible_values})

    def save(self) -> None:  # TODO Next thing to do: Save the results by flattening the SearchTree and then jsonDump
        pass

    def extract(self, gtf_path: str) -> None:
        """
        Method that extracts all genes, proteins, transcripts and exons from a gtf file and sorts them into
        a SearchTree class object.

        :param gtf_path: The absolute path to a gtf file.
        """

        # GTF Boy takes the path and can now be used to stream the gtf line by line.
        gtf_boy: GTFBoy = GTFBoy(gtf_path)

        # Iterate over the GTFBoys GTF file.
        for line in gtf_boy:
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
                    self.gene_assembly.insert_entry(gene)
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
                    # Find the corresponding gene in the SearchTree instance and insert the transcript
                    self.gene_assembly.find(transcript.get_id_gene()).add_entry("transcript", transcript)
                elif feature == "CDS":
                    # Make protein
                    protein: Protein = Protein()
                    # Use Proteins method to construct itself from a split gtf line.
                    protein.from_gtf_line(split_line)
                    # Externally set the taxon id.
                    protein.set_id_taxon(self.taxon_id)
                    # Find the corresponding gene in the SearchTree instance and insert the transcript
                    self.gene_assembly.find(protein.get_id_gene()).add_entry("transcript", protein)
                # elif feature == "exon": # TODO Integrate the use of Exons. For now it does not work since I intended
                #                               exons to only be included in proteins, but some exons will not be part
                #                               of proteins or be part of several
                #    # Make exon
                #    exon: Exon = Exon()
                #    # Use Exons method to construct itself from a split gtf line.
                #    exon.from_gtf_line(split_line)
                #    self.gene_assembly.find(exon.get_id_gene()).add_entry("exon", exon)


def main() -> None:
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens", "9606")
    gene_assembler.update_inclusion_filter("gene_biotype", ["protein_coding"])
    gene_assembler.extract("C:/Users/chris/Desktop/git/root/Homo_sapiens.GRCh38.107.gtf")


if __name__ == "__main__":
    main()
