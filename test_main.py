import json
import os
from typing import Dict, List, Any

from tqdm import tqdm

from Classes.GTFBoy.GTFBoy import GTFBoy
from Classes.ResultBuddy.ExpressionHandling.ExpressionAssembler import ExpressionAssembler
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from Classes.SequenceHandling.Protein import Protein
from Classes.SequenceHandling.Transcript import Transcript
from Classes.WriteGuard.WriteGuard import WriteGuard


def transcript_to_protein_mapper() -> Dict[str, str]:
    transcript_to_protein_map: Dict[str, str] = dict()
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens", "9606")
    gene_assembler.load("C:/Users/chris/Desktop/git/fade_lib_homo_sapiens_107/transcript_data/transcript_set.json")
    for transcript in gene_assembler.get_transcripts():
        if isinstance(transcript, Protein):
            transcript_to_protein_map[transcript.get_id_transcript()] = transcript.get_id()
        elif isinstance(transcript, Transcript):
            transcript_to_protein_map[transcript.get_id()] = ""
    return transcript_to_protein_map


def import_expression_gtf(expression_path: str,
                          expression_name: str,
                          normalization: str,
                          expression_threshold: float,
                          result_path: str,
                          library_path: str) -> None:
    # transcript_to_protein_dict: Dict[str, str] = transcript_to_protein_mapper()
    expression_gtf: GTFBoy = GTFBoy(expression_path)
    # expression_assembler: ExpressionAssembler = ExpressionAssembler(os.path.join(library_path,
    #                                                                              "transcript_data",
    #                                                                              "transcript_set.json"),
    #                                                                 expression_name,
    #                                                                 expression_path,
    #                                                                 normalization,
    #                                                                 True,
    #                                                                 expression_threshold)
    count = 0
    for line in tqdm(expression_gtf,
                     ncols=100,
                     total=expression_gtf.total_lines,
                     desc=expression_name + " GTF extraction progress"):
        if line.startswith("#"):
            continue
        else:
            split_line: List[str] = line.split("\t")
            line_dict: Dict[str, str] = GTFBoy.build_dict(split_line)
            transcript_flag: bool = line_dict["feature"] == "transcript"
            check_list: List[str] = [key + ": " + value for key, value in line_dict.items() if ("ENST" in value or "gene_id" == key) and transcript_flag]
            if len(check_list) == 2:
                count += 1

            # if all([key in line_dict.keys() for key in ["reference_id", "ref_gene_id"]]) and transcript_flag:
            #     line_dict["gene_id"] = line_dict["ref_gene_id"].split(".")[0]
            #     line_dict["transcript_id"] = line_dict["reference_id"].split(".")[0]
            #     # This implicitly checks if the transcript is PROTEIN CODING or NMD bio-typed.
            #     if line_dict["transcript_id"] in transcript_to_protein_dict.keys():
            #         line_dict["protein_id"] = transcript_to_protein_dict[line_dict["transcript_id"]]
            #         expression_assembler.insert_expression_dict(line_dict)
    print(count)
    # expression_assembler.cleanse_assembly()
    # expression_assembler.calc_relative_expression()


def main():
    import_expression_gtf("C:/Users/chris/Desktop/gtfs/ENCFF023EXJ.gtf",
                          "EXJ",
                          "FPKM",
                          0.0,
                          "",
                          "")

    # with open("C:/Users/chris/Desktop/expression_WTC11_ENCFF023EXJ_ENSv107.json", "r") as f:
    #     old_dict: Dict[str, Any] = json.load(f)
#
    # with open("C:/Users/chris/Desktop/git/result/expression/replicates/EXJ.json", "r") as f:
    #     new_dict: Dict[str, Any] = json.load(f)

if __name__ == "__main__":
    main()
