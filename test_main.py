import json
from typing import List, Dict, Any

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():


    with open("C:/Users/chris/Desktop/full_test_files/spice_lib_homo_sapiens_94_147/paths.json", "r") as f:
        paths_dict = json.load(f)
    pass_path: PassPath = PassPath(paths_dict)
    gene_assembler: GeneAssembler = GeneAssembler("homo_sapiens",
                                                  "9606")

    gene_assembler.load(pass_path)

    for gene in gene_assembler.get_genes():
        delete_list = list()
        for entry in gene.get_fas_dict().keys():
            if not entry.startswith("ENS"):
                delete_list.append(entry)
        for entry in delete_list:
            del gene.get_fas_dict()[entry]
        for entry in gene.get_fas_dict().keys():
            for del_entry in delete_list:
                del gene.get_fas_dict()[entry][del_entry]

    gene_assembler.save_fas(pass_path)



def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
