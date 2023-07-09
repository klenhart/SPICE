import json
import os
from typing import List, Dict, Any

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():
    output = []
    condition_list = list(open_json("C:/Users/chris/Desktop/big_test/spice_result_homo_sapiens_94_1ee_tissue_novel/info.json")["expression_imports"]["conditions"].keys())
    for condition_1 in condition_list:
        for condition_2 in condition_list:
            if condition_1 != condition_2:
                if ";".join(sorted([condition_1, condition_2])) not in output:
                    output.append(";".join(sorted([condition_1, condition_2])))

    print(" ".join(output))






def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
