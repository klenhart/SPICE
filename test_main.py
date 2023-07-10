import copy
import json
import os
from typing import List, Dict, Any

from Classes.AssemblyVisualizer.RMSDOptimizer import RMSDOptimizerAlt
from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():
    with open(os.path.join("C:/Users/chris/Desktop/big_test/spice_novlib_homo_sapiens_94_147_tissue", "transcript_data", "transcript_info.json"), "r") as f:
        info_dict = json.load(f)["ENSG00000088682"]

    with open(os.path.join("C:/Users/chris/Desktop/big_test/spice_novlib_homo_sapiens_94_147_tissue", "fas_data", "fas_index.json"), "r") as f:
        filename = json.load(f)["ENSG00000088682"]

    with open(os.path.join("C:/Users/chris/Desktop/big_test/spice_novlib_homo_sapiens_94_147_tissue", "fas_data", "fas_scores", filename), "r") as f:
        distance_dict = json.load(f)["ENSG00000088682"]

    optimizer = RMSDOptimizerAlt(copy.deepcopy(distance_dict),
                                 info_dict,
                                 False,
                                 False)

    optimizer.calc_rmsd()
    print(len(optimizer.matrix))
    print(optimizer.length)
    print(optimizer.get_results().x)
    print(optimizer.objective_function(optimizer.get_results().x))






def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
