import copy
import json
import os
from typing import List, Dict, Any


from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():
    print(len("MAALKSWLSRSVTSFFRYRQCLCVPVVANFKKRCFSELIRPWHKTVTIGFGVTLCAVPIAQ"))



def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
