import json
import os
from typing import List, Dict, Any

from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():
    path: str = ""
    for filename in os.listdir(path):
        new_file: List[str] = list()
        with open(os.path.join(path, filename), "r") as f:
            file = f.read()
        for i, line in enumerate(file.split("\n")):
            if len(line) == 0:
                continue
            elif i == 0:
                new_file.append(line)
            elif line.split(",")[7] == "Yes":
                new_file.append(line)

        file = "\n".join(new_file)
        with open(os.path.join(path, filename), "w") as f:
            f.write(file)


def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
