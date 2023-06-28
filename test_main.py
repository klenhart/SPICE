import json
from typing import List, Dict, Any

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from spice_result import expression


def main():
    print("Hello World!")


def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
