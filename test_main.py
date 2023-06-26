import json
from typing import List, Dict, Any

from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from spice_result import expression


def main():
    argument_parser: ReduxArgParse = ReduxArgParse(["--mode", "--library", "--outdir", "--gtf",
                                                    "--name", "--replicates", "--compared", "--Normalization",
                                                    "--threshold", "--suffix"],
                                                   [str, str, str, str,
                                                    str, str, str, str,
                                                    float, str],
                                                   ["store", "store", "store", "store",
                                                    "store", "store", "store", "store",
                                                    "store", "store"],
                                                   [1, 1, 1, "?",
                                                    "?", "*", "*", "?",
                                                    "?", "?"],
                                                   ["""Mode: Either 'setup', 'expression', 'condition' or 'compare'.
                                                       'setup' creates the result directory system.
                                                       'expression' loads a gtf expression file into results as replicate.
                                                       'condition' combines replicates and calculates the new conditions
                                                       EWFD.
                                                       'compare' generates a comparison between two conditions EWFDs.
                                                       """,
                                                    "Path to the library. Required for all modes.",
                                                    "Parent directory of the result output. Required for all modes.",
                                                    """Path to gtf expression file to be imported.
                                                    Required for 'expression'""",
                                                    """Name of replicate or condition that will be generated. 
                                                    Required for 'expression' and 'condition' modes.""",
                                                    """Names of replicates to be joined into condition.
                                                    Required for 'condition' mode. Must be more than one.""",
                                                    """Names of conditions that shall be compared.
                                                     Required for 'compare' mode.""",
                                                    """Normalization method referencing expression entry in gtf file.
                                                    Usually TPM or FPKM. Required for 'expression' mode.""",
                                                    """Any expression entry below this threshold will be set to 0.0.
                                                    Optional for 'expression' mode. 1.0 by default.""",
                                                    """Suffix of the result directory that will either be manipulated
                                                    or created."""
                                                    ])
    argument_parser.generate_parser()
    argument_parser.execute()
    argument_dict: Dict[str, Any] = argument_parser.get_args()

    if argument_dict["suffix"] is None:
        argument_dict["suffix"] = ""

    if any(argument_dict[key] is None for key in ["outdir", "gtf", "name", "Normalization", "library"]):
        print("""'expression' mode failed. 
        Either 'library', 'outdir', 'gtf', 'Normalization' or 'name' missing from commandline arguments.
        Aborting.""")
    else:
        if argument_dict["threshold"] is None:
            argument_dict["threshold"] = 1.0

        expression(argument_dict["library"][0],
                   argument_dict["outdir"][0],
                   argument_dict["gtf"],
                   argument_dict["name"],
                   argument_dict["Normalization"],
                   argument_dict["threshold"],
                   argument_dict["suffix"])

def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict






if __name__ == "__main__":
    main()
