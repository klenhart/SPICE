import json
import os

from Classes.API.ensembl_mod.LocalEnsembl import LocalEnsembl
from Classes.PassPath.PassPath import PassPath
from Classes.SequenceHandling.GeneAssembler import GeneAssembler


def main():
    local_ensembl: LocalEnsembl = LocalEnsembl("human",
                                               "C:\\Users\\chris\\Desktop\\git\\test01",
                                               "107")

    gene_assembler: GeneAssembler = GeneAssembler(local_ensembl.get_species_name(),
                                                  local_ensembl.get_taxon_id())

    library_name: str = "spice_lib_" + local_ensembl.get_species_name() + "_" + "107"

    with open(os.path.join("C:\\Users\\chris\\Desktop\\git\\test01", library_name, "paths.json"), "r") as f:
        pass_path: PassPath = PassPath(json.load(f))

    gene_assembler.load(pass_path)
    gene_assembler.integrate_fas_json("C:\\Users\\chris\\Desktop\\AKE\\data\\stats\\old_fas.json")
    gene_assembler.save_fas(pass_path)


if __name__ == "__main__":
    main()
