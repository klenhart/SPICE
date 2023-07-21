import copy
import json
import os
from typing import List, Dict, Any


from Classes.PassPath.PassPath import PassPath
from Classes.ReduxArgParse.ReduxArgParse import ReduxArgParse
from Classes.SequenceHandling.GeneAssembler import GeneAssembler
from spice_result import expression


def main():
    x = "h9_definitive_endoderm\;h9_pancreatic_progenitors h9_pancreatic_progenitors\;liver h9_chondrocytes\;h9_pancreatic_progenitors h9_neural_crest\;h9_pancreatic_progenitors h9_osteocytes\;h9_pancreatic_progenitors adrenal_gland\;h9_pancreatic_progenitors h9_pancreatic_progenitors\;ovary h9_pancreatic_progenitors\;k562 h9_pancreatic_progenitors\;hepg2 h9_ESCs\;h9_pancreatic_progenitors colon\;h9_pancreatic_progenitors h9_b_pancreatic_cell\;h9_pancreatic_progenitors h9_pancreatic_progenitors\;lung GM12878\;h9_pancreatic_progenitors h9_pancreatic_progenitors\;heart h9_definitive_endoderm\;h9_pancreatic_progenitors h9_definitive_endoderm\;liver h9_chondrocytes\;h9_definitive_endoderm h9_definitive_endoderm\;h9_neural_crest h9_definitive_endoderm\;h9_osteocytes adrenal_gland\;h9_definitive_endoderm h9_definitive_endoderm\;ovary h9_definitive_endoderm\;k562 h9_definitive_endoderm\;hepg2 h9_ESCs\;h9_definitive_endoderm colon\;h9_definitive_endoderm h9_b_pancreatic_cell\;h9_definitive_endoderm h9_definitive_endoderm\;lung GM12878\;h9_definitive_endoderm h9_definitive_endoderm\;heart h9_pancreatic_progenitors\;liver h9_definitive_endoderm\;liver h9_chondrocytes\;liver h9_neural_crest\;liver h9_osteocytes\;liver adrenal_gland\;liver liver\;ovary k562\;liver hepg2\;liver h9_ESCs\;liver colon\;liver h9_b_pancreatic_cell\;liver liver\;lung GM12878\;liver heart\;liver h9_chondrocytes\;h9_pancreatic_progenitors h9_chondrocytes\;h9_definitive_endoderm h9_chondrocytes\;liver h9_chondrocytes\;h9_neural_crest h9_chondrocytes\;h9_osteocytes adrenal_gland\;h9_chondrocytes h9_chondrocytes\;ovary h9_chondrocytes\;k562 h9_chondrocytes\;hepg2 h9_ESCs\;h9_chondrocytes colon\;h9_chondrocytes h9_b_pancreatic_cell\;h9_chondrocytes h9_chondrocytes\;lung GM12878\;h9_chondrocytes h9_chondrocytes\;heart h9_neural_crest\;h9_pancreatic_progenitors h9_definitive_endoderm\;h9_neural_crest h9_neural_crest\;liver h9_chondrocytes\;h9_neural_crest h9_neural_crest\;h9_osteocytes adrenal_gland\;h9_neural_crest h9_neural_crest\;ovary h9_neural_crest\;k562 h9_neural_crest\;hepg2 h9_ESCs\;h9_neural_crest colon\;h9_neural_crest h9_b_pancreatic_cell\;h9_neural_crest h9_neural_crest\;lung GM12878\;h9_neural_crest h9_neural_crest\;heart h9_osteocytes\;h9_pancreatic_progenitors h9_definitive_endoderm\;h9_osteocytes h9_osteocytes\;liver h9_chondrocytes\;h9_osteocytes h9_neural_crest\;h9_osteocytes adrenal_gland\;h9_osteocytes h9_osteocytes\;ovary h9_osteocytes\;k562 h9_osteocytes\;hepg2 h9_ESCs\;h9_osteocytes colon\;h9_osteocytes h9_b_pancreatic_cell\;h9_osteocytes h9_osteocytes\;lung GM12878\;h9_osteocytes h9_osteocytes\;heart adrenal_gland\;h9_pancreatic_progenitors adrenal_gland\;h9_definitive_endoderm adrenal_gland\;liver adrenal_gland\;h9_chondrocytes adrenal_gland\;h9_neural_crest adrenal_gland\;h9_osteocytes adrenal_gland\;ovary adrenal_gland\;k562 adrenal_gland\;hepg2 adrenal_gland\;h9_ESCs adrenal_gland\;colon adrenal_gland\;h9_b_pancreatic_cell adrenal_gland\;lung GM12878\;adrenal_gland adrenal_gland\;heart h9_pancreatic_progenitors\;ovary h9_definitive_endoderm\;ovary liver\;ovary h9_chondrocytes\;ovary h9_neural_crest\;ovary h9_osteocytes\;ovary adrenal_gland\;ovary k562\;ovary hepg2\;ovary h9_ESCs\;ovary colon\;ovary h9_b_pancreatic_cell\;ovary lung\;ovary GM12878\;ovary heart\;ovary h9_pancreatic_progenitors\;k562 h9_definitive_endoderm\;k562 k562\;liver h9_chondrocytes\;k562 h9_neural_crest\;k562 h9_osteocytes\;k562 adrenal_gland\;k562 k562\;ovary hepg2\;k562 h9_ESCs\;k562 colon\;k562 h9_b_pancreatic_cell\;k562 k562\;lung GM12878\;k562 heart\;k562 h9_pancreatic_progenitors\;hepg2 h9_definitive_endoderm\;hepg2 hepg2\;liver h9_chondrocytes\;hepg2 h9_neural_crest\;hepg2 h9_osteocytes\;hepg2 adrenal_gland\;hepg2 hepg2\;ovary hepg2\;k562 h9_ESCs\;hepg2 colon\;hepg2 h9_b_pancreatic_cell\;hepg2 hepg2\;lung GM12878\;hepg2 heart\;hepg2 h9_ESCs\;h9_pancreatic_progenitors h9_ESCs\;h9_definitive_endoderm h9_ESCs\;liver h9_ESCs\;h9_chondrocytes h9_ESCs\;h9_neural_crest h9_ESCs\;h9_osteocytes adrenal_gland\;h9_ESCs h9_ESCs\;ovary h9_ESCs\;k562 h9_ESCs\;hepg2 colon\;h9_ESCs h9_ESCs\;h9_b_pancreatic_cell h9_ESCs\;lung GM12878\;h9_ESCs h9_ESCs\;heart colon\;h9_pancreatic_progenitors colon\;h9_definitive_endoderm colon\;liver colon\;h9_chondrocytes colon\;h9_neural_crest colon\;h9_osteocytes adrenal_gland\;colon colon\;ovary colon\;k562 colon\;hepg2 colon\;h9_ESCs colon\;h9_b_pancreatic_cell colon\;lung GM12878\;colon colon\;heart h9_b_pancreatic_cell\;h9_pancreatic_progenitors h9_b_pancreatic_cell\;h9_definitive_endoderm h9_b_pancreatic_cell\;liver h9_b_pancreatic_cell\;h9_chondrocytes h9_b_pancreatic_cell\;h9_neural_crest h9_b_pancreatic_cell\;h9_osteocytes adrenal_gland\;h9_b_pancreatic_cell h9_b_pancreatic_cell\;ovary h9_b_pancreatic_cell\;k562 h9_b_pancreatic_cell\;hepg2 h9_ESCs\;h9_b_pancreatic_cell colon\;h9_b_pancreatic_cell h9_b_pancreatic_cell\;lung GM12878\;h9_b_pancreatic_cell h9_b_pancreatic_cell\;heart h9_pancreatic_progenitors\;lung h9_definitive_endoderm\;lung liver\;lung h9_chondrocytes\;lung h9_neural_crest\;lung h9_osteocytes\;lung adrenal_gland\;lung lung\;ovary k562\;lung hepg2\;lung h9_ESCs\;lung colon\;lung h9_b_pancreatic_cell\;lung GM12878\;lung heart\;lung GM12878\;h9_pancreatic_progenitors GM12878\;h9_definitive_endoderm GM12878\;liver GM12878\;h9_chondrocytes GM12878\;h9_neural_crest GM12878\;h9_osteocytes GM12878\;adrenal_gland GM12878\;ovary GM12878\;k562 GM12878\;hepg2 GM12878\;h9_ESCs GM12878\;colon GM12878\;h9_b_pancreatic_cell GM12878\;lung GM12878\;heart h9_pancreatic_progenitors\;heart h9_definitive_endoderm\;heart heart\;liver h9_chondrocytes\;heart h9_neural_crest\;heart h9_osteocytes\;heart adrenal_gland\;heart heart\;ovary heart\;k562 heart\;hepg2 h9_ESCs\;heart colon\;heart h9_b_pancreatic_cell\;heart heart\;lung GM12878\;heart"
    y = x.split(" ")
    super = list()
    print(len(y))
    for entry in y:
        both = entry.split("\\;")
        both.sort()
        if both not in super:
            super.append(both)
    print(len(super))

    final = list()
    for entry in super:
        final.append("\\;".join(entry))
    print(" ".join(final))



def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict


if __name__ == "__main__":
    main()
