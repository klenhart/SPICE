import json
from typing import List, Dict, Any

from Classes.GTFBoy.AnnotationParser import AnnotationParser


def main():
    # anno_list: List[str] = ["C:/Users/chris/Desktop/Alzheimer/ENCSR094NFM/ENCFF129ZDE.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR169YNI/ENCFF270HYK.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR205QMF/ENCFF573QGI.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR257YUB/ENCFF707NQJ.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR462COR/ENCFF561RPD.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR463IDK/ENCFF635HXQ.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR690QHM/ENCFF519SVS.gtf",
    #                         "C:/Users/chris/Desktop/Alzheimer/ENCSR697ASE/ENCFF385SKN.gtf"]
#
    # anno_parse_list = [AnnotationParser([anno_list[0]],
    #                                     [],
    #                                     0.0,
    #                                     name="cookie" + str(i)) for i, path in enumerate(anno_list)]
    # for anno_parse in anno_parse_list:
    #     anno_parse.parse_annotations()
    #     anno_parse.save_json("C:/Users/chris/Desktop/Alzheimer/test/")

    parsed_list: List[str] = ["C:/Users/chris/Desktop/Alzheimer/test/cookie0.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie1.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie2.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie3.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie4.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie5.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie6.json",
                              "C:/Users/chris/Desktop/Alzheimer/test/cookie7.json"]

    anno_dict_list: List[Dict[str, Any]] = [open_json(path) for path in parsed_list]

    for i, dictionary in enumerate(anno_dict_list):
        print(i)
        for j, query_dict in enumerate(anno_dict_list):
            print("\t", j)
            if i < j:
                for key in dictionary.keys():
                    synonyms: List[str] = dictionary[key]["synonyms"]
                    for query_key in query_dict.keys():
                        if query_key != key:
                            for synonym in synonyms:
                                if synonym in query_dict[query_key]["synonyms"]:
                                    print("--->\norigin:", i, key,
                                          "\nquery",  j, query_key,
                                          "\nsynonym:", synonym,
                                          "\nquery_gene: ", query_dict[query_key]["gene_id"],
                                          "\norigin_gene: ", dictionary[key]["gene_id"])


def open_json(path) -> Dict[str, Any]:
    with open(path, "r") as f:
        mea_dict = json.load(f)
    return mea_dict






if __name__ == "__main__":
    main()
