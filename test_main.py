import json

from Classes.GTFBoy.GTFBoy import GTFBoy


def update_tags(tag_list):
    # Search for tags indicating a 3'/5' incomplete transcript.
    if any([tag in GTFBoy.start_incomplete_tags for tag in tag_list]):
        tag_list.append("start_incomplete")
        tag_list.append("incomplete")
    elif any([tag in GTFBoy.end_incomplete_tags for tag in tag_list]):
        tag_list.append("end_incomplete")
        tag_list.append("incomplete")
    else:
        tag_list.append("complete")
    return tag_list


def main():
    path: str = "C:/Users/chris/Desktop/git/spice_lib_homo_sapiens_107/transcript_data/transcript_set.json"
    with open(path, "r") as f:
        mega_dict = json.load(f)
    for gene_id in mega_dict.keys():
        for transcript_id in mega_dict[gene_id]["transcripts"].keys():
            mega_dict[gene_id]["transcripts"][transcript_id]["tags"] = update_tags(mega_dict[gene_id]["transcripts"][transcript_id]["tags"])
    with open(path, "w") as f:
        json.dump(mega_dict, f, indent=4)


if __name__ == "__main__":
    main()
