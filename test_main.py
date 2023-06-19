import json


def main():
    new = "C:/Users/chris/Desktop/spice_novel/test_merge.gtf"
    memory = set()
    problem = set()
    with open(new, "r") as f:
        prev_transcript_id = ""
        for line in f:
            if line.startswith("#"):
                continue
            curr_transcript_id = line.split("\t")[-1].split(";")[1].split(" ")[2][1:-1]
            if curr_transcript_id != prev_transcript_id:
                memory.add(prev_transcript_id)
                if curr_transcript_id in memory:
                    problem.add(curr_transcript_id)
            prev_transcript_id = curr_transcript_id
    print(problem)


if __name__ == "__main__":
    main()
