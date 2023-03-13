from Classes.SequenceHandling.LibraryInfo import LibraryInfo


def manipulate(info):
    info["fade_version"] = 2.0
    print("----")


def main():
    info: LibraryInfo = LibraryInfo("C:\\Users\\chris\\Desktop\\git\\fade_lib_homo_sapiens_107\\info.yaml")
    print(info)
    manipulate(info)
    print(info)


if __name__ == "__main__":
    main()
