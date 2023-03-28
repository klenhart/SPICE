from Classes.SequenceHandling.LibraryInfo import LibraryInfo
import requests
import json
import sys


def manipulate(info):
    info["fade_version"] = 2.0
    print("----")


def main():
    x = [1, 2, 3]
    print(x)
    x.pop(1)
    print(x)

if __name__ == "__main__":
    main()
