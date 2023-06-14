
import requests, sys


def main():

    server = "https://rest.ensembl.org"
    ext = "/sequence/id/94/ENST00000618181?"

    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    print(r.text)


if __name__ == "__main__":
    main()
