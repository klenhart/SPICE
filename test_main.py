import biomart

def main():
     print("Hello World!")
     import requests, sys

     server = "https://rest.ensembl.org"
     ext = "/sequence/id/94/ENST00000617307?"
     #ext = "/sequence/id/ENSP00000482090?"

     r = requests.get(server + ext, headers={"Content-Type": "text/plain"})

     if not r.ok:
         r.raise_for_status()
         sys.exit()
     # 2384
     print(r.text)
     print(len(r.text))


if __name__ == "__main__":
    main()
