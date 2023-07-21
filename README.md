# Spice
## Splicing-based Protein Isoform Comparison Estimator
### Applying the [FAS](https://github.com/BIONF/FAS) algorithm to entire transcript sets and putting differential transcript expression data into a context of function.

Typically gene expression levels are used as a proxy to determine protein composition and thereby functional diversity of cells or tissue samples. This disregards the fact that alternative Splicing (AS) of pre-mRNA can disrupt, delete or insert entire protein domains and thereby change the corresponding proteins functional effect. This means that using only gene expression levels as a proxy for functional diversity is not sufficient. To gain further insight about the functional diversity it is necessary to characterize the function of a protein isoform by relating it to it's domain content. Only a small set of proteins undergo experimental characterization, which is why the remaining proteins must be assessed using bioinformatics tools.

Spice is able to assemble all protein sequences of a species as represented in the ensembl genome browser and quantify the functional similarity between each isoform of each protein coding gene using the Feature Architecture Similarity (FAS) algorithm. FAS makes use of databases that contain protein domains like PFAM or TMHMM.

## Table of Contents
* [Requirements](#requirements)
* [Usage](#usage)
  * [Initialize the Spice library](#initialize-the-spice-library)
  * [Annotation](#annotation)
  * [Generate job arrays](#generate-job-arrays)
  * [Parse domain output](#parse-domain-output)
  * [Apply library to expression data](#apply-library-to-expression-data)
    * [Create result directory](#create-result-directory)
    * [Import expression data](#import-expression-data)
    * [Merge samples into condition](#merge-samples-into-condition)
    * [Compare conditions](#compare-conditions)
* [Contact](#contact)

## Requirements

Spice has been implemented for python 3.7. To run the scripts the following python modules are required:

```
greedyFAS
pyranges
requests
json
argparse
itertools
plotly
pathlib
tqdm
numpy
matplotlib
scipy
yaml
pandas
gzip
shutil
```

Since Spice makes use of the FAS algorithm, it should be installed beforehand. Check out the [FAS](https://github.com/BIONF/FAS) repository for further instruction on setting it up.

## Usage

In order to simply get all Spice scripts in this GitHub repository use this set of commands in the terminal:

```
mkdir SPICE
cd SPICE
git clone https://github.com/chrisbluemel/SPICE
```

**IMPORTANT**: Running these scripts without access to a processing cluster will be more or less impossible. You can try to calculate 200k FAS comparisons on your personal computer, but it may take a very long time. Also be aware that the helper scripts to generate job arrays are only capable of generating SLURM job arrays.

### Initialize the Spice library

To initialize the library, use this command:

```
python spice_library.py \
--species human \
--release 107 \
--outdir parent/directory/of/the/library
```

There is further arguments that can be passed to spice_library.py. Descriptions can be accessed by using the --help argument.

### Annotation

After having initialized the library, we need to annotate them. This example was written when the most recent ensembl release was 107. You could have any higher release number. (Greetings from the past!) Run this command: 

```
fas.doAnno \
-i /path/to/spice_lib_homo_sapiens_107_1ee/transcript_data/transcript_set.fasta \
-o /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/ \
-t /path/to/annoTools/ \
-n annotations \
--cpus 16 \
&& \
python \
get_domain_importance.py \
-i /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/annotations.json \
-o /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/ \
&& \
python \
restructure_anno.py \
-i /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/annotations.json \
-o /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/
```

You can use as many or as few CPUs as you have access to, but remember that annotation takes a lot of time and it will take longer, if you have less workers.

**IMPORTANT** The output name under the argument -n must be annotations.

### Generate job arrays

Now that all sequences are collected and annotated we can do the FAS Scoring. It is recommended to run fas.run once per gene and not all genes at once. Create a SLURM job array that references the gene_ids.txt, which was created in the directory /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/. This script will generate all necessary SLURM job-arrays.

```
python \
FASJobAssistant.py \
--Lib_dir /path/to/spice_lib_homo_sapiens_107_1ee/ \
--memory 2 \
--partitions partition1 partition2 etc \
--dir_fas /path/to/FAS/bin/ \
--fas_mode run \
--outdir /path/to/directory/that/shall/contain/job/arrays/
```

Run all of these job-arrays.

### Parse domain output

Once all FAS runs have finished which can take a few days (mostly due to a few very large proteins like TITIN) only one last script needs to be run:

```
python \
parse_domain_out.py \
-f /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/forward.domains \
-r /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/reverse.domains \
-m /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/annotations_map.json \
-o /path/to/spice_lib_homo_sapiens_107_1ee/fas_data/
```


### Apply library to expression data

The library is now finished. Now we want to make use of it. Grand-Trumpet uses the FAS scores for all protein coding   and nonsense-mediated-decay isoforms of a gene and combines them with the expression levels of each isoform to generate a Expression Weighted Functional Disturbance (EWFD) value for each isoform. The EWFD assumes how much of the functional diversity of the gene does not represented transcripts functional effect. 

To calculate the EWFD scores for any number of samples you need to create two text files. One that contains the names of the samples (one name per row) and one that contains the corresponding paths to the expression GTF files.

#### Create result directory

To initialize a result directory from the Spice library execute this script:

```
python \
spice_result.py \
-m setup \
-l /path/to/spice_lib_homo_sapiens_107/ \
-o /path/to/result/parent/directory/
```

#### Import expression data
To import expression gtf files use this command:

```
python \
spice_result.py \
-m expression \
-l /path/to/spice_lib_homo_sapiens_107_1ee \
-o /path/to/parent/directory/of/result \
-n sample1 \
-g /path/to/expression.gtf \
-N FPKM \
-t 1.0
```

The samples will automatically be assumed as single-replicate conditions.

#### Merge samples into condition

To merge several already imported samples into a condition use this command:

```
python \
spice_result.py \
-m condition \
-l /path/to/spice_lib_homo_sapiens_107_1ee \
-o /path/to/parent/directory/of/result \
-n conditionName \
-r sample1 sample2 sample3
```

#### Compare conditions

To compare several already generated conditions use this command:

```
python \
spice_result.py \
-m compare \
-l /path/to/spice_lib_homo_sapiens_107_1ee \
-o /path/to/parent/directory/of/result \
-c condition1;condition2 condition1;condition3 condition2;condition3
```

## Contact

Christian Blümel christian.bluemel@stud.uni-frankfurt.de
