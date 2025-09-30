# Spice
## Splicing-based Protein Isoform Comparison Estimator
### Applying the [FAS](https://github.com/BIONF/FAS) algorithm to entire transcript sets and putting differential transcript expression data into a context of function.

Typically gene expression levels are used as a proxy to determine protein composition and thereby functional diversity of cells or tissue samples. This disregards the fact that alternative Splicing (AS) of pre-mRNA can disrupt, delete or insert entire protein domains and thereby change the corresponding proteins functional effect. This means that using only gene expression levels as a proxy for functional diversity is not sufficient. To gain further insight about the functional diversity it is necessary to characterize the function of a protein isoform by relating it to it's domain content. Only a small set of proteins undergo experimental characterization, which is why the remaining proteins must be assessed using bioinformatics tools.

Spice is able to assemble all protein sequences of a species as represented in the ensembl genome browser and quantify the functional similarity between each isoform of each protein coding gene using the Feature Architecture Similarity (FAS) algorithm. FAS makes use of databases that contain protein domains like PFAM or TMHMM.

## Table of Contents
- [Spice](#spice)
  - [Splicing-based Protein Isoform Comparison Estimator](#splicing-based-protein-isoform-comparison-estimator)
    - [Applying the FAS algorithm to entire transcript sets and putting differential transcript expression data into a context of function.](#applying-the-fas-algorithm-to-entire-transcript-sets-and-putting-differential-transcript-expression-data-into-a-context-of-function)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Initialize the Spice library](#initialize-the-spice-library)
      - [Create result directory](#create-result-directory)
      - [Import expression data](#import-expression-data)
      - [Merge samples into condition](#merge-samples-into-condition)
      - [Compare conditions](#compare-conditions)
  - [Contact](#contact)

## Installation

In order to simply get all Spice scripts in this GitHub repository use this set of commands in the terminal:

```bash
# Go to the folder where SPICE should be downloaded
cd your/src/folder
git clone https://github.com/klenhart/SPICE
# Go to SPICE folder
cd SPICE
# Create Conda Environment
conda create -n spice_env
conda activate spice_env
pip install -e .

```
If you want to be able to run SPICE modules from any working directory (as long as the `spice_env` environment is active), you can add a path configuration file (.pth) inside your conda environment:
```bash
echo "/path/to/SPICE" > /home/user/miniconda3/envs/spice_env/lib/python3.13/site-packages/spice.pth
```

This tells Python to always include your local SPICE folder in its import path whenever `spice_env` is active.
Without this step, you would need to either:
- run scripts from inside the SPICE folder, or
- set `PYTHONPATH=/path/to/SPICE` manually when running commands.


### Initialize the Spice library

There's now a pipeline available for generating the library. Please visit https://github.com/felixhaidle/spice-nf.


#### Create result directory

To initialize a result directory from the Spice library execute this script:

```bash
python \
spice_result.py \
-m setup \
-l /path/to/spice_lib_homo_sapiens_107/ \
-o /path/to/result/parent/directory/
```

#### Import expression data
To import expression gtf files use this command:

```bash
python -m spice_result \
--mode expression \
--library /path/to/spice_lib_homo_sapiens_107_1ee \
--outdir /path/to/parent/directory/of/result \
--name sample1 \
--gtf /path/to/expression.gtf \
--Normalization FPKM \
--threshold 0.0
```

The samples will automatically be assumed as single-replicate conditions.

**IMPORTANT**: The --threshold argument is deprecated, please use --threshold 0.0.
#### Merge samples into condition

To merge several already imported samples into a condition use this command:

```bash
python -m spice_result \
--mode condition \
--library /path/to/spice_lib_homo_sapiens_107_1ee \
--outdir /path/to/parent/directory/of/result \
--name conditionName \
--replicates sample1 sample2 sample3
```

#### Compare conditions

To compare several already generated conditions use this command:

```bash
python -m spice_result.py \
--mode compare \
--library /path/to/spice_lib_homo_sapiens_107_1ee \
--outdir /path/to/parent/directory/of/result \
--compared condition1;condition2 condition1;condition3 condition2;condition3
```

## Contact

Christian Blümel christian.bluemel@stud.uni-frankfurt.de
Katharina Lenhart k.lenhart@outlook.com
