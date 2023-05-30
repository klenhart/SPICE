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
    * [Identify genes of interest](#identify-genes-of-interest)
    * [Visualize comparison](#visualize-comparison)
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

Here is an example for the name path text file:


```
sample1
sample2
sample3
```

And here an example for the expression path text file:

```
path/to/sample1_replicate1.gtf
path/to/sample2_replicate1.gtf
path/to/sample3_replicate1.gtf
```

You can then calculate the Movement for all expression files in the text files by running this command:

```
python fas_movement.py \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
--namepath /path/to/name.txt \
--expressionpath /path/to/expressionpath.txt 
```
The Movement scores will be saved as tsv files within the libraries FAS_polygon directory


#### Calculate comparison between pair of samples

If you have calculated Movement scores for at least two samples you can generate a comparison between them. At the moment comparisons between more than two samples are not possible. Simply use this command:

```
python fas_compare.py \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
--input path/to/polygonFAS_sample1.tsv path/to/polygonFAS_sample2.tsv
```

The input parameter takes paths to the files generated in the previous step. The fas_compare.py script generates another tsv file in the pictures directory of the library. For the example it would be called **sample1xsample2.tsv**. Additionally in will generate a sorted and filtered version of the output file, which will be important for the next step.

#### Identify genes of interest

Now the raw output of the comparison does not give a lot of insight about which comparisons might have yielded something interesting. You can ofcourse just visualize any gene that you are interested in - maybe you even have information from previous analysis of the data that you wanna get another perspective on. But to identify genes of interest independent of any other information we need the sorted and filtered version of the file. It has been automatically generated in the last step. It only contains genes for which any of the two samples had at least two isoforms and the roots-mean-square-deviation (RMSD) between Movement scores was higher than 0 and lower than 1. Additionally the file is sorted by unscaled RMSD first, and then by scaled RMSD. You can now manually scroll through the file from top to bottom. Look for entries where the scaled and unscaled RMSD are not equal, but still both greater than 0.5. 

#### Visualize comparison

To visualize a comparison you are interested in, simply use this command:
```
python fas_visualize.py \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
--gene ENSG00000000003 \
--path /path/to/comparison.tsv
--outFormat svg
```

The scaled aswell as the unscaled version of the comparison will be output into the pictures directory of the library.

## Contact

Christian Blümel christian.bluemel@stud.uni-frankfurt.de
