# Spice
## Splicing-based Protein Isoform Comparison Estimator
### Applying the [FAS](https://github.com/BIONF/FAS) algorithm to entire transcript sets and putting differential transcript expression data into a context of function.

Typically gene expression levels are used as a proxy to determine protein composition and thereby functional diversity of cells or tissue samples. This disregards the fact that alternative Splicing (AS) of pre-mRNA can disrupt, delete or insert entire protein domains and thereby change the corresponding proteins functional effect. This means that using only gene expression levels as a proxy for functional diversity is not sufficient. To gain further insight about the functional diversity it is necessary to characterize the function of a protein isoform by relating it to it's domain content. Only a small set of proteins undergo experimental characterization, which is why the remaining proteins must be assessed using bioinformatics tools.

Grand-Trumpet is able to assemble all protein sequences of a genome as represented in the ensembl genome browser and quantify the functional similarity between each isoform of each protein coding gene using the Feature Architecture Similarity (FAS) algorithm. FAS makes use of databases that contain all known protein domains like PFAM.

Additionally the pipeline is able to apply the quantification of functional similarity to a set of expression data, thereby generating a map of functional diversity among all genes while also being able to visualize the difference of functional diversity between two given samples.


## Table of Contents
* [Requirements](#requirements)
* [Usage](#usage)
  * [Generate FAS library](#generate-fas-library)
    * [Initialize the FAS library and download local ensembl release](#initialize-the-fas-library-and-download-local-ensembl-release)
    * [Install FAS](#install-fas)
    * [Collect sequences](#collect-sequences)
    * [Annotation](#annotation)
    * [Generate job arrays](#generate-job-arrays)
    * [Run FAS](#run-fas)
    * [Concatenate FAS output](#concatenate-fas-output)
  * [Apply library to expression data](#apply-library-to-expression-data)
    * [Calculate comparison between pair of samples](#calculate-comparison-between-pair-of-samples)
    * [Identify genes of interest](#identify-genes-of-interest)
    * [Visualize comparison](#visualize-comparison)
* [Contact](#contact)

## Requirements

Grand-Trumpet has been implemented for python 3.7. To run the scripts the following python modules are required:

```
argparse
gzip
itertools
json
plotly
pyranges
requests
shutil
urllib
```

Since Grand-Trumpet makes use of the FAS algorithm, it should be installed beforehand. Check out the [FAS](https://github.com/BIONF/FAS) repository for further instruction on setting it up.

## Usage

In order to simply get all Grand-Trumpet scripts in this GitHub repository use this set of commands in the terminal:

```
mkdir grand-trumpet
cd grand-trumpet
git clone https://github.com/chrisbluemel/grand-trumpet
```

**IMPORTANT**: Running these scripts without access to a processing cluster will be more or less impossible. You can try to calculate one million FAS comparisons on your home computer, but do not say I didn't warn you. Also be aware that the helper scripts to generate job arrays are only capable of generate SLURM job arrays.

### Generate FAS library

First decide on a location for your FAS library. When finished the library can take up several GB of space.

#### Initialize the FAS library and download local ensembl release

Now we need to download a local ensembl assembly for the required species. This will also create the entire directory structure and initialize the config.tsv file. As an output choose the directory that is supposed to contain the finished library. This is what an example command for initialiazing a library for human would look like: 

```
python fas_lib.py -l \
--output /parent/directory/of/the/library \
--species human
```

The path the **config.tsv** file, which can be found at
```
/parent/directory/of/the/library/homo_sapiens/release-num/config.tsv
```
will be necessary as an input for all further steps of the pipeline.


#### Collect sequences

Now we will collect the actual protein sequences from ensembl by using this command:

```
python fas_lib.py \
--config /parent/directory/of/the/library/homo_sapiens/release-num/config.tsv
```

This step has significant RAM requirements and should not be run on any old potato, but rather using a decently powerful machine or even a processing cluster. The latter will be necessary for the later steps anyway. For the human genome this one is done in about three hours.

**IMPORTANT**: Sometimes it might happen that while you are collecting sequences your internet connection drops out or the ensembl servers have a major hiccup. For this reason Grand-Trumpet remembers what it was doing, if it crashes during sequence collection. Just run the command again and the requests will continue right where they were interrupted.

#### Annotation

After having collected all the sequences, we need to annotate them. This example was written when the most recent ensembl release was 107. You could have any higher release number. (Greetings from the past!) Run this command: 

```
fas.doAnno \
-i /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/isoforms.fasta \
-o /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/annotation/ \
-n isoforms \
--cpus 16
```

You can use as many or as few CPUs as you have access to, but remember that annotation takes a lot of time and it will take longer, if you have less workers.

**IMPORTANT** The output name under the argument -n must be the same as the name of the fasta-file. This is why in my case I used isoforms as the -n argument. 

#### Generate job arrays

Now that all sequences are collected and annotated we can do the FAS Scoring. It is recommended to run fas.run once per gene and not all genes at once. Create a SLURM job array that references the gene_ids.txt, which was created in the directory /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/. This is an example job array to calculate the FAS-Scores for the first 1000 genes in the gene_ids.txt: 

```
#!/bin/bash

#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name="fas_human"
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
#SBATCH --array=1-1000

gene=$(awk FNR==$SLURM_ARRAY_TASK_ID "/parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/gene_ids.txt")
python fas_handler.py \
--maketsv \
--gene $gene \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
&& \
fas.run \
--seed /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/isoforms.fasta \
--query /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/isoforms.fasta \
--annotation_dir /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/annotation/ \
--out_dir /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/FAS_buffer/ \
--bidirectional \
--pairwise /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/tsv_buffer/$gene.tsv \
--out_name $gene \
--tsv \
--phyloprofile /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/phyloprofile_ids.tsv \
--domain \
--empty_as_1 \
; \
python fas_handler.py \
--remove \
--gene $gene \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
```
Since creating 20 of those files by hand to wade through all roughly 20000 protein coding genes in the human genome is pretty tedious, I created a helper script that automatically generates one job array for every 1000 protein coding genes in the library. For human this might turn into about 20 job arrays, for other species it might be less or more depending on the gene count. To automatically generate the arrays use this command:

```
python fas_bashAssist.py \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
--partitions all ni nini ninini
--memory 2
```

The finished job arrays can be found in in the library in this directory:

```
/parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/SLURM/
```


#### Run FAS

Now that we have either created the job arrays by hand or by using the fas_bashAssist.py script, we can start running FAS. Remember that many processing clusters have limits on how many jobs you can commit at once. I made Grand-Trumpet work in such a way that it usually automatically prevents the case of too many files being generated at once, which could overburden a file system easily. During the FAS runs I did not implement such a fail save. This is why you should do the next step several times before having finished to run FAS on all genes.


#### Concatenate FAS Output

Run this script several times. It concatenates all FAS output that is currently done into one distance_master.phyloprofile file. Use this command:


```
python fas_handler.py \
--join \
--config /parent/directory/of/the/library/FAS_library/homo_sapiens/release-107/config.tsv \
```

### Apply library to expression data

The library is now finished. Now we want to make use of it. Grand-Trumpet uses the FAS scores for all protein coding  isoforms of a gene and combines them with the expression levels of each isofrom to generate a Movement score for each isoform. The movement score assumes how much of the functional diversity of the isoform is represented in the transcript composition of the gene. 

To calculate the Movement scores for any number of samples you need to make two textfiles. One that contains the names of the samples (one name per row) and one that contains the corresponding paths to the expression GTF files. Several expression files can be joined together by writing them into the same line and seperating them with a semicolon (;). At the moment the expression files will be joined by calculating the mean of all FPKM values between expression files.

Here is an example for the name path text file:


```
sample1
sample2
sample3
```

And here an example for the expression path text file:

```
path/to/sample1_replicate1.gtf;path/to/sample1_replicate2.gtf
path/to/sample2_replicate1.gtf
path/to/sample3_replicate1.gtf;path/to/sample3_replicate2.gtf;path/to/sample3_replicate3.gtf
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
