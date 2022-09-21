# How to FAS
 
1. Assemble all sequence that FAS scores should be calculated for in a fasta file. (You will need the headers for each sequence)
2. Annotate the entire fasta file using. Note that the output file must have the same name as the input file for FAS to recognize it automatically.

> fas.doAnno \
> -i sequences.fasta \
> -o path/to/annotation/directory/ \
> -t directory/containing/pfam/smart/coils2/and/such/ \
> -n sequences \
> --cpus 16

3. Generate a tsv file with all the pairings you want to calculate. Requires two columns containing the fasta headers of the sequences that FAS scores should be calulcated for.
4. Generate a tsv file with two columns. The first column contains all fasta headers and the second column the the ncbi taxon id, for example ncbi9606 for human.
5. Run fas.run like this:

> fas.run \
> --seed sequences.fasta \
> --query sequences.fasta \
> --annotation_dir path/to/annotation/directory \
> --out_dir /path/to/arbitrary/output/directory \
> --pairwise path/to/tsv/containing/pairings.tsv \
> --out_name arbitrary_output_name \
> --tsv \
> --phyloprofile path/to/tsv/containg/headers/and/taxon_id.tsv \
> --domain \
> --empty_as_1 \

6. Done! This can take quite a while, but you can do all pairings in a single FAS run like this, by generating a file that contains only the pairings you want.