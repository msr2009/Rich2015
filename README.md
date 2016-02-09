# Rich2016

These are scripts from 

<b>Matthew S Rich, Celia Payen, Alan F Rubin, Giang Ong, Monica R. Sanchez, Nozomu Yachie, Maitreya J Dunham, and Stanley Fields. "The effects of cis-regulatory mutations in the SUL1 gene on yeast sulfate-limited fitness" Genetics (submitted)</b>

Included in this repository are:

`enrich/` -- this folder contains scripts from the Enrich software package which were used to all initial analyses of barcode counts. These scripts are preliminary versions of the Enrich software published in [1], which have more complex normalization and regression functionality.

`plot_heatmap_enrich2.py` -- this script takes output of enrich's barcode frequency analysis and creates the plots found in the Figure 2.

`6mut_frequencies.py` -- this script genotypes paired-end reads from combinatorial mutations experiments. The output of this script is then analyzed using a simple bash workflow:
	
	sort OUTPUT | uniq -c | awk '{OF="\t"}{print $2, $1}'

which creates files to use with 6mut_analysis.R.

`6mut_analysis.R` -- R script to perform analysis of combinatorial mutations experiments. The workflow of this analysis is the same as that found in `enrich/.`  

`enumerate_sequences.py` -- this script enumerates all single mutations of a given sequence, and was used to create the <em>in silico</em> sequences used for cryptic transcription factor binding site analysis (Figure X).

`fimo_newsites.py` -- this script parses FIMO output to determine the statistic (log2(p_wt)/log2(p_mut)) used in Figure S5 for every matching motif.

`p_diff-v-fitness.R` -- this script creates plot found in Figure S5, including coloring points based on their annotations from this paper.

`fastq_tools.py` -- A small library of Python utilities to work with FASTQ files.

`RICH2016_FILE-S2.bz2` -- A bzipped file containing all raw data files (tab-delimited tables of data for each promoter variant; the output of Enrich)

Raw reads used in the analyses in this manuscript can be found at SRA Bioproject ID PRJNA273419.

[1] -- Alan F Rubin, Terry Speed, and Douglas M Fowler. in prep.
