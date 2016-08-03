# PEGASUS

Pegasus Version 1.3  
Priyanka Nakka, Ramachandran Lab, Brown University 

June 23, 2016

Code is modified from VEGAS source code (Liu et al AJHG 2010) 

Questions: Contact priyanka_nakka@brown.edu 

Requirements:

Perl 5

R (version 3.0.2 or higher)

Plink 1.07 (1.9 beta 3, 7 Jun is also okay) 

Please note the following:
1. R and PLINK must be installed and in your PATH for this program to run. 
2. This program requires the R packages corpcor and CompQuadForm.  If they are missing, it will return an error.

This directory (PEGASUS) contains the following: 

1. pegasus.pl - source code for PEGASUS program. 

2. 1kgEUR - directory containing genotype data for the 1000 Genomes (Phase 3) EUR super-population and gene lists needed for this script to function.  Please do not delete this directory even if you are using your own reference genotypes.
3. glist-hg19 - hg19 gene list from the UCSC genome browser

4. example.txt - example input pvalues file

5. example.out - example output file

Usage: 
Commands to run this program should have the following format: 

./pegasus.pl [input pvalues file] [parameters]

The input p-values file is a tab-separated table of SNP names and corresponding p-values. For an example file of p-values, please see example.txt.  The one required parameter is to specify which reference genotype data or custom LD file should be used for LD calculations.  If using the provided 1000 Genomes EUR data, type: "-pop 1kgEUR". If using your own genotype data to calculate LD, type "-custom" followed by the path to and prefix of your PLINK bed files. Alternatively, if you would like to use your own LD file containing pre-calculated values for correlation between SNPs in every gene in the dataset, type “-ld-file” followed by the absolute path to your LD file in PLINK format. Note that correlation values should be r and not r^2. (See below for further explanation and sample commands.) 
  
This program uses the following flags for parameters:
Required:

-custom [path + prefix of PLINK bed files for custom genotypes to use for calculating LD, must be in bed format]     specify path and name of plink bed/bim/fam files with custom genotypes to calculate LD matrices
OR

-pop 1kgEUR	use given 1000 Genomes EUR data as reference genotypes to calculate LD
OR 

-ld-file [absolute path + name of pre-calculated LD file]	specify absolute path and name of pre-calculated LD file in default PLINK format (see below) containing correlation values (r) between SNPs in each gene to be evaluated.

Accepted format for custom LD file: 
CHR_A         BP_A           SNP_A  CHR_B         BP_B           SNP_B            R 

     1       792429       rs3094315      1       819185       rs4040617     0.872127 

     1       792429       rs3094315      1       825852       rs2980300     0.804592 

Please note that processing the custom LD file can slow down the PEGASUS program considerably for large genome-wide SNP datasets.  Parallelizing gene score computation by running each chromosome separately is highly recommended with this option. 

Optional:
-out [path + prefix of outfile]	specify name for the .out output file with gene scores

-chr [# between 1 and 23]  compute gene scores for given chromosome only

-upper [# of bp downstream of gene to be included ex. 30000]	    This flag can only be used with custom genotypes or pre-calculated LD files.  The default is 50kb. 

-lower [# of bp upstream of gene to be included ex. 30000]	    This flag can only be used with custom genotypes or pre-calculated LD files.  The default is 50kb. 

The following is a sample command using the given example files: 

./pegasus.pl example.txt -pop 1kgEUR -out test

The resulting output file test.out should match the given example.out file.

To parallelize the process for large datasets, it is useful to run each chromosome separately using the -chr flag as in the following example:

./pegasus.pl example.txt -pop 1kgEUR -chr 1 -out example_chr1

Troubleshooting:
- When using custom genotypes, please make sure the length of the .bim file is the same as the length of the p-values file (they should contain the same markers).

- "ERROR: Could not place marker for left/right window edge" - this is a plink error message that occurs when there are not enough markers in a given bp range and can be ignored.

- Make sure the directory 1kgEUR is in the same directory as pegasus.pl — this directory contains files needed for the program to function. 