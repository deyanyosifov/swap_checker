# Swap_checker

<div align="justify">
An R function that can detect sample swaps or sample mislabelling from tNGS (targeted next generation sequencing) data

## Installation

#### Prerequisites
* A working installation of R (version 4.2.0 or more recent) on a computer with a Linux or Windows operating system. (Theoretically MacOS should be possible, too, but I haven't had the chance to test whether it works.)
* The following R packages have to be installed: `maftools` (obligatory) and `pheatmap` (optional, for generating graphs).

#### Installation
Download the zip archive of all files in the repository by clicking on "Code" and then on "Download ZIP" on the GitHub page of the Swap_checker project or by following this direct download link: https://github.com/deyanyosifov/Swap_checker/archive/refs/heads/master.zip. Unzip the archive, this action will create a new directory named `swap_checker-master` in the current directory. You can rename the new directory to `Swap_checker` or whatever other name you choose and move it to a convenient place on your computer. For the purpose of this manual, we will assume that your installation is located in the directory `Swap_checker` in your home folder on a Linux machine, i.e. `~/Swap_checker`. If your situation is different, just replace the `~/Swap_checker` part in the further instructions with the actual path to your installation.

## Usage

#### Preparation of the data
Your data (coordinate-sorted and indexed bam files) should be located in one or several directories.

#### Usage
* First, load the swap_checker() function into the R environment by entering the command `source("~/Swap_checker/Swap_checker.R")`. Modify accordingly if Swap_checker.R is located in another directory.
* Then, call the swap_checker() function using the following syntax: `swap_checker(argument1 = 'value1', argument2 = 'value2', ...)`
* `swap_checker('help')` will display help.

Alternatively, you can use the function directly from a Linux terminal with the help of Rscript:
`Rscript -e \"source('your/installation/path/Swap_checker.R'); swap_checker(argument1 = 'value1', argument2 = 'value2', ...)`

#### Arguments
* bam.dir              Directory or list of directories that contain the BAM files to be processed. The BAM files must be coordinate-sorted and indexed. Both relative and absolute paths are accepted. A list of directories should have the format `c(path/to/directory/with/BAMs/1','/path/to/directory/with/BAMs/2')`. Note that in this case the separate directories have to be enclosed in single quotation marks but the list as a whole not. If not specified, swap_checker will look for BAM files in the current working directory of the R environment. This will not work if swap_checker() is called from a terminal via Rscript - in this case the bam.dir argument must always be specified, e.g. `swap_checker(bam.dir = '.')` if the BAM files are located in the current working directory.
* bam.chrom.notation   Specifies the chromosome notation used in the BAM files. Possible values: `'chr'` (default, to be used when chromosomes are named "chr1", "chr2", etc.) and `'nochr'` (to be used when chromosomes are named "1", "2", etc.).
* bam.pattern          This is an optional argument that helps to limit the analysis to a subset of the available BAM files, based on file name patterns. For example, if bam.dir contains both unsorted and sorted BAM files as different intermediate outputs of an NGS processing pipeline, running swap_checker() will result in error because the function accepts only coordinate-sorted and indexed BAM files. However, if file names of the sorted BAM files end in "sorted.bam", `swap_checker(bam.pattern = '*.sorted.bam$')` will analyze only the sorted BAM files without errors. The value of the "bam.pattern" argument should be a regular expression like in the previous simple example. Regular expressions are a powerful tool and their use allows very precise selection of BAM files that have to be analyzed. Inexperienced users who would wish to use them in more complicated scenarios might want to consult a quick guide to regular expression syntax, such as https://www.rexegg.com/regex-quickstart.html.
* snp.list             BED file containing a custom list of SNPs. You can use the SNPs_selected.bed file as a template for your own list. Make sure that you provide genomic coordinates according to the same genome version that was used for the alignment. If not specified, the program will use a list of 28 SNPs in genes frequently mutated in CLL (SNPs_selected.bed) and will assume that your BAM files are aligned to hg19.
* nthreads             Number of threads to use. Each BAM file will be launched on a separate thread. Works only on Unix and macOS. Default: 4
* samples              Argument with two possible values: `'single'` (default) and `'multiple'`. Specifies whether all BAM files originate from separate individuals ('single'), in which case swap_checker() will be able to identify only unexpectedly matching samples, or there are more than one sample per subject ('multiple'), which will allow swap_checker() to also check whether samples labelled as originating from the same individual match as expected. The latter will work only if BAM file names follow an uniform pattern and contain identifiers for the subjects (see argument \"subject.identifier\").
* subject.identifier   Specifies the part of BAM file names that contains the unique subject identifier. Relevant only when `"samples = 'multiple'"` (see above). This argument accepts one of six preset values or a custom regular expression. The preset values are 'NUMERICID-XXX', 'MIXEDID-XXX', 'XXX-NUMERICID-XXX', 'NUMERICID_XXX', 'MIXEDID_XXX' and 'XXX_NUMERICID_XXX'. 'NUMERICID-XXX' (default) tells the swap_checker() function that file names begin with a numeric sequence that is unique for each subject, followed by a hyphen and the rest of the file name. 'MIXEDID-XXX' is similar but allows alphanumeric characters in the subject identifier. Accordingly, 'XXX-NUMERICID-XXX' instructs the swap_checker() function to look for the first digits-only part of the file name that is surrounded by hyphens and to accept it as unique subject identifier. The next three options are analogous to the already described but are suited for the case when file name parts are separated by underscores instead of hyphens. If none of the preset values corresponds to the file naming convention for the BAM files in bam.dir, you can define the subject-identifying part of the file name by using a regular expression, e.g. `'^(\\\\d+)-.*$'` (this is equivalent to 'NUMERICID-XXX'). Inexperienced users might want to consult a quick guide to regular expression syntax, such as https://www.rexegg.com/regex-quickstart.html.
* results.dir          Directory in which the results will be saved. By default, the results will be saved in the current working directory.
* suffix               Optional suffix that can be appended to the main part of result file names. A time stamp suffix is always added automatically and if an optional suffix is specified, it will appear between the main part of the file name and the time stamp. For example, `swap_checker(suffix = 'ProjectX')` will result in the creation of result files with names like "Pairwise_concordance_ProjectX_20250531_174936.csv" and "SNP_readcounts_ProjectX_20250531_174936.csv".
* plot                 Optional argument with two possible values: `'no'` (default) or `'yes'`. If set to 'yes', `swap_checker()` will produce a correlation plot. Such a plot can be useful only for relatively small cohorts, otherwise it will be too crowded and unreadable if printed on a single page. If sample names (taken from BAM file names) are too long and take up most of the plotting space, you can shorten them by removing redundant file name parts using the `"remove.name.part"` argument (see below).
* remove.name.part     Optional argument related to the plotting function (see above). Its value should be a regular expression defining a part of the sample name that should be removed. For example, if all BAM files end in "_L001_R1_bwa.deduplicated.bam", `swap_checker(plot = 'yes', remove.name.part = '_L001.*$')` will produce a correlation plot and the displayed sample names will not contain the redundant part from the BAM file names.

#### Output files
* Pairwise_concordance_[suffix]_[time stamp].csv       This is the main results file, a table in CSV format with the following columns:
                                                       - "X_bam" and "Y_bam" contain the sample names (taken from the input BAM files) arranged so that each sample in the first column is matched against any of the remaining samples in the second column;
                                                       - "concordant_snps": number of the concordant SNPs;
                                                       - "discordant_snps": number of the discordant SNPs;
                                                       - "discordant_snps_detailed": genomic coordinates of the discordant SNPs;
                                                       - "fract_concordant_snps": fraction of the concordant SNPs;
                                                       - "cor_coef": Pearson's correlation coefficient;
                                                       - "XY_possibly_paired": estimation whether samples originate from the same individual - "Yes" or "No". Positive findings are listed first.
* SNP_readcounts_[suffix]_[time stamp].csv             This is a detailed table with the number of the reads that support the reference or the alternative allele at each of the assessed SNP positions in each sample. Columns:
                                                       - "BAM": sample name (taken from the respective input BAM file);
                                                       - "loci": genomic coordinates of the SNP;
                                                       - "ref_rc": read counts supporting the reference allele;
                                                       - "alt_rc": read counts supporting the alternative allele;
                                                       - "vaf": variant allele fraction of the alternative allele.
* Expected_concordant_pairs_[suffix]_[time stamp].csv  This is an optional subset of the main table, created if "samples" was set to 'multiple'. It contains only the rows for combinations of samples that are expected to be concordant because they should originate from the same individual.
* Pairwise_concordance_[suffix]_[time stamp].png       This is an optional correlation plot (created if "plot" was set to 'yes') based on the correlation coefficients from the main results table.
