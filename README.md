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
First load the swap_checker() function into the R environment by entering the command source("~/Swap_checker/Swap_checker.R"). Modify accordingly if Swap_checker.R is located in another directory.
Then, call the swap_checker() function with the respective arguments (see below): swap_checker(argument1 = 'value1', argument2 = 'value2', ...)
swap_checker('help') will display help.

Alternatively, you can use the function directly from a Linux terminal with the help of Rscript: Rscript -e \"source('your/installation/path/Swap_checker.R'); swap_checker(argument1 = 'value1', argument2 = 'value2', ...)

Arguments:
bam.dir              Directory or list of directories that contain the BAM files to be processed. The BAM files must be coordinate-sorted\n",
                     and indexed. Both relative and absolute paths are accepted. A list of directories should have the format\n",
                     c('/path/to/directory/with/BAMs/1','/path/to/directory/with/BAMs/2'). Note that in this case the separate directories\n",
                     have to be enclosed in single quotation marks but the list as a whole not.\n",
                     If not specified, swap_checker will look for BAM files in the current working directory of the R environment. This will\n",
                     not work if swap_checker is called from a terminal via Rscript - in this case the bam.dir argument must always be\n",
                     specified, e.g. swap_checker(bam.dir = '.') if the BAM files are located in the current working directory.\n\n",
bam.chrom.notation   Specifies the chromosome notation used in the BAM files. Possible values: 'chr' (default, to be used when chromosomes\n",
                     are named \"chr1\", \"chr2\", etc.) and 'nochr' (to be used when chromosomes are named \"1\", \"2\", etc.).\n\n",
bam.pattern          This is an optional argument that helps to limit the analysis to a subset of the available BAM files, based on file\n",
                     name patterns. For example, if bam.dir contains both unsorted and sorted BAM files as different intermediate outputs\n",
                     of an NGS processing pipeline, running swap_checker() will result in error because the function accepts only\n",
                     coordinate-sorted and indexed BAM files. However, if file names of the sorted BAM files end in \"sorted.bam\",\n",
                     swap_checker(bam.pattern = '*.sorted.bam$') will analyze only the sorted BAM files without errors. The value of the\n",
                     \"bam.pattern\" argument should be a regular expression like in the previous simple example. Regular expressions are a\n",
                     powerful tool and their use allows very precise selection of BAM files that have to be analyzed. Inexperienced users\n",
                     who would wish to use them in more complicated scenarios might want to consult a quick guide to regular expression\n",
                     syntax, such as https://www.rexegg.com/regex-quickstart.html.\n\n",
snp.list             BED file containing a custom list of SNPs. You can use the SNPs_selected.bed file as a template for your own list. Make\n",
                     sure that you provide genomic coordinates according to the same genome version that was used for the alignment.\n",
                     If not specified, the program will use a list of 28 SNPs in genes frequently mutated in CLL (SNPs_selected.bed) and\n",
                     will assume that your BAM files are aligned to hg19.\n\n",
nthreads             Number of threads to use. Each BAM file will be launched on a separate thread. Works only on Unix and macOS. Default: 4\n\n",
samples              Argument with two possible values: 'single' (default) and 'multiple'. Specifies whether all BAM files originate from\n",
                     separate individuals ('single'), in which case swap_checker() will be able to identify only unexpectedly matching\n",
                     samples, or there are more than one sample per subject ('multiple'), which will allow swap_checker() to also check\n",
                     whether samples labelled as originating from the same individual match as expected. The latter will work only if BAM\n",
                     file names follow an uniform pattern and contain identifiers for the subjects (see argument \"subject.identifier\").\n\n",
subject.identifier   Specifies the part of BAM file names that contains the unique subject identifier. Relevant only when\n",
                     \"samples = 'multiple'\" (see above). This argument accepts one of six preset values or a custom regular expression.\n",
                     The preset values are 'NUMERICID-XXX', 'MIXEDID-XXX', 'XXX-NUMERICID-XXX', 'NUMERICID_XXX', 'MIXEDID_XXX' and\n",
                    'XXX_NUMERICID_XXX'. 'NUMERICID-XXX' (default) tells the swap_checker() function that file names begin with a numeric\n",
                     sequence that is unique for each subject, followed by a hyphen and the rest of the file name. 'MIXEDID-XXX' is similar\n",
                     but allows alphanumeric characters in the subject identifier. Accordingly, 'XXX-NUMERICID-XXX' instructs the\n",
                     swap_checker() function to look for the first digits-only part of the file name that is surrounded by hyphens and to\n",
                     accept it as unique subject identifier. The next three options are analogous to the already described but are suited\n",
                     for the case when file name parts are separated by underscores instead of hyphens. If none of the preset values\n",
                     corresponds to the file naming convention for the BAM files in bam.dir, you can define the subject-identifying part of\n",
                     the file name by using a regular expression, e.g. '^(\\\\d+)-.*$' (this is equivalent to 'NUMERICID-XXX'). Inexperienced\n",
                     users might want to consult a quick guide to regular expression syntax, such as\n",
                     https://www.rexegg.com/regex-quickstart.html.\n\n",
results.dir          Directory in which the results will be saved. By default, the results will be saved in the current working directory.\n\n",
suffix               Optional suffix that can be appended to the main part of result file names. A time stamp suffix is always added\n",
                     automatically and if an optional suffix is specified, it will appear between the main part of the file name and the\n",
                     time stamp. For example, swap_checker(suffix = 'ProjectX') will result in the creation of result files with names like\n",
                     \"Pairwise_concordance_ProjectX_20250531_174936.csv\" and \"SNP_readcounts_ProjectX_20250531_174936.csv\"\n\n",
plot                 Optional argument with two possible values: 'no' (default) or 'yes'. If set to 'yes', swap_checker() will produce a\n",
                     correlation plot. Such a plot can be useful only for relatively small cohorts, otherwise it will be too crowded and\n",
                     unreadable if printed on a single page. If sample names (taken from BAM file names) are too long and take up most of\n",
                     the plotting space, you can shorten them by removing redundant file name parts using the \"remove.name.part\" argument\n",
                     (see below).\n\n",
remove.name.part     Optional argument related to the plotting function (see above). Its value should be a regular expression defining\n",
                     a part of the sample name that should be removed. For example, if all BAM files end in \"_L001_R1_bwa.deduplicated.bam,\n",
                     swap_checker(plot = 'yes', remove.name.part = '_L001.*$') will produce a correlation plot and the displayed sample\n",
                     names will not contain the redundant part from the BAM file names.\n\n\n",
Output files:
Pairwise_concordance_[suffix]_[time stamp].csv       This is the main results file, a table in CSV format with the following columns:\n",
                                                     - \"X_bam\" and \"Y_bam\" contain the sample names (taken from the input BAM files),\n",
                                                     arranged so that each sample in the first column is matched against any of the\n",
                                                     remaining samples in the second column;\n",
                                                     - \"concordant_snps\": number of the concordant SNPs;\n",
                                                     - \"discordant_snps\": number of the discordant SNPs;\n",
                                                     - \"discordant_snps_detailed\": genomic coordinates of the discordant SNPs;\n",
                                                     - \"fract_concordant_snps\": fraction of the concordant SNPs;\n",
                                                     - \"cor_coef\": Pearson's correlation coefficient;\n",
                                                     - \"XY_possibly_paired\": estimation whether samples originate from the same\n",
                                                     individual - \"Yes\" or \"No\". Positive findings are listed first.\n\n",
SNP_readcounts_[suffix]_[time stamp].csv             This is a detailed table with the number of the reads that support the reference or\n",
                                                     the alternative allele at each of the assessed SNP positions in each sample. Columns:\n",
                                                     - \"BAM\": sample name (taken from the respective input BAM file);\n",
                                                     - \"loci\": genomic coordinates of the SNP;\n",
                                                     - \"ref_rc\": read counts supporting the reference allele;\n",
                                                     - \"alt_rc\": read counts supporting the alternative allele;\n",
                                                     - \"vaf\": variant allele fraction of the alternative allele.\n\n",
Expected_concordant_pairs_[suffix]_[time stamp].csv  This is an optional subset of the main table, created if \"samples\" was set to\n",
                                                     'multiple'. It contains only the rows for combinations of samples that are expected\n",
                                                     to be concordant because they should originate from the same individual.\n\n",
Pairwise_concordance_[suffix]_[time stamp].png       This is an optional correlation plot (created if \"plot\" was set to 'yes') based on\n",
                                                     the correlation coefficients from the main results table.\n\n")

#### Preparation of the script
Some variables in the Swap_checker.R file will have to be adjusted to suit your data before you run the script. This can be done most conveniently in RStudio.
* Line 16 – this is the place to specify the list of SNPs that will be used for the analysis. If left unchanged, the built-in list will be used which is specific for our custom NGS panel covering various genes of relevance for chronic lymphocytic leukaemia. Prepare your own list of SNPs following the format of the `SNPs_selected.bed` file and save it in the same directory as the other files of the program. Then modify line 16 of the `Swap_checker.R` file to point to the file with your own list of SNPs. Alternatively, you can directly edit the `SNPs_selected.bed` file and leave line 16 of the `Swap_checker.R` file unmodified.  
* Line 19 – Replace the placeholders `“/path/to/folder/with/BAMs/1”` and `“/path/to/folder/with/BAMs/2”` with the actual paths to the directories containing your bam files, keeping the quotation marks. You can add as many directories as you want inside the brackets, separating them with commas. If all your bam files are in one directory, remove the comma after the path to it, as well as the second placeholder (`“/path/to/folder/with/BAMs/2”`). The `pattern="*.bam$"` part can be used to select only a particular type of BAM files from the data director(y/ies). For example, if you have both "normal" and realigned BAM files in the same directory, you can perform the analysis only on the realigned files if they have a respective uniform file name ending, using a pattern like e.g. `*.realigned.bam$` instead of the default `*.bam$` pattern that will select all bam files.
* Line 21 – Here you have the option to limit the analysis to a subset of the bam files in your director(y/ies), based on their file names, i.e. only bam files that contain a specified text string in their names will be considered. To make use of this possibility, uncomment the line (remove the `# ` at its beginning and replace `CLL` with the actual text that has to be present in the file names. You can use more complicated rules for selection of files if you are familiar with regular expressions.
* Line 24 – Replace `getwd()` with the path (absolute or relative) to the directory where the results should be saved. Otherwise the results will be saved in the directory where the script file is located, potentially overwriting older results without warning. The chosen directory must be pre-existing.
* Line 27 – Set the `ncores` parameter to the number of CPU cores that are available on your system.
* Line 30 – the pairwise concordance table will be saved under the name `Pairwise_concordance.csv`. It is a good idea to modify this generic name by adding to it a descriptor of the analyzed study/samples, thus avoiding the confusion that might ensue when you analyze several groups of samples and get several different results files with the same names in different directories.
* Line 31 – the SNP readcounts for reference and alternative alleles for all samples will be saved in the file `SNP_readcounts.csv`. It is a good idea to modify this generic name by adding to it a descriptor of the analyzed study/samples, thus avoiding the confusion that might ensue when you analyze several groups of samples and get several different results files with the same names in different directories.
* Line 38 – relevant if you want to execute the optional part of the script that checks whether samples that should be concordant are really concordant (e.g. multiple samples from the same individual). This can be tested automatically only if file names contain identifiers for the subjects.  File names should follow a specific and uniform pattern, i.e. the subject identifier should be present in a definable part of each file name. In this way, the subject identifier can be extracted using a regular expression so that samples with the same identifier can be compared to each other. The example regular expression at line 38 `^.*-(\\d+)-.*$` assumes that the file names follow the pattern "something-number-something" where the number uniquely identifies each subject/patient. Modify the regular expression to fit your most probably different naming convention. A good quick guide to regular expressions: https://www.rexegg.com/regex-quickstart.html
* Line 47 – relevant if you want to execute the optional part of the script that checks whether samples that should be concordant are really concordant (see above). The correlation table of all pairs that should be concordant will be saved in the file `Expected_concordant_pairs.csv`. It is a good idea to modify this generic name by adding to it a descriptor of the analyzed study/samples, thus avoiding the confusion that might ensue when you analyze several groups of samples and get several different results files with the same names in different directories.
* Line 52 – relevant if you want to execute the optional part of the script that plots a graphical depiction of the correlations between any two samples of the cohort (this can be useful only for relatively small cohorts, otherwise the plot will be too crowded and unreadable if printed on one page). At this line, you can enter a regular expression between the quotation marks to find and remove repetitive parts of sample names. For example, all of our file names end in `_L001_R1_001.cutadapt.alignMEM.dedupOptPIC.realignGATK.reads.bam`, reflecting all pipeline stages that led to the generation of the final bam file, and entering `_L.*` as a regular expression will remove this whole unnecessary part of the names so that labels can fit on the graph. Modify according to your needs. A good quick guide to regular expressions: https://www.rexegg.com/regex-quickstart.html
* Line 73 – relevant if you want to execute the optional part of the script that plots a graphical depiction of the correlations between any two samples of the cohort (see above). The graph will be saved under the file name `Pairwise_concordance.png`. It is a good idea to modify this generic name by adding to it a descriptor of the analyzed study/samples, thus avoiding the confusion that might ensue when you analyze several groups of samples and get several different results files with the same names in different directories.

#### Execution of the script
After changing the parameters as necessary, start the script. There are different ways to do this but the most convenient way is to do it from within RStudio. Select the part of the script that you want to execute (from the beginning to line 32 if you only want to check for unexpected pairing of samples, from the beginning to line 48 if you also want to check whether expected pairs are paired as they should be, or all of the text if in addition you also want to get a graph of the correlation coefficients) and then click on Run.
Alternatives without RStudio:
* You can navigate in the terminal to the directory where Swap Checker has been installed and then issue the command `Rscript Swap_checker.R`. This will execute the whole script.
* You can start R in a terminal and then enter the command `source("~/Swap_checker/Swap_checker.R")` (if Swap Checker is installed in a different directory on your computer, you will have to modify the path accordingly).
It will take from several minutes (tens of samples) to a few hours (thousands of samples) for the script to execute depending on the speed of your computer, the sizes of the dataset and of the SNP list, and the extent of the analysis (only matching samples or also checking whether expected pairs are really concordant ± producing a graph of the correlation coefficients).
