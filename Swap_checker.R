## Title: Swap Checker - an R script for detection of NGS sample swaps and mislabelling based on (mis)matching SNPs
## File name: Swap_checker.R
## Version: 1.0.0
## Date: 2024-04-19
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: University Hospital Ulm, Germany, 2024
## License: GNU General Public License v3 (GPL-3), https://www.gnu.org/licenses/gpl-3.0.html


### Load necessary packages
library(maftools) # Load the maftools package
source("sample_swaps_modified.R") # Load the modified sampleSwaps function from the sample_swaps_modified.R file

### Specify a bed file with a list of SNPs to be used
snp.list = "SNPs_selected.bed"

### Prepare a list of BAM files to be scanned for possible sample swaps
bams <- bam.list <- list.files(path=c("/path/to/folder/with/BAMs/1", "/path/to/folder/with/BAMs/2"),pattern="*.bam$", full.names = TRUE) # Enter the path(s) to one or several directories with indexed bam files (i.e. the bai files must also be present). The pattern can be used to select only a particular type of BAM files from the data directory. For example, if you have both "normal" and realigned BAM files in the same directory, you can perform the analysis only on the realigned files if they have a respective uniform file name ending, using a pattern like e.g. pattern="*.realigned.bam$"
## Next line is optional if you need to limit your analysis to only a subset of the bam files in your folder(s). You can use regular expressions for precise selection. Uncomment the line to make it work.
# bams <- grep("CLL", bam.list, value = TRUE) # This example will select only BAM files that contain "CLL" in their file names. Uncomment the line to activate it and modify it according to your needs.

### Set directory for saving results
res.dir <- getwd() # Replace "getwd()" with the path (absolute or relative) to the directory where the results should be saved. Otherwise the results will be saved in the directory where the script file is located, potentially overwriting older results without warning. The chosen directory must be pre-existing and its name must be put within quotation marks.

### Check for sample swaps
res = sampleSwaps(bams = bams, ncores = 4, snps = snp.list) # Set the "ncores" parameter according to the CPU cores that are available on your system.

### Save results
write.csv(res$pairwise_comparison, file = file.path(res.dir, "Pairwise_concordance.csv"), row.names = FALSE) # Save the pairwise concordance table
write.csv(res$SNP_readcounts, file = file.path(res.dir, "SNP_readcounts.csv"), row.names = FALSE) # Save SNP readcounts for ref and alt alleles for all samples

### Optional step - Check whether samples that should be concordant are concordant indeed
## Samples that originate from the same subject should be concordant with each other. This can be tested automatically 
## if file names follow a specific pattern and an identifier for the subject is present as a part of each file name.
pc <- res$pairwise_comparison
sample.names <- unique(c(pc$X_bam, pc$Y_bam))
subjects <- gsub("^.*-(\\d+)-.*$", "\\1", sample.names) # This line serves to automatically identify the subjects from the sample names. In this specific case, it is assumed that the sample names follow the pattern "something-number-something" where the number uniquely identifies each subject. Modify the regular expression to fit your most probably different naming convention.
subjects <- unique(subjects)
matched <- data.frame() 
for (i in 1:length(subjects)) {
  a <- pc[grep(subjects[i], pc$X_bam),]
  b <- a[grep(subjects[i], a$Y_bam),]
  matched <- rbind(matched,b)
} 
matched <- matched[order(matched$cor_coef, decreasing = TRUE),]
write.csv(matched, file = file.path(res.dir, "Expected_concordant_pairs.csv"), row.names = FALSE) # Save the pairwise concordance table

### Optional step - graphical depiction of correlation between all pairs of samples (it is useful only for relatively small cohorts)
library(pheatmap) # Load the pheatmap package

regex.pattern <- "" # You can enter a regular expression between the quotation marks to find and remove repetitive parts of sample names. For example, all of our file names end in "_L001_R1_001.cutadapt.alignMEM.dedupOptPIC.realignGATK.reads.bam" and entering "_L.*" here will remove this whole unnecessary part of the names so that labels can fit on the graph. Modify according to your needs.

pc <- res$pairwise_comparison
pc <- pc[,c(1,2,7)]
pc$X_bam <- gsub(regex.pattern, "", pc$X_bam)
pc$Y_bam <- gsub(regex.pattern, "", pc$Y_bam)
pc1 <- cbind(pc$X_bam, pc$Y_bam, pc$cor_coef)
pc2 <- cbind(pc$Y_bam, pc$X_bam, pc$cor_coef)
pc3 <- data.frame(unique(c(pc$X_bam, pc$Y_bam)), unique(c(pc$X_bam, pc$Y_bam)), rep(1, length(unique(c(pc$X_bam, pc$Y_bam)))), stringsAsFactors = FALSE)
pc3 <- cbind(pc3[,1], pc3[,2], pc3[,3])
pc <- rbind(pc1, pc2, pc3)
pc <- as.data.frame(pc, stringsAsFactors = FALSE)
pc$V3 <- as.numeric(pc$V3)
colnames(pc)[1:2] <- c("X_bam", "Y_bam")
cor_table <- reshape(pc, idvar = "X_bam", timevar = "Y_bam", direction = "wide")
colnames(cor_table) <- gsub("(V3.)", "", colnames(cor_table))
sample.names <- cor_table$X_bam
rownames(cor_table) <- sample.names
cor_table <- cor_table[,match(rownames(cor_table), colnames(cor_table))]

pheatmap(cor_table, breaks = seq(0, 1, 0.01), show_rownames = TRUE,
         filename = file.path(res.dir, "Pairwise_concordance.png"), angle_col = 90) # Display correlations as a heatmap and save it to a file
