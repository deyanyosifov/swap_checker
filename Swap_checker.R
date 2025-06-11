## Title: Swap Checker - an R function for detection of NGS sample swaps and mislabelling based on (mis)matching SNPs
## File name: Swap_checker.R
## Version: 1.1.0
## Date: 2025-05-28
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: University Hospital Ulm, Germany, 2025
## License: GNU General Public License v3 (GPL-3), https://www.gnu.org/licenses/gpl-3.0.html

swap_checker <- function(..., bam.dir = getwd(), snp.list = cll.snp.list, results.dir = getwd(), subject.identifier = "NUMERICID-XXX", nthreads = 4, bam.pattern="*.bam$", suffix = "", samples = "single", plot = "no", remove.name.part = "", bam.chrom.notation = "chr") {
  
    ## Help message
    args <- list(...)
    if (length(args) == 1 && is.character(args[[1]]) && args[[1]] == "help") {
    cat(
"\nSwap_checker - an R function for detection of NGS sample swaps and mislabelling based on (mis)matching SNPs\n\n",
"Version: 1.1.0\n\n",
"Usage:\n\n",
"  If called from within an R environment:\n",
"    swap_checker(argument1 = 'value1', argument2 = 'value2', ...)\n\n",
"  If called directly from terminal (requires Rscript):\n",
"    Rscript -e \"source('your/installation/path/Swap_checker.R'); swap_checker(argument1 = 'value1', argument2 = 'value2', ...)\"\n\n\n",
"Arguments:\n\n",
"  bam.dir              Directory or list of directories that contain the BAM files to be processed. The BAM files must be coordinate-sorted\n",
"                       and indexed. Both relative and absolute paths are accepted. A list of directories should have the format\n",
"                       c('/path/to/directory/with/BAMs/1','/path/to/directory/with/BAMs/2'). Note that in this case the separate directories\n",
"                       have to be enclosed in single quotation marks but the list as a whole not.\n",
"                       If not specified, swap_checker will look for BAM files in the current working directory of the R environment. This will\n",
"                       not work if swap_checker is called from a terminal via Rscript - in this case the bam.dir argument must always be\n",
"                       specified, e.g. swap_checker(bam.dir = '.') if the BAM files are located in the current working directory.\n\n",
"  bam.chrom.notation   Specifies the chromosome notation used in the BAM files. Possible values: 'chr' (default, to be used when chromosomes\n",
"                       are named \"chr1\", \"chr2\", etc.) and 'nochr' (to be used when chromosomes are named \"1\", \"2\", etc.).\n\n",
"  bam.pattern          This is an optional argument that helps to limit the analysis to a subset of the available BAM files, based on file\n",
"                       name patterns. For example, if bam.dir contains both unsorted and sorted BAM files as different intermediate outputs\n",
"                       of an NGS processing pipeline, running swap_checker() will result in error because the function accepts only\n",
"                       coordinate-sorted and indexed BAM files. However, if file names of the sorted BAM files end in \"sorted.bam\",\n",
"                       swap_checker(bam.pattern = '*.sorted.bam$') will analyze only the sorted BAM files without errors. The value of the\n",
"                       \"bam.pattern\" argument should be a regular expression like in the previous simple example. Regular expressions are a\n",
"                       powerful tool and their use allows very precise selection of BAM files that have to be analyzed. Inexperienced users\n",
"                       who would wish to use them in more complicated scenarios might want to consult a quick guide to regular expression\n",
"                       syntax, such as https://www.rexegg.com/regex-quickstart.html.\n\n",
"  snp.list             BED file containing a custom list of SNPs. You can use the SNPs_selected.bed file as a template for your own list. Make\n",
"                       sure that you provide genomic coordinates according to the same genome version that was used for the alignment.\n",
"                       If not specified, the program will use a list of 28 SNPs in genes frequently mutated in CLL (SNPs_selected.bed) and\n",
"                       will assume that your BAM files are aligned to hg19.\n\n",
"  nthreads             Number of threads to use. Each BAM file will be launched on a separate thread. Works only on Unix and macOS. Default: 4\n\n",
"  samples              Argument with two possible values: 'single' (default) and 'multiple'. Specifies whether all BAM files originate from\n",
"                       separate individuals ('single'), in which case swap_checker() will be able to identify only unexpectedly matching\n",
"                       samples, or there are more than one sample per subject ('multiple'), which will allow swap_checker() to also check\n",
"                       whether samples labelled as originating from the same individual match as expected. The latter will work only if BAM\n",
"                       file names follow an uniform pattern and contain identifiers for the subjects (see argument \"subject.identifier\").\n\n",
"  subject.identifier   Specifies the part of BAM file names that contains the unique subject identifier. Relevant only when\n",
"                       \"samples = 'multiple'\" (see above). This argument accepts one of six preset values or a custom regular expression.\n",
"                       The preset values are 'NUMERICID-XXX', 'MIXEDID-XXX', 'XXX-NUMERICID-XXX', 'NUMERICID_XXX', 'MIXEDID_XXX' and\n",
"                       'XXX_NUMERICID_XXX'. 'NUMERICID-XXX' (default) tells the swap_checker() function that file names begin with a numeric\n",
"                       sequence that is unique for each subject, followed by a hyphen and the rest of the file name. 'MIXEDID-XXX' is similar\n",
"                       but allows alphanumeric characters in the subject identifier. Accordingly, 'XXX-NUMERICID-XXX' instructs the\n",
"                       swap_checker() function to look for the first digits-only part of the file name that is surrounded by hyphens and to\n",
"                       accept it as unique subject identifier. The next three options are analogous to the already described but are suited\n",
"                       for the case when file name parts are separated by underscores instead of hyphens. If none of the preset values\n",
"                       corresponds to the file naming convention for the BAM files in bam.dir, you can define the subject-identifying part of\n",
"                       the file name by using a regular expression, e.g. '^(\\\\d+)-.*$' (this is equivalent to 'NUMERICID-XXX'). Inexperienced\n",
"                       users might want to consult a quick guide to regular expression syntax, such as\n",
"                       https://www.rexegg.com/regex-quickstart.html.\n\n",
"  results.dir          Directory in which the results will be saved. By default, the results will be saved in the current working directory.\n\n",
"  suffix               Optional suffix that can be appended to the main part of result file names. A time stamp suffix is always added\n",
"                       automatically and if an optional suffix is specified, it will appear between the main part of the file name and the\n",
"                       time stamp. For example, swap_checker(suffix = 'ProjectX') will result in the creation of result files with names like\n",
"                       \"Pairwise_concordance_ProjectX_20250531_174936.csv\" and \"SNP_readcounts_ProjectX_20250531_174936.csv\"\n\n",
"  plot                 Optional argument with two possible values: 'no' (default) or 'yes'. If set to 'yes', swap_checker() will produce a\n",
"                       correlation plot. Such a plot can be useful only for relatively small cohorts, otherwise it will be too crowded and\n",
"                       unreadable if printed on a single page. If sample names (taken from BAM file names) are too long and take up most of\n",
"                       the plotting space, you can shorten them by removing redundant file name parts using the \"remove.name.part\" argument\n",
"                       (see below).\n\n",
"  remove.name.part     Optional argument related to the plotting function (see above). Its value should be a regular expression defining\n",
"                       a part of the sample name that should be removed. For example, if all BAM files end in \"_L001_R1_bwa.deduplicated.bam,\n",
"                       swap_checker(plot = 'yes', remove.name.part = '_L001.*$') will produce a correlation plot and the displayed sample\n",
"                       names will not contain the redundant part from the BAM file names.\n\n\n",
"Output files:\n\n",
"  Pairwise_concordance_[suffix]_[time stamp].csv       This is the main results file, a table in CSV format with the following columns:\n",
"                                                       - \"X_bam\" and \"Y_bam\" contain the sample names (taken from the input BAM files),\n",
"                                                       arranged so that each sample in the first column is matched against any of the\n",
"                                                       remaining samples in the second column;\n",
"                                                       - \"concordant_snps\": number of the concordant SNPs;\n",
"                                                       - \"discordant_snps\": number of the discordant SNPs;\n",
"                                                       - \"discordant_snps_detailed\": genomic coordinates of the discordant SNPs;\n",
"                                                       - \"fract_concordant_snps\": fraction of the concordant SNPs;\n",
"                                                       - \"cor_coef\": Pearson's correlation coefficient;\n",
"                                                       - \"XY_possibly_paired\": estimation whether samples originate from the same\n",
"                                                       individual - \"Yes\" or \"No\". Positive findings are listed first.\n\n",
"  SNP_readcounts_[suffix]_[time stamp].csv             This is a detailed table with the number of the reads that support the reference or\n",
"                                                       the alternative allele at each of the assessed SNP positions in each sample. Columns:\n",
"                                                       - \"BAM\": sample name (taken from the respective input BAM file);\n",
"                                                       - \"loci\": genomic coordinates of the SNP;\n",
"                                                       - \"ref_rc\": read counts supporting the reference allele;\n",
"                                                       - \"alt_rc\": read counts supporting the alternative allele;\n",
"                                                       - \"vaf\": variant allele fraction of the alternative allele.\n\n",
"  Expected_concordant_pairs_[suffix]_[time stamp].csv  This is an optional subset of the main table, created if \"samples\" was set to\n",
"                                                       'multiple'. It contains only the rows for combinations of samples that are expected\n",
"                                                       to be concordant because they should originate from the same individual.\n\n",
"  Pairwise_concordance_[suffix]_[time stamp].png       This is an optional correlation plot (created if \"plot\" was set to 'yes') based on\n",
"                                                       the correlation coefficients from the main results table.\n\n")
    return(invisible(NULL))
  } else {
    
    valid_samples <- c("single", "multiple")
    
    if (!samples %in% valid_samples) {
      stop(
        sprintf(
          "Invalid value for argument 'samples': '%s'.\nValid options are: %s.",
          samples,
          paste(shQuote(valid_samples), collapse = " or ")
        ),
        call. = FALSE
      )
    }
    
    valid_plot <- c("no", "yes")
    
    if (!plot %in% valid_plot) {
      stop(
        sprintf(
          "Invalid value for argument 'plot': '%s'.\nValid options are: %s.",
          plot,
          paste(shQuote(valid_plot), collapse = " or ")
        ),
        call. = FALSE
      )
    }
    
    valid_bam.chrom.notation <- c("chr", "nochr")
    
    if (!bam.chrom.notation %in% valid_bam.chrom.notation) {
      stop(
        sprintf(
          "Invalid value for argument 'bam.chrom.notation': '%s'.\nValid options are: %s.",
          bam.chrom.notation,
          paste(shQuote(valid_bam.chrom.notation), collapse = " or ")
        ),
        call. = FALSE
      )
    }
    
    
    ## Load necessary packages
    library(maftools) 
    
    ## Give a warning if the user did not specify a custom list of SNPs
    if (is.data.frame(snp.list)) {
      cat("\033[31mSwap_checker will use the list of SNPs provided with the program and will assume that your BAM files are aligned to hg19. If this is not what you want, use the snp.list argument to import your own list of SNPs from a BED file and make sure that the coordinates used in the BED file are according to the genome version that was used for the alignment.\033[0m\n", fill = TRUE)
    }
    
    ## Prepare a list of BAM files to be scanned for possible sample swaps
    bams <- bam.list <- list.files(path=bam.dir,pattern=bam.pattern, full.names = TRUE)

    ### Check for sample swaps
    res = sampleSwaps(bams = bams, ncores = nthreads, snps = snp.list, bam.chrom.notation = bam.chrom.notation)
    
    ### Save results
    now <- Sys.time()
    write.csv(res$pairwise_comparison, file = file.path(results.dir, paste0("Pairwise_concordance_", paste0(suffix, "_"), format(now, "%Y%m%d_%H%M%S"), ".csv")), row.names = FALSE) # Save the pairwise concordance table
    write.csv(res$SNP_readcounts, file = file.path(results.dir, paste0("SNP_readcounts_", paste0(suffix, "_"), format(now, "%Y%m%d_%H%M%S"), ".csv")), row.names = FALSE) # Save SNP readcounts for ref and alt alleles for all samples
    
    ### Delete temporary files
    unlink(gsub("\\.bam$", "_nucleotide_counts.tsv",bams))
    
    ### Optional step - Check whether samples that should be concordant are concordant indeed
    ## Samples that originate from the same subject should be concordant with each other. This can be tested automatically 
    ## if file names follow a specific pattern and an identifier for the subject is present as a part of each file name.
    if (samples == "multiple") {
      if (subject.identifier == "XXX-NUMERICID-XXX") {
        sub.id = "^.*-(\\d+)-.*$"
      } else if (subject.identifier == "NUMERICID-XXX") {
        sub.id = "^(\\d+)-.*$"
      } else if (subject.identifier == "MIXEDID-XXX") {
        sub.id = "^(.*?)-.*$"
      } else if (subject.identifier == "XXX_NUMERICID_XXX") {
        sub.id = "^.*_(\\d+)_.*$"
      } else if (subject.identifier == "NUMERICID_XXX") {
        sub.id = "^(\\d+)_.*$"
      } else if (subject.identifier == "MIXEDID_XXX") {
        sub.id = "^(.*?)_.*$"
      } else {
        sub.id = subject.identifier
      }
      pc <- res$pairwise_comparison
      sample.names <- unique(c(pc$X_bam, pc$Y_bam))
      subjects <- gsub(sub.id, "\\1", sample.names)
      subjects <- unique(subjects)
      matched <- data.frame() 
      for (i in 1:length(subjects)) {
        a <- pc[grep(gsub("\\(.*\\)", subjects[i], sub.id), pc$X_bam),]
        b <- a[grep(gsub("\\(.*\\)", subjects[i], sub.id), a$Y_bam),]
        matched <- rbind(matched,b)
      } 
      matched <- matched[order(matched$cor_coef, decreasing = TRUE),]
      write.csv(matched, file = file.path(results.dir, paste0("Expected_concordant_pairs_", paste0(suffix, "_"), format(now, "%Y%m%d_%H%M%S"), ".csv")), row.names = FALSE) # Save the pairwise concordance table
    }
    
    ### Optional step - graphical depiction of correlation between all pairs of samples (it is useful only for relatively small cohorts)
    if (plot == "yes") {
      library(pheatmap) # Load the pheatmap package
      
      pc <- res$pairwise_comparison
      pc <- pc[,c(1,2,7)]
      pc$X_bam <- gsub(remove.name.part, "", pc$X_bam)
      pc$Y_bam <- gsub(remove.name.part, "", pc$Y_bam)
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
               filename = file.path(results.dir, paste0("Pairwise_concordance_", paste0(suffix, "_"), format(now, "%Y%m%d_%H%M%S"), ".png")), angle_col = 90) # Display correlations as a heatmap and save it to a file
    }
  }
}
### The sampleSwaps() function was modified from the original one published in the maftools package by Anand Mayakonda. The original license is reproduced below followed by the modified function.
## MIT License
## Copyright (c) 2018 Anand Mayakonda
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.

sampleSwaps = function(bams = NULL, bam.chrom.notation = bam.chrom.notation, min_depth = 30, ncores = 4 , snps = "", ...){
  
  if(length(bams) < 2){
    stop("Needs 2 or more BAM files!")
  }
  
  if (is.character(snps)) {
    snps = data.table::fread(input = snps, sep = "\t")
  }
  
  if(bam.chrom.notation == "nochr") {
    snps$chr = gsub(pattern = "chr", replacement = '', x = snps$chr, fixed = TRUE)
  }
  
  rc = bamreadcounts(bam = bams, loci = snps, nthreads = ncores, ...)
  
  snps[, id := paste0(chr, ":", start)]
  
  cat("Summarizing allele frequency table..\n")
  rc_af = lapply(rc, function(x){
    xrc = merge(x, snps[,.(id, ref, alt)], by.x = 'loci', by.y = 'id')
    
    xrc_acounts = apply(xrc, 1, function(x){
      
      ref_allele = x['ref']
      ref_rc = 0
      
      ref_rc = switch (ref_allele,
                       A = x['A'],
                       T = x['T'],
                       G = x['G'],
                       C = x['C'])
      
      alt_allele = x['alt']
      alt_rc = 0
      
      alt_rc = switch (alt_allele,
                       A = x['A'],
                       T = x['T'],
                       G = x['G'],
                       C = x['C'])
      
      vaf_tbl = data.table::data.table(ref_rc = as.numeric(ref_rc), alt_rc = as.numeric(alt_rc), loci = x["loci"])
      vaf_tbl[, vaf := alt_rc/(ref_rc+alt_rc)]
      vaf_tbl
    })
    xrc = merge(xrc, data.table::rbindlist(l = xrc_acounts), by = "loci")
    xrc[,.(loci, ref_rc, alt_rc, vaf)]
  })
  
  rc_bind = data.table::rbindlist(l = rc_af, idcol = "sample")
  rc_bind = rc_bind[!is.nan(vaf)]
  if(nrow(rc_bind) == 0){
    stop("Zero SNPs to analyze!")
  }
  rc_bind = rc_bind[, total := ref_rc + alt_rc][total > min_depth]
  
  rc_df = data.table::dcast(data = rc_bind, loci ~ sample, value.var = 'vaf', fill = NA)
  data.table::setDF(x = rc_df, rownames = rc_df$loci)
  rc_df$loci = NULL
  rc_df = rc_df[complete.cases(rc_df),]
  
  cat("Performing pairwise comparison..\n")
  sample_matches = parallel::mclapply(seq_along(rc_af), function(idx){
    
    x = rc_af[[idx]]
    if(idx == length(rc_af)){
      return(NULL)
    }else{
      rest_samps = seq_along(rc_af)[(idx+1):(length(rc_af))]
    }
    #rest_samps = setdiff(seq_along(rc_af), idx)
    cat("Comparing", names(rc_af)[idx], "against rest..\n")
    
    x_compare = lapply(rest_samps,function(rest_idx){
      y = rc_af[[rest_idx]]
      
      concordant_snps = lapply(seq_along(1:nrow(y)), function(row_idx){
        fisher.test(x = matrix(c(x[row_idx, ref_rc], x[row_idx, alt_rc], y[row_idx, ref_rc], y[row_idx, alt_rc]), nrow = 2, ncol = 2))$p.value
      })
      concordant_snps_detailed = unlist(lapply(concordant_snps, function(x) ifelse(x < 0.01, no = "concordant", yes = "discordant"))) # vector of the status of each SNP
      discordant <- paste(x$loci[which(concordant_snps_detailed == "discordant")], collapse = ",") # Positions of discordant SNPs
      concordant_snps = table(concordant_snps_detailed)
      if (length(concordant_snps) == 1) {
        concordant_snps[2] <- 0
        names(concordant_snps) <- c("concordant", "discordant")
      }
      cor_coef = cor.test(x$vaf, y$vaf)$estimate
      
      data.table::data.table(
        X_bam = names(rc_af)[idx],
        Y_bam = names(rc_af)[rest_idx],
        concordant_snps = concordant_snps[[1]],
        discordant_snps = concordant_snps[[2]],
        discordant_snps_detailed = discordant, # additional column with positions of discordant SNPs
        fract_concordant_snps = prop.table(concordant_snps)[[1]],
        cor_coef = cor_coef
      )
    })
    
    x_compare = data.table::rbindlist(l = x_compare)
    #>70% concordant SNPs and a >0.9 cor.coef suggest high possibility of matched samples
    x_compare[, XY_possibly_paired := ifelse(test = fract_concordant_snps >= 0.70 & cor_coef >= 0.90, yes = "Yes", no = "No")]
    x_compare
  }, mc.cores = ncores)
  
  sample_matches = sample_matches[1:(length(sample_matches)-1)]
  
  pos_matches = lapply(sample_matches, function(sample_pair){
    if(nrow(sample_pair) > 0){
      unique(unlist(sample_pair[XY_possibly_paired == "Yes", .(X_bam, Y_bam)], use.names = FALSE))
    }
  })
  pos_matches = pos_matches[which(lapply(pos_matches, function(x) length(x) > 0) == TRUE)]
  
  cat("Done!\n")
  
  list(
    AF_table = invisible(rc_df),
    SNP_readcounts = data.table::rbindlist(
      l = rc_af,
      idcol = "BAM",
      use.names = TRUE,
      fill = TRUE
    ),
    pairwise_comparison = data.table::rbindlist(sample_matches)[order(XY_possibly_paired, decreasing = TRUE)],
    BAM_matches = pos_matches
  )
}

### List of SNPs covered by our CLL panel (hg19 coordinates)
cll.snp.list <- data.table::data.table(
    chr = c("chr1", "chr1", "chr2", "chr4", "chr6", "chr7", "chr9", 
            "chr11", "chr11", "chr11", "chr12", "chr16", "chr16", "chr16", 
            "chr16", "chr16", "chr16", "chr16", "chr16", "chr17", "chr18", 
            "chr19", "chr22", "chr22", "chr22", "chr22", "chr22", "chrX"),
    start = c(115256669L, 150552042L, 111907691L, 153332380L, 44228162L, 124469267L, 139391636L, 64039175L, 102207851L, 
              108225483L, 25378456L, 81819768L, 81888152L, 81902797L, 81914583L, 81927217L, 81941319L, 81954789L, 81971403L, 
              7579801L, 60985879L, 49459104L, 18226521L, 23054916L, 23055539L, 23241705L, 23247082L, 100611285L), 
    ref = c("G", "G", "T", "T", "G", "A", "G", "G", "G", "C", "T", "T", "A", "A", "C", "G", "C", "C", "T", "G", "T", 
            "A", "C", "G", "G", "G", "C", "T"), 
    alt = c("A", "A", "C", "G", "A", "T", "A", "A", "A", "T", "C", "C", "G", "G", "T", "C", "T", "G", "C", "C", "C", 
            "G", "T", "T", "T", "A", "G", "C"), 
    rsid = c("rs969273", "rs16837903", "rs724710", "rs12644477", "rs2233437", "rs10263573", "rs2229974", "rs2286615", "rs1055088", 
             "rs664982", "rs76433096", "rs1143685", "rs1143686", "rs12445580", "rs11865395", "rs4889430", "rs1143689", "rs12446127", 
             "rs1071644", "rs1642785", "rs1801018", "rs1805419", "rs738094", "rs5759408", "rs1985791", "rs62218846", "rs2009433", 
             "rs3747288")
    ) 
