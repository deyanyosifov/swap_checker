## Title: Modified sampleSwaps function
## File name: sample_swaps_modified.R
## Version: 1.0.0
## Date: 2024-04-19
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: Anand Mayakonda,  2018 and Deyan Yordanov Yosifov, 2024
## License: GNU General Public License v3 (GPL-3), https://www.gnu.org/licenses/gpl-3.0.html
## The original sampleSwaps function is part of the maftools package by Anand Mayakonda, distributed under the MIT license. The original copyright notice is provided below.

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

#' Identify sample swaps and similarities - modified version to accommodate for the case when there are no discordant SNPs
#' @description Given a list BAM files, the function genotypes known SNPs and identifies potentially related samples. For the source of SNPs, see reference
#' @param bams Input bam files. Required.
#' @param build reference genome build. Default "hg19". Can be hg19 or hg38
#' @param prefix Prefix to add or remove from contig names in SNP file. If BAM files are aligned GRCh37/38 genome, use prefix `chr` to `add`
#' @param add If prefix is used, default is to add prefix to contig names in SNP file. If FALSE prefix will be removed from contig names.
#' @param min_depth Minimum read depth for a SNP to be considered. Default 30.
#' @param ncores Default 4. Each BAM file will be launched on a separate thread. Works only on Unix and macOS.
#' @param ... Additional arguments passed to \code{\link{bamreadcounts}}
#' @return a list with results summarized
#' @export
#' @references Westphal, M., Frankhouser, D., Sonzone, C. et al. SMaSH: Sample matching using SNPs in humans. BMC Genomics 20, 1001 (2019). https://doi.org/10.1186/s12864-019-6332-7
#'
sampleSwaps = function(bams = NULL, build = "hg19", prefix = NULL, add = TRUE, min_depth = 30, ncores = 4 , snps = "", ...){
  
  if(length(bams) < 2){
    stop("Needs 2 or more BAM files!")
  }
  
  #build = match.arg(arg = build, choices = c("hg19", "hg38"))
  
  #if(build == "hg19"){
    #snps = system.file("extdata", "hg19_smash_snps.tsv.gz", package = "maftools")
    #snps = "inst/extdata/hg19_smash_snps.tsv.gz"
  #}else{
    #snps = "inst/extdata/hg38_smash_snps.tsv.gz"
    #snps = system.file("extdata", "hg38_smash_snps.tsv.gz", package = "maftools")
  #}
  snps = data.table::fread(input = snps, sep = "\t")
  
  if(!is.null(prefix)){
    if(add){
      snps$chr = paste(prefix, snps$chr, sep = '')
    }else{
      snps$chr = gsub(pattern = prefix, replacement = '', x = snps$chr, fixed = TRUE)
    }
  }
  
  rc = bamreadcounts(bam = bams, loci = snps, nthreads = ncores, ...)
  
  snps[, id := paste0(chr, ":", start)]
  
  cat("Summarizing allele frequncy table..\n")
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
  
  pos_mathces = lapply(sample_matches, function(sample_pair){
    if(nrow(sample_pair) > 0){
      unique(unlist(sample_pair[XY_possibly_paired == "Yes", .(X_bam, Y_bam)], use.names = FALSE))
    }
  })
  pos_mathces = pos_mathces[which(lapply(pos_mathces, function(x) length(x) > 0) == TRUE)]
  
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
    BAM_matches = pos_mathces
  )
}