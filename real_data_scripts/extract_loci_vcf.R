# Extract loci from a regions file from a given exposure VCF

### IMPORTS

library(optparse)
library(dplyr)
library(data.table)
library(MendelianRandomization)
source("utils.R")


### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--exposure_vcf"), type="character", default='NONE',
              help="Exposure summary statistics vcf file", metavar="character"),
  make_option(c("--merged_regions"), type="character", default='NONE',
              help="File with regions to extract", metavar="character"),
  make_option(c("--out"), type="character", default='selected_snps.rds',
              help="Output file to write", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# select SNP with 'peak' Zscore, either max or mean across all traits
select_snp_peak_univar = function(all_reg, y_gwas) {
  all_ss = list()
  for(sumstats in all_reg) {
    # exclude SNPs not in y_gwas
    subset_indices = sumstats$ID %in% y_gwas$ID
    subset_ids = sumstats$ID[subset_indices]
    zscores = sumstats$Z[subset_indices]
    if(length(zscores) == 0)  next

    # select SNP with max absolute Z-score, append to overall set of sumstats
    idx = which.max(abs(zscores))
    sel_snp = sumstats[sumstats$ID == subset_ids[idx], ]
    all_ss = rbind(all_ss, sel_snp)
  }
  return(all_ss)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # read in VCF and merged regions file
  dat = read_vcf(opt$exposure_vcf)
  merged_regions = read.table(opt$merged_regions)
  names(merged_regions) = c('CHR', 'START', 'END', 'LENGTH', 'MAXZ')

  # Create a combined filter condition for all regions
  combined_condition <- merged_regions %>%
    mutate(CHR = as.numeric(CHR), START = as.numeric(START), END = as.numeric(END)) %>%
    rowwise() %>%
    mutate(condition = list(dat$CHR == CHR & dat$POS >= START & dat$POS <= END)) %>%
    ungroup() %>%
    pull(condition) %>%
    Reduce(`|`, .)

  # filter and save results
  dat = dat %>% filter(combined_condition)
  saveRDS(dat, file = opt$out)
}
