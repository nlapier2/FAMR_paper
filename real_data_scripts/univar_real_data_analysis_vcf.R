# Extract summary statistics and select SNPs to use as instruments for univarite MR

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
  make_option(c("--outcome_vcf"), type="character", default='NONE',
              help="Outcome summary statistics vcf file", metavar="character"),
  make_option(c("--regions"), type="character", default='NONE',
              help="File with regions to extract", metavar="character"),
  make_option(c("--old_col_order"), type="logical", action="store_true", default=FALSE,
              help="Use older region file column order (e.g. for INTERVAL)"),
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
  # read in summary statistics
  expo_dat = read_vcf(opt$exposure_vcf)
  outc_dat = read_vcf(opt$outcome_vcf)
  outc_dat = outc_dat %>% select(ID, BETA, SE)  # remove unnecessary fields
  names(outc_dat) = c('ID', 'BETA.Y', 'SE.Y')  # rename BETA/SE to avoid confusion
  # read regions file
  regions = read.table(opt$regions)
  names(regions) = c('CHR', 'START', 'END', 'LENGTH', 'MAXZ')
  if(opt$old_col_order)  names(regions) = c('CHR', 'POS', 'logp', 'START', 'END')

  # loop through and collect sumstats within all regions
  all_reg = list()
  for(idx in 1:nrow(regions)) {
    # extract the SNPs in the current region, append to overall dataframe
    dat_reg = expo_dat[expo_dat$CHR == regions[['CHR']][idx]
                  & expo_dat$POS >= regions[['START']][idx]
                  & expo_dat$POS <= regions[['END']][idx] , ]
    all_reg[[idx]] = dat_reg
  }

  # select snps via different methods
  sel_ss = select_snp_peak_univar(all_reg, outc_dat)

  # add in y betas/stderrs for selected SNPs
  sel_ss = inner_join(sel_ss, outc_dat, by='ID')

  # run MR analysis, write results
  mr_res = mr_ivw(mr_input(bx = as.numeric(sel_ss$BETA), bxse = as.numeric(sel_ss$SE),
                           by = as.numeric(sel_ss$BETA.Y), byse = as.numeric(sel_ss$SE.Y)))
  # saveRDS(mr_res, opt$out)
  mr_vals = as.numeric(c(mr_res$Estimate, mr_res$StdError, mr_res$Pvalue))
  write.table(t(mr_vals), file = opt$out, row.names = FALSE, col.names = FALSE, append = TRUE)
}
