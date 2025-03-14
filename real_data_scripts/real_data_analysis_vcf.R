# Script for real data analysis for UKBB data, given pre-computed loci and
#   an outcome VCF file
# Used for methods that rely on instrument selection

### IMPORTS

library("optparse")
library(data.table)
library(qvalue)
source("run_methods.R")


### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--sumstat_dir"), type="character", default='./',
              help="Directory of summary statistic files", metavar="character"),
  make_option(c("--ld_dir"), type="character", default='./',
              help="Directory of corresponding LD matrices", metavar="character"),
  make_option(c("--outcome_vcf"), type="character", default='NONE',
              help="VCF for the outcome trait", metavar="character"),
  make_option(c("--sel_method"), type="character", default='pruned',
              help="SNP selection method to use.", metavar="character"),
  make_option(c("--prune_thresh"), type="numeric", default=0.1, metavar="numeric",
              help="Threshold to prune SNPs at for pruned/clumped criterion."),
  make_option(c("--N"), type="numeric", default=100000,
              help="Sample size of exposure data.", metavar="numeric"),
  make_option(c("--methods"), type="character", default='all',
              help="Comma-separated list of methods to run (default: all)", metavar="character"),
  make_option(c("--out"), type="character", default='results.rds',
              help="Output results RDS file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# helper function for LD pruning SNPs
ld_prune_ss = function(sumstats, ld, prune_thresh) {
  if(length(ld) == 1)  return(c(1))  # if only one SNP, no pruning
  # Sort betas in descending order, reshuffle LD matrix
  max_abs_betas <- apply(abs(sumstats$betas), 1, max)
  sorted_ss <- data.frame(max_abs_betas = max_abs_betas, ord = 1:nrow(sumstats$betas))
  sorted_ss <- sorted_ss[order(-sorted_ss$max_abs_betas), ]
  ld <- abs(ld[sorted_ss$ord, sorted_ss$ord])

  # Prune SNPs based on LD threshold in the upper triangle of sorted LD matrix
  to_prune <- rep(FALSE, ncol(ld))
  for (idx in 2:ncol(ld)) {
    if (to_prune[idx]) next  # Skip already pruned SNPs
    vals <- ld[1:(idx-1), idx]  # upper triangle
    if (any(vals[!to_prune[1:(idx-1)]] > prune_thresh)) {
      to_prune[idx] <- TRUE
    }
  }
  keep_idx <- sorted_ss$ord[!to_prune]  # original order
  return(sort(keep_idx))
}


# select SNPs to use as instruments
select_snps = function(indir, ld_dir, criterion, y_ids, prune_thresh) {
  all_ss = list()
  ss_fnames = list.files(path = indir, pattern = '*.rds', full.names = T)
  ld_fnames = list.files(path = ld_dir, pattern = '*.rds', full.names = T)
  for(f_idx in 1:length(ss_fnames)) {
    # read in data and subset SNPs by those available for the outcome
    sumstats = readRDS(ss_fnames[f_idx])
    ld = as.matrix(readRDS(ld_fnames[f_idx]))
    subset_indices = sumstats$pos$ID %in% y_ids
    sub_ss = subset_sumstats(sumstats, subset_indices)
    abs_zscores = abs(as.matrix(sub_ss$Z))
    if(length(abs_zscores) == 0)  next
    ld = ld[subset_indices, subset_indices, drop=F]

    # select SNP based on one of several possible criteria
    if(criterion == 'pruned' || criterion == 'clumped') {  # simply LD prune variants
      idx = ld_prune_ss(sub_ss, ld, prune_thresh)
    } else if(criterion == 'max') {  # get SNP with max abs(Z) score across all traits
      maxloc = which.max(abs_zscores)
      idx = maxloc %% nrow(abs_zscores)
      if(idx == 0)  idx = nrow(abs_zscores)  # if evenly divides, it's in last row
    } else {  # get SNP with max abs(Z) averaged across traits
      abs_z_means = apply(abs_zscores, 1, mean)
      idx = as.numeric(which.max(abs_z_means))
    }
    sel_snp = subset_sumstats(sub_ss, idx)
    all_ss = merge_sumstats(all_ss, sel_snp)  # add onto overall sumstats
  }
  return(all_ss)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # read in outcome data
  outc_dat = read_vcf(opt$outcome_vcf)
  outc_dat = outc_dat %>% select(ID, BETA, SE)  # remove unnecessary fields

  # select SNPs based on max or mean criterial, then add in outcome trait
  sel_ss = select_snps(opt$sumstat_dir, opt$ld_dir, opt$sel_method,
                       outc_dat$ID, opt$prune_thresh)
  outc_dat_sub = outc_dat[match(sel_ss$pos$ID, outc_dat$ID), ]

  # run methods
  methods = unlist(strsplit(tolower(opt$methods), ','))
  raw_res = run_methods(methods, G = c(), X = c(), Y = c(), Z = c(),
              out = opt$out, verbose = TRUE, write.res = FALSE, N = opt$N,
              x_betas = as.matrix(sel_ss$betas), x_stderrs = as.matrix(sel_ss$stderrs),
              y_betas = outc_dat_sub$BETA, y_stderrs = outc_dat_sub$SE)

  # write structured results
  res = list()
  exposure_names = colnames(sel_ss$betas)
  n_expo = ncol(sel_ss$betas)
  for(met in names(raw_res)) {
    res[[met]] = data.frame('name' = exposure_names,
                            'beta' = raw_res[[met]][1:n_expo],
                            'se' = raw_res[[met]][(n_expo+1):(n_expo*2)],
                            'pvalue' = raw_res[[met]][(2*n_expo+1):(n_expo*3)])

    # either compute q-value or set to 1-PIP for methods that give PIPs
    if(grepl('susie', met) || grepl('mrash', met) || grepl('varbvs', met)
       || grepl('brms', met) || grepl('vebboost', met)) {
      res[[met]]$qvalue = 1 - res[[met]]$pvalue
    } else {
      res[[met]]$qvalue = qvalue(res[[met]]$pvalue, lambda=0)$qvalues
    }
  }
  saveRDS(res, opt$out)
}
