# simple helper script to run FAMR package on real data
library(FAMR)
library(optparse)

option_list = list(
  make_option(c("-i", "--in_dir"), type="character", default='./',
              help="Directory of summary statistics files to use", metavar="character"),
  make_option(c("--ld_dir"), type="character", default='./',
              help="Directory of LD matrices for loci", metavar="character"),
  make_option(c("--y_gwas_file"), type="character", default='NONE',
              help="File with GWAS results for outcome phenotype.", metavar="character"),
  make_option(c("--fa_method"), type="character", default='gfa',
              help="Factor analysis method to use", metavar="character"),
  make_option(c("-N", "--num_samples"), type="numeric", default=10000,
              help="Sample size of exposure data.", metavar="numeric"),
  make_option(c("--susieL"), type="numeric", default=30,
              help="SuSiE 'L' setting (max number of effect variables)", metavar="numeric"),
  make_option(c("--n_iter"), type="numeric", default=30,
              help="Number of iterations to run FAMR EM loop", metavar="numeric"),
  make_option(c("--idcol"), type="numeric", default=1,
              help="Column number with SNP IDs in y_gwas_file.", metavar="numeric"),
  make_option(c("--betacol"), type="numeric", default=2,
              help="Column number with SNP betas in y_gwas_file.", metavar="numeric"),
  make_option(c("--secol"), type="numeric", default=3,
              help="Column number with SNP std errors in y_gwas_file.", metavar="numeric"),
  make_option(c("--no_header"), type="logical", action="store_true", default=FALSE,
              help="Use if sumstat file has no header.", metavar="logical"),
  make_option(c("--prune_thresh"), type="numeric", default=0.5, metavar="numeric",
              help="Threshold for LD pruning (default = 0.5)"),
  make_option(c("--fa_prune_thresh"), type="numeric", default=0.1, metavar="numeric",
              help="Threshold for LD pruning for factor analysis (default = 0.1)"),
  make_option(c("--final_prune_thresh"), type="numeric", default=0.1, metavar="numeric",
              help="Threshold for LD pruning for final regressions (default = 0.1)"),
  make_option(c("--annihilate_factors"), type="logical", action="store_true", default=FALSE,
              help="Project factors out of X & Y instead of including as covariates"),
  make_option(c("--out"), type="character", default='res_',
              help="Prefix of output and temporary files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

res = famr_susie(opt$in_dir, opt$ld_dir, opt$y_gwas_file, fa_method=opt$fa_method, 
                 num_samples=opt$num_samples, n_iter=opt$n_iter, susieL=opt$susieL, 
                 prune_thresh=opt$prune_thresh, fa_prune_thresh=opt$fa_prune_thresh, 
                 final_prune_thresh=opt$final_prune_thresh, 
                 annihilate=opt$annihilate_factors, idcol=opt$idcol, 
                 betacol=opt$betacol, secol=opt$secol, header=!opt$no_header) 
saveRDS(res, opt$out)
