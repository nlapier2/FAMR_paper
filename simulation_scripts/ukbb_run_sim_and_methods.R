# samples UK Biobank genotypes and causal SNPs, simulates phenotypes,
#   runs MVMR methods on those, and writes the results

### IMPORTS

source("run_methods.R")
#source("famr.R")
source("sim_mvmr_sparse.R")
library(FAMR)
library(pgenlibr)


### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("-i", "--in_dir"), type="character", default='./',
              help="Directory of plink loci to use as genotypes", metavar="character"),
  make_option(c("-o", "--out_prefix"), type="character", default='res_',
              help="Prefix of output and temporary files", metavar="character"),
  make_option(c("--ld_dir"), type="character", default='./',
              help="Directory of LD matrices for loci", metavar="character"),
  make_option(c("--two_sample"), type="logical", action='store_true', default=FALSE,
              help="Use this flag to active two-sample mode", metavar="logical"),
  make_option(c("--criterion"), type="character", default='mvsusie',
              help="Criterion for SNP selection (mvsusie/max/mean)", metavar="character"),
  make_option(c("--num_sims"), type="numeric", default=1,
              help="Number of simulation replicates to run", metavar="numeric"),
  make_option(c("-N", "--num_samples"), type="numeric", default=20000,
              help="Number of samples to subset from loci", metavar="numeric"),
  make_option(c("--num_loci"), type="numeric", default=200,
              help="Number of loci to implant causal SNPs in", metavar="numeric"),
  make_option(c("--snps_per_locus"), type="numeric", default=1000,
              help="SNPs to retain per locus (incl. non-causal)", metavar="numeric"),
  make_option(c("--max_causal"), type="numeric", default=3,
              help="Maximum num of causal SNPs per locus", metavar="numeric"),
  make_option(c("-K", "--num_exposures"), type="numeric", default=1,
              help="Number of exposures / risk factors to simulate", metavar="numeric"),
  make_option(c("-J", "--num_confounders"), type="numeric", default=1,
              help="Number of correlated pleiotropic confounders to simulate", metavar="numeric"),
  make_option(c("-H", "--num_uncor"), type="numeric", default=1,
              help="Number of uncorrelated pleiotropic confounders to simulate", metavar="numeric"),
  make_option(c("-L", "--num_med"), type="numeric", default=1,
              help="Number of non-confounding G-X mediators to simulate", metavar="numeric"),
  make_option(c("--beta"), type="character", default="0.2", metavar="numeric",
              help="Comma separated list of exposure-outcome effects, e.g '0,0.5,1'"),
  make_option(c("--gamma_gx"), type="numeric", default=0.5,
              help="Size of effect of genotypes on exposures", metavar="numeric"),
  make_option(c("--gamma_gy"), type="numeric", default=0.2,
              help="Size of effect of genotypes on outcome", metavar="numeric"),
  make_option(c("--gamma_gz"), type="numeric", default=0.2,
              help="Size of effect of genotypes on confounders correlated with X", metavar="numeric"),
  make_option(c("--gamma_gu"), type="numeric", default=0.0,
              help="Size of effect of genotypes on confounders uncorrelated with X", metavar="numeric"),
  make_option(c("--gamma_gw"), type="numeric", default=0.0,
              help="Size of effect of genotypes on non-confounding G-X mediators", metavar="numeric"),
  make_option(c("--gamma_zx"), type="numeric", default=0.2,
              help="Size of effect of confounders on exposures", metavar="numeric"),
  make_option(c("--gamma_zy"), type="numeric", default=0.2,
              help="Size of effect of correlated pleiotropic confounders on outcome", metavar="numeric"),
  make_option(c("--gamma_uy"), type="numeric", default=0.0,
              help="Size of effect of uncorrelated pleiotropic confounders on outcome", metavar="numeric"),
  make_option(c("--gamma_wx"), type="numeric", default=0.0,
              help="Size of effect of non-confounding G-X mediators on X", metavar="numeric"),
  make_option(c("--gamma_gx2"), type="numeric", default=0.0,
              help="Size of effect of G on Xs downstream of Z", metavar="numeric"),
  make_option(c("--gamma_yx2"), type="numeric", default=0.0,
              help="Size of effect of reverse-causal effect of Y on downstream Xs", metavar="numeric"),
  make_option(c("--gamma_x1z"), type="numeric", default=0.0,
              help="Size of effect of Xs upstream of Z on Z", metavar="numeric"),
  make_option(c("--phi_gx"), type="numeric", default=0.5,
              help="Density of effect of genotypes on exposures", metavar="numeric"),
  make_option(c("--phi_gy"), type="numeric", default=0.5,
              help="Density of effect of genotypes on outcome", metavar="numeric"),
  make_option(c("--phi_gz"), type="numeric", default=0.5,
              help="Density of effect of genotypes on confounders correlated with X", metavar="numeric"),
  make_option(c("--phi_gu"), type="numeric", default=0.5,
              help="Density of effect of genotypes on confounders uncorrelated with X", metavar="numeric"),
  make_option(c("--phi_gw"), type="numeric", default=0.5,
              help="Density of effect of genotypes on non-confounding G-X mediators", metavar="numeric"),
  make_option(c("--phi_zx"), type="numeric", default=0.5,
              help="Density of effect of confounders on exposures", metavar="numeric"),
  make_option(c("--phi_wx"), type="numeric", default=0.5,
              help="Density of effect of non-confounding G-X mediators on X", metavar="numeric"),
  make_option(c("--phi_gx2"), type="numeric", default=0.5,
              help="Density of effect of G on Xs downstream of Z", metavar="numeric"),
  make_option(c("--phi_yx2"), type="numeric", default=0.5,
              help="Density of effect of reverse-causal effect of Y on downstream Xs", metavar="numeric"),
  make_option(c("--phi_x1z"), type="numeric", default=0.5,
              help="Density of effect of Xs upstream of Z on Z", metavar="numeric"),
  make_option(c("--mu"), type="numeric", default=0,
              help="Set nonzero mean effect sizes (directional pleiotropy)", metavar="numeric"),
  make_option(c("--methods"), type="character", default='all',
              help="Comma-separated list of methods to run (default: all)", metavar="character"),
  make_option(c("--susieL"), type="numeric", default=10,
              help="SuSiE 'L' setting (max number of effect variables", metavar="numeric"),
  make_option(c("--upstream_pct"), type="numeric", default=0.0,
              help="Fraction of Xs upstream of factors", metavar="numeric"),
  make_option(c("--reverse_cause"), action="store_true", default=FALSE,
              help="Whether to have Y 'reverse cause' some Xs", metavar="logical"),
  make_option(c("--pip_thresh"), type="numeric", default=0.1,
              help="PIP threshold for SNP selection using mvSusie", metavar="numeric"),
  make_option(c("--famr_n_iter"), type="numeric", default=30,
              help="Number of iterations to run FAMR EM loop", metavar="numeric"),
  make_option(c("--famr_prune"), type="numeric", default=1.0,
              help="LD pruning threshold for FAMR", metavar="numeric"),
  make_option(c("--fa_prune_thresh"), type="numeric", default=1.0, metavar="numeric",
              help="Threshold for LD pruning for factor analysis (default = 1.0 [none])"),
  make_option(c("--final_prune_thresh"), type="numeric", default=1.0, metavar="numeric",
              help="Threshold for LD pruning for final regressions (default = 1.0 [none])"),
  make_option(c("--sel"), type="numeric", metavar="numeric", default=0,
              help="Perform SNP selection based only on specified exposure (TEST ONLY)"),
  make_option(c("--skip_selection"), type="logical", action="store_true", default=FALSE,
              help="Skip SNP selection -- use true causal SNPs"),
  make_option(c("--univar_allsnp"), type="logical", action="store_true", default=FALSE,
              help="Also run univariate MR with SNPs selected for all exposures"),
  make_option(c("--nome"), type="logical", action="store_true", default=FALSE,
              help="Enforce no measurement error in SNP effect estimates"),
  make_option(c("--annihilate_factors"), type="logical", action="store_true", default=FALSE,
              help="Project factors out of exposures and outcome instead of including as covariates"),
  make_option(c("--debug"), type="logical", action="store_true", default=FALSE,
              help="Activate debug mode -- writes a lot more output")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
N = opt$num_samples
K = opt$num_exposures
half = floor(N/2)
beta = as.numeric(unlist(strsplit(opt$beta, ",")))
methods = unlist(strsplit(tolower(opt$methods), ','))
run_famr = sum(grepl('famr', methods)) > 0
out = opt$out_prefix
iid_fname = paste0(out, '.TEMP.iids.txt')
merge_list_fname = paste0(out, '.TEMP.mergelist')
merged_pgen = paste0(out, '.TEMP.merged')


# pick a subset of loci and samples to use in simulations
subset_loci_and_samples = function(locidir, ld_dir, out, num_samp, num_loci) {
  all_loci = list.files(path = locidir, pattern = '*.pvar')
  all_ld = list.files(path = ld_dir, pattern = '*.ld.rds')
  if(num_loci < length(all_loci)) {
    selected_idx = sample(1:length(all_loci), num_loci)
  } else {
    selected_idx = 1:length(all_loci)
  }
  selected_loci = all_loci[selected_idx]
  selected_ld = all_ld[selected_idx]
  # selected_loci = sample(all_loci, num_loci)
  pvar = selected_loci[1]
  psam_file = paste0(locidir, substr(pvar, 1, nchar(pvar)-3), 'sam')
  psam_data = read.table(psam_file, header=T, sep=' ')
  all_iids = psam_data[,1]
  if(num_samp < length(all_iids)) {
    selected_iids = sort(sample(all_iids, num_samp))
  } else {
    selected_iids = sort(all_iids)
  }
  write.table(selected_iids, file = iid_fname,
              sep=' ', row.names=F, col.names=F, quote=F)
  return(list('loci' = selected_loci, 'ld' = selected_ld))
}


# randomly (uniformly) pick up to max_causal SNPs per locus to be causal; subset plink files
select_causal_snps = function(locidir, selected_loci, out, max_causal, max_snps=1000) {
  out_loci = c()
  all_samp_snps = c()  # all SNPs downsampled prior to selecting causal SNPs
  all_caus_snps = c()  # all causal SNPs (across all loci)
  all_samp_idx = c()
  for(locus in 1:length(selected_loci)) {
    # filter out individuals, subset SNPs by MAF threshold
    locus_name = paste0(locidir, selected_loci[locus])
    prefix = substr(locus_name, 1, nchar(locus_name)-5)
    pgen_out_allsnp = paste0(out, '.TEMP.allsnp.locus', locus)
    system(paste0('plink2 --pfile ', prefix, ' --keep ', iid_fname,
                  ' --make-pgen --maf 0.05 --out ', pgen_out_allsnp))
    locus_data = read.table(paste0(pgen_out_allsnp, '.pvar'))
    all_snps = locus_data[,3]
    system(paste0('rm ', out, '.TEMP.allsnp.locus*'))

    # select and write causal SNPs
    num_causal = floor(runif(1, 1, max_causal+1))  # pick num causal SNPs
    if(max_snps < length(all_snps)) {  # to save time, downsample SNPs
      samp_idx = sort(sample(1:length(all_snps), max_snps))
      samp_snps = all_snps[samp_idx]
      # samp_snps = sample(all_snps, max_snps)
    } else {
      samp_idx = 1:length(all_snps)
      samp_snps = all_snps
    }
    all_samp_idx = rbind(all_samp_idx, samp_idx)
    if(num_causal < length(samp_snps)) {
      selected_snps = sample(samp_snps, num_causal)
    } else {
      selected_snps = samp_snps
    }
    selected_snps = samp_snps[samp_snps %in% selected_snps]  # original order
    all_caus_snps = c(all_caus_snps, selected_snps)
    snps = paste0(out, '.TEMP.snps.txt')
    write.table(selected_snps, file = snps, sep=' ', row.names=F, col.names=F, quote=F)
    pgen_out = paste0(out, '.TEMP.locus', locus)
    system(paste0('plink2 --pfile ', prefix, ' --keep ', iid_fname,
                  ' --extract ', snps, ' --make-pgen --out ', pgen_out))
    out_loci = c(out_loci, pgen_out)
    all_samp_snps = rbind(all_samp_snps, samp_snps)
  }
  # write merged file with all selected causal SNPs
  if(length(out_loci) == 1) {
    system(paste0('plink2 --pfile ', out_loci[1], ' --export A-transpose --out ', merged_pgen))
  } else {
    tmp = out_loci[2:length(out_loci)]
    write.table(tmp, file = merge_list_fname, sep='\n', row.names=F, col.names=F, quote=F)
    system(paste0('plink2 --pfile ', out_loci[1], ' --pmerge-list ', merge_list_fname,
                  ' --indiv-sort none --export A-transpose --out ', merged_pgen))
  }
  system(paste0('rm ', out, '.TEMP.locus*'))
  return(list('snps' = all_samp_snps, 'caus_snps' = all_caus_snps, 'idx' = all_samp_idx))
}


# helper function to extract SNPs and MAF filter on plink file
extract_snps_maf_filter = function(locidir, sel_loci, all_samp_snps, idx, all_samp_idx=c()) {
  updated_idx = all_samp_idx
  snps = paste0(out, '.TEMP.snps.txt')
  samp_snps = all_samp_snps[idx,]
  write.table(samp_snps, file = snps, sep=' ', row.names=F, col.names=F, quote=F)
  plink_fname = paste0(locidir, sel_loci[idx])
  prefix = substr(plink_fname, 1, nchar(plink_fname)-5)
  pgen_out_allsnp = paste0(out, '.TEMP.allsnp.locus', idx)
  system(paste0('plink2 --pfile ', prefix, ' --keep ', iid_fname, ' --extract ', snps,
                ' --make-pgen --maf 0.05 --out ', pgen_out_allsnp))
  if(length(all_samp_idx) != 0) {  # need to update indices of SNPs after MAF filter
    new_snps = read.table(paste0(pgen_out_allsnp, '.pvar'))$V3
    updated_idx = all_samp_idx[samp_snps %in% new_snps]
    updated_idx = updated_idx[!is.na(updated_idx)]
  }
  return(list('fname' = pgen_out_allsnp, 'new_idx' = updated_idx))
}


# read in LD and SNP names and genotypes
read_ld_and_genos = function(plink_fname, two_samp) {
  pgen = NewPgen(paste0(plink_fname, '.pgen'))
  n_samp = GetRawSampleCt(pgen)
  snp_names = read.table(paste0(plink_fname, '.pvar'))$V3
  genos = ReadList(pgen, 1:length(snp_names))
  colnames(genos) = snp_names
  if(two_samp) {
    genos1 = genos[1:half, ]
    genos2 = genos[(half+1):n_samp, ]
    n_samp = half
    R2 = cor(genos2)
  } else {
    genos1 = genos
    genos2 = c()
    R2 = c()
  }
  R = cor(genos1)
  return(list('R' = R, 'R2' = R2, 'snp_names' = snp_names, 'N' = n_samp,
              'genos1' = genos1, 'genos2' = genos2))
}


# select SNPs as instruments using mvSusie or mvSusie-RSS
select_snp_mvsusie = function(do_rss, npheno, genos=c(), phenos=c(),
                              zscores=c(), ld=c(), N=0, susieL=10, pip_thresh=0.1,
                              sel=0) {
  if(!do_rss) {  # individual-level mode
    if(sel) {  # only select SNPs based on given exposure (testing runs only)
      res = susieR::susie(genos, as.numeric(phenos[,sel]), L=susieL)
    } else {
      res = mvsusie(genos, phenos, L=susieL, prior_variance=diag(npheno))
    }
  } else {  # summary statistics mode
    if(sel) {
      res = susieR::susie_rss(zscores[,sel], ld, N, L=susieL)
    } else {
      res = mvsusie_rss(as.matrix(zscores), ld, N, L=susieL, prior_variance=diag(npheno))
    }
  }
  pips = as.numeric(res$pip)
  # select by either pip threshold or max pip
  # selected_rows = pips > pip_thresh  # select by pip threshold
  selected_rows = 1:length(pips) %in% as.numeric(apply(res$alpha, 1, which.max))
  return(which(selected_rows))
}


# select SNP with 'peak' Zscore, either max or mean across all traits
select_snp_peak = function(zscores, criterion, sel=0, thresh=0.0) {
  if(sel) {  # only select SNPs based on given exposure (testing runs only)
    row = which.max(abs(zscores[,sel]))
    maxval = max(abs(zscores[,sel]))
  }
  else if(criterion == 'max') {  # get SNP with max abs(Z) score across all traits
    maxloc = which.max(abs(as.matrix(zscores)))
    maxval = max(abs(as.matrix(zscores)))
    numrows = dim(zscores)[1]
    row = maxloc %% numrows
    if(row == 0)  row = numrows  # if evenly divides, it's in last row, not "row 0"
  } else {  # get SNP with max abs(Z) averaged across traits
    abs_z_means = apply(abs(zscores), 1, mean)
    row = as.numeric(which.max(abs_z_means))
    maxval = max(abs_z_means)
  }
  # only select SNP if Z-score is above the threshold
  if(maxval > thresh) {
    return(row)
  } else {
    return(c())
  }
}


# generate betahats, stderrs, and zscores for a locus using RSS model,
#   given the true betas
run_rss = function(R, R2, snp_names, N, caus_snps, true_betas, nome) {
  # find causal SNP positions in locus, then fill out true_betas matrix with 0s
  #   for non-causal SNPs
  if(length(R2) == 0) {  # one sample mode
    R2 = R
  }
  M = length(snp_names)
  K = dim(true_betas)[2]
  causal_idx = match(caus_snps, snp_names)
  locus_true_betas = true_betas[!is.na(causal_idx),]  # only causal SNPs in this locus
  causal_idx = causal_idx[!is.na(causal_idx)]
  full_true_betas = matrix(0, M, K)
  full_true_betas[causal_idx,] = locus_true_betas
  S = diag(rep(1 / sqrt(N), M))

  # draw betahats from RSS model, fix stderrs = 1 / sqrt(N)
  if(nome) {  # if user specifies no measurement error, use true betas
    betahats_x = full_true_betas[,1:(K-1)]
  } else {
    mu_x = S %*% R %*% solve(S) %*% full_true_betas[,1:(K-1)]
    sigma_x = S %*% R %*% S
    betahats_x = apply(mu_x, 2, function(col) mvrnorm(1, col, sigma_x))
  }
  stderrs_x = matrix(1 / sqrt(N), M, K-1)
  zscores_x = betahats_x / stderrs_x

  # also draw for the outcome Y, which may come from a different sample; thus R2
  if(nome) {
    betahats_y = as.matrix(full_true_betas[,K])
  } else {
    mu_y = S %*% R2 %*% solve(S) %*% as.matrix(full_true_betas[,K])
    sigma_y = S %*% R2 %*% S
    betahats_y = mvrnorm(1, mu_y, sigma_y)
  }
  stderrs_y = rep(1 / sqrt(N), length(betahats_y))
  zscores_y = betahats_y / stderrs_y
  return(list('betas' = betahats_x, 'stderrs' = stderrs_x, 'Z' = zscores_x,
              'zscores_y' = zscores_y, 'betas_y' = betahats_y,
              'stderrs_y' = stderrs_y, 'pos' = snp_names))
}


# helper function to call snp selection for a locus and append results
select_and_append = function(res, univar, two_samp, locus_data, oracle, do_rss,
                             phenos, sumstats, susieL, pip_thresh, prune_thresh,
                             sel, criterion) {
  # if oracle (true effect SNPs) provided, use those; else run SNP selection
  if(length(oracle) > 0) {
    idx = which(colnames(locus_data$genos1) %in% oracle$all)
  } else {
    # run (mv)Susie or peak SNP selection -- return index of selected SNP(s)
    if(criterion == 'all') {
      idx = 1:ncol(locus_data$genos1)
    } else if(criterion == 'pruned' || criterion == 'clumped') {
      idx = FAMR:::ld_prune_famr(sumstats, locus_data$R, prune_thresh)
    } else if(grepl('susie', criterion)) {
      idx = select_snp_mvsusie(do_rss=do_rss, npheno=ncol(phenos), genos=locus_data$genos1,
                               phenos=phenos, zscores=sumstats$Z, ld=locus_data$R,
                               N=locus_data$N, susieL=susieL, pip_thresh=pip_thresh,
                               sel=sel)
    } else {
      idx = select_snp_peak(sumstats$Z, criterion, sel=sel)
    }
  }

  # append locus variables to the variables across all SNPs
  if(univar) {
    ci = as.character(sel)
    res$univ_G1[[ci]] = cbind(res$univ_G1[[ci]], locus_data$genos1[,idx])
    if(two_samp)  res$univ_G2[[ci]] = cbind(res$univ_G2[[ci]], locus_data$genos2[,idx])
  } else {
    res$G1 = cbind(res$G1, locus_data$genos1[,idx])
    if(two_samp)  res$G2 = cbind(res$G2, locus_data$genos2[,idx])
    res$names = c(res$names, locus_data$snp_names[idx])

    sumstats = subset_sumstats(sumstats, idx)
    res$sumstats = merge_sumstats(res$sumstats, sumstats)
  }
  return(res)
}


# generates summary statistics, writes them for FAMR,
#   and runs user-specified snp selection and famr filtering on all loci
gen_sumstats_sel_snps = function(phenos, caus_snps, true_betas, locidir, sel_loci,
                                 all_samp_snps, all_samp_idx, criterion='max',
                                 susieL=10, pip_thresh=0.1, prune_thresh=0.1,
                                 two_samp=F, univar=F, oracle=c(), nome=F) {
  res = list('G1' = c(), 'G2' = c(), 'univ_G1' = list(), 'univ_G2' = list(),
             'names' = c(), 'sumstats' = list(), 'ss_famr' = list(), 'R' = list())
  for(i in 1:ncol(phenos)) {
    res$univ_G1[[as.character(i)]] = c()
    res$univ_G2[[as.character(i)]] = c()
  }
  do_rss = !(criterion == 'mvsusie' || criterion == 'susie')  # flag for sumstats mode

  for(i in 1:length(sel_loci)) {
    # filter out individuals and subset SNPs by MAF threshold;
    #   get LD & genotypes; simulate summary statistics with RSS model
    filt = extract_snps_maf_filter(locidir, sel_loci, all_samp_snps, i, unique(all_samp_idx[i,]))
    locus_data = read_ld_and_genos(filt$fname, two_samp)
    system(paste0('rm ', filt$fname, '*'))
    sumstats = run_rss(locus_data$R, locus_data$R2, locus_data$snp_names,
                       locus_data$N, caus_snps, true_betas, nome)

    # if oracle SNP set provided, filter all other SNPs
    if(length(oracle) > 0) {  # use oracle if provided
      idx = which(colnames(locus_data$genos1) %in% oracle$all)
      ss_oracle = subset_sumstats(sumstats, idx)
      R_oracle = locus_data$R[idx, idx, drop=F]
      # append data for FAMR
      ss_oracle$locus_idx = rep(i, nrow(ss_oracle$betas))
      res$ss_famr = merge_sumstats(res$ss_famr, ss_oracle)
      res$R[[i]] = R_oracle
    } else {
      # append data for FAMR
      sumstats$locus_idx = rep(i, nrow(sumstats$betas))
      res$ss_famr = merge_sumstats(res$ss_famr, sumstats)
      res$R[[i]] = locus_data$R
    }

    # now run SNP selection and append the results for this locus to overall results
    res = select_and_append(res, F, two_samp, locus_data, oracle, do_rss,
                            phenos, sumstats, susieL, pip_thresh, prune_thresh,
                            0, criterion)
    # if running univariate MR, repeat SNP selection for each exposure
    if(univar) {
      for(i in 1:ncol(phenos)) {
        res = select_and_append(res, T, two_samp, locus_data, oracle, do_rss,
                                phenos, sumstats, susieL, pip_thresh, prune_thresh,
                                i, criterion)
      }
    }
  }
  return(res)
}


# compute true genetic effects on X and Y based on true parameters
compute_true_effts = function(par) {
  x_effts = par$theta$gx2 + (par$theta$gz %*% par$theta$zx)
  y_effts = par$theta$gy + (par$theta$gz %*% par$theta$zy) + x_effts %*% par$theta$xy
  comb_effts = cbind(x_effts, y_effts)
  gen_cor = cor(comb_effts)
  return(list('x' = x_effts, 'y' = y_effts, 'comb' = comb_effts, 'cor' = gen_cor))
}


# compute precision, recall, and F1 for SNP selection methods
compute_selection_metrics = function(true_names, pred_names, n_true, n_sel) {
  if(length(pred_names) == 0 || length(true_names) == 0) {
    return(list('precision' = 0, 'recall' = 0, 'f1' = 0))
  }
  tp = sum(pred_names %in% true_names)
  if(tp == 0) {
    return(list('precision' = 0, 'recall' = 0, 'f1' = 0))
  }
  fp = length(pred_names) - tp
  fn = length(true_names) - tp
  precision = tp / (tp + fp)
  recall = tp / (tp + fn)
  f1 = 2 * precision * recall / (precision + recall)

  print(paste0('Number of true causal SNPs: ', n_true))
  print(paste0('Number of selected SNPs: ', n_sel))
  print(paste0('Precision of SNP selection: ', precision))
  print(paste0('Recall of SNP selection: ', recall))
  print(paste0('F1 Score of SNP selection: ', f1))
  return(list('precision' = precision, 'recall' = recall, 'f1' = f1))
}


# simulate phenotypes based on user parameters and selected causal SNPs
sim_phenos_from_G = function(G_orig, opt) {
  sim_res = run_sim(G=G_orig, N=N, M=ncol(G_orig), K=K, J=opt$num_confounders,
                    H=opt$num_uncor, L=opt$num_med, beta=beta,
                    gamma_gx=opt$gamma_gx, gamma_gy=opt$gamma_gy, gamma_gz=opt$gamma_gz,
                    gamma_gu=opt$gamma_gu, gamma_gw=opt$gamma_gw,gamma_zx=opt$gamma_zx,
                    gamma_zy=opt$gamma_zy, gamma_uy=opt$gamma_uy, gamma_wx=opt$gamma_wx,
                    gamma_gx2=opt$gamma_gx2, gamma_yx2=opt$gamma_yx2, gamma_x1z=opt$gamma_x1z,
                    phi_gx=opt$phi_gx, phi_gz=opt$phi_gz, phi_gy=opt$phi_gy,
                    phi_gx2=opt$phi_gx2, phi_yx2=opt$phi_yx2, phi_x1z=opt$phi_x1z,
                    phi_zx=opt$phi_zx, phi_wx=opt$phi_wx, mu=opt$mu,
                    upstream_pct=opt$upstream_pct, reverse_cause=opt$reverse_cause)
  Z = sim_res$phe$Z
  X = sim_res$phe$X
  Y = sim_res$phe$Y
  true_effects = compute_true_effts(sim_res)
  Ztilde = sim_res$phe$Ztilde
  Xtilde = sim_res$phe$Xtilde  #G_orig %*% true_effects$x
  print('Genetic correlations of phenotypes:')
  print(true_effects$cor)

  # split samples in two if two sample mode is activated
  if(opt$two_sample) {
    # G is halved within select_snps
    X1 = X[1:half, ]
    X2 = X[(half+1):N, ]
    Y = Y[(half+1):N, ]  # only need second half (used in second stage)
    Z = Z[1:half, ]  # only need first half (used in first stage)
    Ztilde = Ztilde[(half+1):N, ]  # only need second half (oracle mode)
  } else {
    X1 = X
    X2 = c()
  }
  return(list('sim_res' = sim_res, 'X1' = X1, 'X2' = X2, 'Y' = Y, 'Z' = Z,
              'Ztilde' = Ztilde, 'Xtilde' = Xtilde, 'theta_gy' = sim_res$theta$gy,
              'true_effects' = true_effects))
}


# if user opts to skip snp selection, provide true causal SNP names
get_oracle = function(skip_selection, sim_res, caus_snps) {
  oracle = list()
  if(skip_selection) {
    all_idx = c()
    psi_x = sim_res$theta$gz %*% sim_res$theta$zx
    for(i in 1:K) {
      # get SNPs that affect exposure i
      effect_snps = sort(unique(c(which(sim_res$theta$gx2[,i] != 0), which(psi_x[,i] != 0))))
      oracle[[as.character(i)]] = caus_snps[effect_snps]
      all_idx = c(all_idx, effect_snps)
    }
    oracle$all = caus_snps
  }
  return(oracle)
}


# wrapper function to run non-FAMR methods, including univariate and sel options
run_methods_wrapper = function(methods, res, opt, X1, X2, Y, Z, Ztilde, theta_gy, true_effects) {
  mvmr_methods = methods
  univ_methods = mvmr_methods[which(grepl('univar', mvmr_methods))]
  if(length(univ_methods) > 0) {
    mvmr_methods = mvmr_methods[-which(grepl('univar', mvmr_methods))]
  }
  # already ran famr, so remove from list of methods
  famr_methods = which(grepl('famr', mvmr_methods))
  if(length(famr_methods) > 0) {
    mvmr_methods = mvmr_methods[-which(grepl('famr', mvmr_methods))]
  }
  if(length(mvmr_methods) > 0) {
    if(opt$nome && opt$skip_selection) {  # NOME: provide true x_betas and y_betas
      run_methods(mvmr_methods, res$G1, X1, Y, Z, out, G2=res$G2, X2=X2, Ztilde=Ztilde,
                  theta_gy=theta_gy, susieL=opt$susieL, verbose=TRUE,
                  x_betas=true_effects$x, y_betas=as.numeric(true_effects$y),
                  x_stderrs=(true_effects$x*0+0.01),
                  y_stderrs=as.numeric(true_effects$y*0+0.01))
    } else {
      run_methods(mvmr_methods, res$G1, X1, Y, Z, out, G2=res$G2, X2=X2, Ztilde=Ztilde,
                  theta_gy=theta_gy, susieL=opt$susieL, verbose=TRUE,
                  x_betas=res$sumstats$betas, y_betas=res$sumstats$betas_y,
                  x_stderrs=res$sumstats$stderrs, y_stderrs=res$sumstats$stderrs_y)
    }
  }
  # if we want results with SNP selection for a specific exposure, generate that
  if(opt$sel) {
    sel_methods = sapply(mvmr_methods, function(x) gsub('(_[^_]+)?$', '.sel\\1', x))
    if(length(mvmr_methods) > 0) {
      G1 = res$univ_G1[[as.character(opt$sel)]]
      G2 = res$univ_G2[[as.character(opt$sel)]]
      run_methods(sel_methods, G1, X1, Y, Z, out, G2=G2, X2=X2, Ztilde=Ztilde,
                  theta_gy=theta_gy, susieL=opt$susieL, verbose=TRUE)
    }
  }

  # run univariate SNP selection and methods, if needed
  if(length(univ_methods) > 0) {
    all_res = list()
    for(m in univ_methods)  all_res[[m]] = c()

    # run univariate MR on each exposure
    for(idx in 1:K) {
      G1 = res$univ_G1[[as.character(idx)]]
      G2 = res$univ_G2[[as.character(idx)]]
      idx_res = run_methods(univ_methods, G1, as.matrix(X1[,idx]), Y, Z, out,
                            G2=G2, X2=as.matrix(X2[,idx]),
                            susieL=opt$susieL, verbose=T, write.res=F)
      for(m in names(idx_res))  all_res[[m]] = rbind(all_res[[m]], idx_res[[m]])
    }

    # optionally run univariate MR on each exposure using SNPs selected for all exposures
    if(opt$univar_allsnp) {
      univ_methods_allsnp = sapply(
        univ_methods, function(x) gsub('(_[^_]+)?$', '.allsnp\\1', x))
      for(m in univ_methods_allsnp)  all_res[[m]] = c()
      for(idx in 1:K) {
        idx_res = run_methods(univ_methods_allsnp, res$G1, as.matrix(X1[,idx]),
                              Y, Z, out, G2=res$G2, X2=as.matrix(X2[,idx]),
                              susieL=opt$susieL, verbose=T, write.res=F)
        for(m in names(idx_res))  all_res[[m]] = rbind(all_res[[m]], idx_res[[m]])
      }
    }

    # write results
    for(m in names(all_res)) {
      all_res_m = as.numeric(all_res[[m]])
      write.table(t(all_res_m), file = paste0(out, m, '.txt'),
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}


# wrapper function to run FAMR (RSS mode)
run_famr_rss = function(ss_famr, R, annih_fact, prune_thresh, fa_prune_thresh,
                        final_prune_thresh, nome, theta_gz=c(), caus_snps=c()) {
  dat_to_return = data.frame()  # we return the GFA dat for accuracy checks
  # get only the FAMR methods, then extract the FA method used for each version
  famr_methods = methods[grepl('famr', methods)]
  fa_methods = sapply(famr_methods, function(x) unlist(strsplit(x, '_'))[2])

  # read data from input files (done once)
  sumstats = gather_sumstats_from_dat(ss_famr, R, prune_thresh = prune_thresh,
                                      fa_prune_thresh = fa_prune_thresh, 
                                      y_gwas = c(), oracle_mode = F)
  # if(debug)  saveRDS(sumstats, paste0(out, '_sumstats.rds'))
  if('famr_oracle' %in% famr_methods || 'famr_given' %in% famr_methods) {
    # no pruning in oracle mode -- don't want to prune true causal SNPs
    oracle_ss = gather_sumstats_from_dat(ss_famr, R, prune_thresh = 1.0,
                                        fa_prune_thresh = 1.0, 
                                        y_gwas = c(), oracle_mode = T)
  }

  # write oracle factors if appropriate
  oracle_fname = paste0(out, '.TEMP_oracle_factors.rds')
  if('famr_oracle' %in% famr_methods || 'famr_given' %in% famr_methods) {
    if(length(theta_gz) > 0 && max(abs(theta_gz)) != 0) {  # if there are any true factors
      factor_ss = list('betas' = theta_gz, 'pos' = caus_snps)
      factor_ss$stderrs = factor_ss$betas * 0 + (1 / sqrt(N))
      factor_ss$Z = factor_ss$betas / factor_ss$stderrs
      saveRDS(factor_ss, oracle_fname)
    }
  }

  # for each FA method, generate factors, add or annihilate them, and run FAMR
  for(fa in fa_methods) {
    if(!is.na(fa) && (fa == 'oracle' || fa == 'given') && length(theta_gz) > 0) {
      ss = oracle_ss
    } else {
      ss = sumstats
    }

    if(!is.na(fa) && (fa == 'oracle' || fa == 'given') && length(theta_gz) > 0) {
      # read in oracle factors
      factor_ss = generate_factors(fa, ss$ss_for_fa, N = N, 
                                   given_factors = oracle_fname, full_ss = ss)
    } else {
      # generate factors
      factor_ss = generate_factors(fa, ss$ss_for_fa, N = N, 
                                   given_factors = 'NONE', full_ss = ss)
      # if(debug)  saveRDS(factor_ss, paste0(out, '_', fa, '_factor_ss.rds'))
    }
    # annihilate factors if requested
    if(annih_fact) {
      ss = annihilate_factors(ss, factor_ss)
      # factor_ss = list('n_factors' = 0)
    }

    # precompute Z-scores and LD for FAMR, then run it and save results
    dat = learn_wts_precomp_merge(ss, factor_ss, N=N, nome=nome, 
                                  prune_thresh=final_prune_thresh)
    names(dat$sumstats$Z) = paste0('expo', 1:ncol(dat$sumstats$Z))
    if(!is.na(fa) && (fa == 'oracle' || (length(dat_to_return) == 0 && fa == 'gfa'))) {
      dat_to_return = dat
    }
    # famr_res = run_famr_rss_all(dat, opt$susieL, opt$famr_n_iter, 
    #                             K, N, annih=annih_fact)
    # write_famr_res(out, famr_res, K, fa, for_sim=T)
    famr_res = run_modified_ctwas(dat, opt$susieL, opt$famr_n_iter, 
                                K, N, annih=annih_fact)
    FAMR:::write_famr_res(out, famr_res, K, fa, for_sim=T)
  }
  return(dat_to_return)
}


# check the accuracy of the learned weights for exposures and factors
assess_famr_weights = function(famr_res, Xtilde, Ztilde, locidir, sel_loci,
                               all_samp_snps, all_samp_idx, two_samp=F) {
  Xtilde_pred = 0
  Ztilde_pred = 0
  K = ncol(Xtilde)

  for(i in 1:length(sel_loci)) {
    # filter out individuals and subset SNPs by MAF threshold;
    #   get LD & genotypes; simulate summary statistics with RSS model
    filt = extract_snps_maf_filter(locidir, sel_loci, all_samp_snps, i, unique(all_samp_idx[i,]))
    locus_data = read_ld_and_genos(filt$fname, two_samp)
    system(paste0('rm ', filt$fname, '*'))

    # get FAMR SNPs and their weights for this locus
    match_idx = na.omit(match(locus_data$snp_names, famr_res$prior$sumstats$all_ss$pos))
    snp_names = famr_res$prior$sumstats$all_ss$pos[match_idx]
    snp_wts = as.matrix(famr_res$prior$sumstats$all_ss$weights[match_idx, ])
    match_idx_fa = na.omit(match(locus_data$snp_names, famr_res$prior$sumstats$pos))
    snp_names_fa = famr_res$prior$sumstats$pos[match_idx_fa]
    snp_wts_fa = as.matrix(famr_res$prior$sumstats$weights[match_idx_fa, ])
    if(length(match_idx) == 1) {
      snp_wts = t(snp_wts)
      snp_wts_fa = t(snp_wts_fa)
    }
    # now fill in all weights for this locus, with zeros for SNPs filtered by FAMR
    idx_wts = matrix(0, length(locus_data$snp_names), ncol(snp_wts))
    idx_of_snps = c(na.omit(match(snp_names, locus_data$snp_names)))
    idx_wts[idx_of_snps, ] = snp_wts
    idx_wts_fa = matrix(0, length(locus_data$snp_names), ncol(snp_wts_fa))
    idx_of_snps_fa = c(na.omit(match(snp_names_fa, locus_data$snp_names)))
    idx_wts_fa[idx_of_snps_fa, ] = snp_wts_fa

    # now get the predicted values for Xtilde and Ztilde for this locus
    #   based on the given weights, and add to the running total
    pred_wts = locus_data$genos1 %*% idx_wts
    pred_wts_fa = locus_data$genos1 %*% idx_wts_fa
    idx_Xtilde_pred = pred_wts[, 1:K]
    if(K != ncol(pred_wts_fa)) {  # if factors were inferred
      idx_Ztilde_pred = pred_wts_fa[, (K+1):ncol(pred_wts_fa)]
    } else {
      idx_Ztilde_pred = 0
    }
    Xtilde_pred = Xtilde_pred + idx_Xtilde_pred
    Ztilde_pred = Ztilde_pred + idx_Ztilde_pred
  }

  # now compute r^2 of predicted X/Ztilde with true X/Ztilde, write results
  if(two_samp) {
    Xtilde = Xtilde[1:half,]
    Ztilde = Ztilde[1:half,]
  }
  Xtilde_cor = cor(Xtilde, Xtilde_pred) ** 2
  if(length(Ztilde_pred) > 1) {  # if factors were inferred
    Ztilde_cor = cor(Ztilde, Ztilde_pred) ** 2
  } else {
    Ztilde_cor = 0
  }
  dat = list('Xtilde_cor' = Xtilde_cor, 'Ztilde_cor' = Ztilde_cor)
  saveRDS(dat, paste0(out, 'XZtilde_cor.rds'))
}


if (sys.nframe() == 0) {
for (sim in 1:opt$num_sims) {
  # subset loci and samples
  sel = subset_loci_and_samples(opt$in_dir, opt$ld_dir, out, N, opt$num_loci)
  # pick causal SNPs
  all_samp = select_causal_snps(opt$in_dir, sel$loci, out, opt$max_causal,
                                max_snps=opt$snps_per_locus)
  caus_snps = all_samp$caus_snps
  # read in genotypes, undo plink re-ordering
  genos_raw = read.table(paste0(merged_pgen, '.traw'), sep='\t', header=T)
  G_orig = t(as.matrix(sapply(genos_raw[,7:dim(genos_raw)[2]], as.numeric)))
  if(nrow(G_orig) == 1)  G_orig = t(G_orig)  # need to re-flip if 1 SNP
  pvar = fread(paste0(merged_pgen, '.pvar'))
  ord = match(caus_snps, pvar$ID)
  G_orig = G_orig[, ord]

  # Simulate phenotypes given the picked causal SNPs
  phe = sim_phenos_from_G(G_orig, opt)
  if(opt$debug) {
    phe$caus_snps = caus_snps
    saveRDS(phe, paste0(out, '_phe.rds'))
  }

  # generate SNP summary statistics, select SNPs for input to traditional methods,
  #   write sumstats for FAMR
  any_univar = (sum(sapply(methods, function(x) grepl('univar', x))) > 0)
  if(opt$sel)  any_univar=T  # will gather SNPs based only on selected expo
  oracle = get_oracle(opt$skip_selection, phe$sim_res, caus_snps)
  res = gen_sumstats_sel_snps(phe$X1, caus_snps, phe$true_effects$comb,
                              opt$in_dir, sel$loci, all_samp$snp, all_samp$idx,
                              criterion=tolower(opt$criterion), susieL=opt$susieL,
                              pip_thresh=opt$pip_thresh, prune_thresh=opt$famr_prune,
                              two_samp=opt$two_sample, univar=any_univar,
                              oracle=oracle, nome=opt$nome)
  # evaluate how well SNP selection did
  sel_met = compute_selection_metrics(caus_snps, res$names, ncol(G_orig), ncol(res$G1))

  # run methods and save results
  famr_res = run_famr_rss(res$ss_famr, res$R, opt$annihilate_factors,
                          opt$famr_prune, opt$fa_prune_thresh, opt$final_prune_thresh,
                          opt$nome, theta_gz=phe$sim_res$theta$gz, caus_snps=caus_snps)
  run_methods_wrapper(methods, res, opt, phe$X1, phe$X2, phe$Y, phe$Z,
                      phe$Ztilde, phe$theta_gy, phe$true_effects)

  # assess accuracy of learned famr weights
  if(opt$debug && length(famr_res) != 0)   {
    assess_famr_weights(famr_res, phe$Xtilde, phe$Ztilde, opt$in_dir, sel$loci,
                        all_samp$snp, all_samp$idx, two_samp=opt$two_sample)
  }

  # write out true correlations between phenotypes
  if(opt$debug) {
    if(!opt$two_sample) {
      cor_xz = cor(cbind(phe$Xtilde, phe$Ztilde))
    } else {
      cor_xz = cor(cbind(phe$Xtilde[1:half,], phe$Ztilde[1:half,]))
    }
    saveRDS(cor_xz, paste0(out, '.trueR.rds'))
  }

  # clean up temporary files
  system(paste0('rm ', out, '.TEMP*'))
}
}
