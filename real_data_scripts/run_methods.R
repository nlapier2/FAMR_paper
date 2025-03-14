# Given data, runs MR methods on it

library(mvsusieR)
library(optparse)
source("methods.R")


# check whether enough information has been provided
check_run_inputs = function(methods, G, X, Y, x_betas, x_stderrs, y_betas, y_stderrs, N) {
  if(length(methods) == 0) {
    stop('Error: you must provide at least one method to run.')
  }

  # various checks that if we are in summary statistics mode that we have
  #   adequate information and are not running individual-level methods
  if(length(G) == 0 || length(X) == 0 || length(Y) == 0) {
    if(length(x_betas) == 0 || length(x_stderrs) == 0 ||
       length(y_betas) == 0 || length(y_stderrs) == 0 || N == 0) {
      stop('Error: you must provide all of G, X, and Y or all of x/y betas & stderrs and N.')
    }
    if(sum(grepl('cate', methods)) > 0  || sum(grepl('2sls', methods)) > 0  ||
       sum(grepl('tsls', methods)) > 0) {
      stop('Error: cannot run Cate or two stage least squares unless G, X, Y provided.')
    }
    for(i in 1:length(methods)) {
      second_stage = unlist(strsplit(methods[i], '_'))[1]
      if(second_stage == 'susie' || second_stage == 'mrash') {
        stop('Error: cannot run individual-level susie or mr.ash unless G, X, Y provided.')
      }
    }
  }

}


# helper function to generate covariates for downstream analysis given inferred factors
gen_covariates = function(factors, x_betas, x_stderrs, G, G2, N=10000) {
  if(length(factors) == 0) {
    ztildehat = c()
    aug_x_betas = x_betas
    aug_x_stderrs = x_stderrs
  } else {
    if(length(G) == 0 && length(G2) == 0) {
      ztildehat = c()
    } else {
      if(length(G2) == 0) {
        N = dim(G)[1]
        ztildehat = G %*% factors
      } else {
        N = dim(G2)[1]
        ztildehat = G2 %*% factors
      }
    }
    M = dim(x_betas)[1]
    J = dim(factors)[2]
    aug_x_betas = cbind(x_betas, factors)
    imp_stderr = matrix(rep(1/sqrt(N), M*J), M, J)
    aug_x_stderrs = cbind(x_stderrs, imp_stderr)
  }
  return(list('ztildehat' = ztildehat, 'factors' = factors,
              'aug_x_betas' = aug_x_betas, 'aug_x_stderrs' = aug_x_stderrs))
}


# like gen_covariates, but specifically for the oracle case
gen_covariates_oracle = function(x_betas, x_stderrs, G, G2, Z, Ztilde, theta_gy) {
  oracle = list()
  if(length(G2) == 0) {  # one sample setting
    GthetaY = G %*% theta_gy
  } else {  # two/three sample setting
    GthetaY = G2 %*% theta_gy
  }
  if(max(Ztilde) == 0 && min(Ztilde) == 0) {
    if(max(GthetaY) == 0 && min(GthetaY) == 0) {
      confounder = c()
    } else {
      confounder = GthetaY
    }
  } else if(max(GthetaY) == 0 && min(GthetaY) == 0) {
    confounder = Ztilde
  } else {
    confounder = cbind(Ztilde, GthetaY)
  }
  regr_res = get_betas_stderrs(G, Z, Z[,1])
  z_betas = regr_res$x_betas
  z_stderrs = regr_res$x_stderrs
  oracle$ztildehat = confounder
  oracle$factors = z_betas
  oracle$aug_x_betas = cbind(x_betas, z_betas)
  oracle$aug_x_stderrs = cbind(x_stderrs, z_stderrs)
  return(oracle)
}


# prune SNPs above certain LD based on R^2 threshold
ld_prune = function(G, G2, x_betas, x_stderrs, y_betas, y_stderrs, prior_weights,
                    thresh, verbose) {
  if(length(G) == 0) {
    print('Warning: LD pruning currently not supported in summary statistics mode.')
    return(list(G = G, G2 = G2, x_betas = x_betas, x_stderrs = x_stderrs,
         y_betas = y_betas, y_stderrs = y_stderrs, prior_weights = prior_weights))
  }
  if(verbose) {
    print(paste0('LD pruning variables at R^2 threshold of: ', thresh))
  }
  exclude = c()
  mat = cbind(G, G2)
  r2mat = cor(mat) ** 2
  for(i in 1:dim(r2mat)[1]) {
    for(j in 1:(i-1)) {
      if(0 < j && j < i && r2mat[i,j] > thresh) {
        exclude = c(exclude, i)
      }
    }
  }
  exclude = unique(exclude)

  if(length(exclude) > 0) {
    G = G[,-exclude]
    G2 = G2[,-exclude]
    if(length(x_betas) > 0) {
      x_betas = x_betas[-exclude, ]
      x_stderrs = x_stderrs[-exclude, ]
    }
    if(length(y_betas) > 0) {
      y_betas = y_betas[-exclude]
      y_stderrs = y_stderrs[-exclude]
    }
    if(length(prior_weights) > 0) {
      prior_weights = prior_weights[-exclude]
    }
  }
  if(verbose) {
    print(paste0('Pruned ', length(exclude), ' SNPs; ', dim(G)[2], ' left.'))
  }
  return(list(G = G, G2 = G2, x_betas = x_betas, x_stderrs = x_stderrs,
              y_betas = y_betas, y_stderrs = y_stderrs,
              prior_weights=prior_weights, exclude=exclude))
}


# simple helper method to print info about inferred factors (in verbose mode)
print_factor_info = function(factors, method, x_betas) {
  n_factors = dim(factors)[2]
  if(!is.null(n_factors)) {
    print(paste0('Number of factors detected by ', method, ': ', n_factors))
    print('Correlations with x_betas:')
    print(cor(factors, x_betas))
  } else {
    print(paste0('Number of factors detected by ', method, ': 0'))
  }
  print('')
}


# run first stage factor analysis
run_first_stage = function(methods, G, G2, X, x_betas, x_stderrs, N,
                           Z=c(), Ztilde=c(), theta_gy=c(), verbose=FALSE) {
  fs = list('gfa' = list(), 'rpca' = list(), 'spca' = list(), 'cate' = list(), 'oracle' = list())
  if(verbose) { print('Performing first-stage factor analysis...') }

  if(sum(grepl('cate', methods)) > 0 ) {
    factors = run_cate(G, X)
    fs$cate = gen_covariates(factors, x_betas, x_stderrs, G, G2, N)
    if(verbose) { print_factor_info(fs$cate$factors, 'Cate', x_betas) }
  }
  if(sum(grepl('rpca', methods)) > 0 ) {
    factors = run_rpca(x_betas, x_stderrs)
    fs$rpca = gen_covariates(factors, x_betas, x_stderrs, G, G2, N)
    if(verbose) { print_factor_info(fs$rpca$factors, 'RPCA', x_betas) }
  }
  if(sum(grepl('spca', methods)) > 0 ) {
    factors = run_spca(x_betas, x_stderrs)
    fs$spca = gen_covariates(factors, x_betas, x_stderrs, G, G2, N)
    if(verbose) { print_factor_info(fs$spca$factors, 'SPCA', x_betas) }
  }
  if(sum(grepl('gfa', methods)) > 0 ) {
    factors = run_gfa(x_betas, x_stderrs, N)
    fs$gfa = gen_covariates(factors, x_betas, x_stderrs, G, G2, N)
    if(verbose) { print_factor_info(fs$gfa$factors, 'GFA', x_betas) }
  }
  if(sum(grepl('oracle', methods)) > 0 ) {
    fs$oracle = gen_covariates_oracle(x_betas, x_stderrs, G, G2, Z, Ztilde, theta_gy)
  }
  return(fs)
}


# load appropriate first stage results into standard variable names
load_first_stage = function(first_stage, fs, x_betas, x_stderrs) {
  this.fs = list(ztildehat = c(), aug_x_betas = x_betas, aug_x_stderrs = x_stderrs)
  if(!is.na(first_stage)) {
    if(first_stage == 'cate') {
      this.fs$ztildehat = fs$cate$ztildehat
      this.fs$aug_x_betas = fs$cate$aug_x_betas
      this.fs$aug_x_stderrs = fs$cate$aug_x_stderrs
    } else if(first_stage == 'rpca') {
      this.fs$ztildehat = fs$rpca$ztildehat
      this.fs$aug_x_betas = fs$rpca$aug_x_betas
      this.fs$aug_x_stderrs = fs$rpca$aug_x_stderrs
    } else if(first_stage == 'spca') {
      this.fs$ztildehat = fs$spca$ztildehat
      this.fs$aug_x_betas = fs$spca$aug_x_betas
      this.fs$aug_x_stderrs = fs$spca$aug_x_stderrs
    } else if(first_stage == 'gfa') {
      this.fs$ztildehat = fs$gfa$ztildehat
      this.fs$aug_x_betas = fs$gfa$aug_x_betas
      this.fs$aug_x_stderrs = fs$gfa$aug_x_stderrs
    } else if(first_stage == 'oracle') {
      this.fs$ztildehat = fs$oracle$ztildehat
      this.fs$aug_x_betas = fs$oracle$aug_x_betas
      this.fs$aug_x_stderrs = fs$oracle$aug_x_stderrs
    }
  }
  return(this.fs)
}


run_second_stage = function(second_stage, G, G2, X, X2, Y, x_pred, this.fs, ss,
                            corg, susieL, N, K, prior_weights, prior_sigma2) {
  # extract some variables for simplicity / readability
  x_betas = this.fs$aug_x_betas
  x_stderrs = this.fs$aug_x_stderrs
  y_betas = ss$y_betas
  y_stderrs = ss$y_stderrs
  ztildehat = this.fs$ztildehat
  mr_methods = c('ivw', 'robust', 'median', 'egger', 'lasso')
  if(grepl('.univar', second_stage)) {  # strip off univariate designator
    second_stage = unlist(strsplit(second_stage, ".univar"))[1]
  } else if(grepl('.sel', second_stage)) {  # strip off SNP selection designator
    second_stage = unlist(strsplit(second_stage, ".sel"))[1]
  }

  # run second stage methods
  if(second_stage == '2sls' || second_stage == 'tsls') {
    res = run_2sls(G, X, Y, covariates=ztildehat, G2=G2, X2=X2, x_pred=x_pred)
  } else if(second_stage == '2sls.as' || second_stage == 'tsls.as') {
    weights = as.matrix(as.matrix(lm(X ~ G)$coefficients)[2:(dim(G)[2]+1),])
    res = run_allele_score(G, X, Y, weights, covariates=ztildehat, G2=G2, X2=X2)
  } else if(second_stage %in% mr_methods) {
    res = run_mr(second_stage, x_betas, x_stderrs, y_betas, y_stderrs, corg, 1:K)
  } else if(second_stage == 'susie') {
    res = run_susie_mrash(X, Y, G=G, L=susieL, ztildehat=ztildehat, G2=G2,
                          x_pred=x_pred, method='susie',
                          prior_weights=prior_weights, prior_sigma2=prior_sigma2)
  } else if(second_stage == 'mrash' || second_stage == 'mr.ash') {
    res = run_susie_mrash(X, Y, G=G, ztildehat=ztildehat, G2=G2,
                          x_pred=x_pred, method='mrash')
  } else if(second_stage == 'varbvs') {
    res = run_susie_mrash(X, Y, G=G, ztildehat=ztildehat, G2=G2,
                          x_pred=x_pred, method='varbvs', L=susieL,
                          prior_weights=prior_weights, prior_sigma2=prior_sigma2)
  } else if(second_stage == 'susie.ss') {
    res = run_susie_mrash(x_betas, y_betas, L=susieL, indices=1:K, method='susie')
  } else if(second_stage == 'mrash.ss' || second_stage == 'mr.ash.ss') {
    res = run_susie_mrash(x_betas, y_betas, indices=1:K, method='mrash')
  } else if(second_stage == 'varbvs.ss') {
    res = run_susie_mrash(x_betas, y_betas, indices=1:K, method='varbvs')
  } else if(second_stage == 'vebboost' || second_stage == 'veb.boost') {
    res = run_vebboost(G, X, Y, ztildehat=ztildehat, G2=G2, x_pred=x_pred)
  } else if(second_stage == 'brms.hs') {
    res = run_brms(x_betas, y_betas, 1:K, prior='horseshoe')
  } else if(second_stage == 'brms') {
    res = run_brms(x_betas, y_betas, 1:K, prior='normal')
  } else if(second_stage == 'grapple' || second_stage == 'GRAPPLE') {
    res = run_grapple(G, X, Y, 1:K, x_betas, x_stderrs, y_betas, y_stderrs)
  } else if(second_stage == 'bma' || second_stage == 'zuber') {
    res = run_bma(G, X, Y, 1:K, x_betas, x_stderrs, y_betas, y_stderrs)
  } else if(second_stage == 'cML' || second_stage == 'cml') {
    mat = cbind(x_betas, y_betas)
    rho_mat = cor(mat)
    res = run_cml(x_betas, x_stderrs, y_betas, y_stderrs, rho_mat, N, 1:K, dp=F)
  } else if(second_stage == 'cML.DP' || second_stage == 'cml.dp') {
    mat = cbind(x_betas, y_betas)
    rho_mat = cor(mat)
    res = run_cml(x_betas, x_stderrs, y_betas, y_stderrs, rho_mat, N, 1:K, dp=T)
  } else {
    return(c())
  }
  return(res)
}


run_methods = function(methods, G=c(), X=c(), Y=c(), Z=c(), out='results_',
                       Ztilde=c(), theta_gy=c(), G2=c(), X2=c(),
                       x_betas=c(), x_stderrs=c(), y_betas=c(), y_stderrs=c(),
                       prior_weights=c(), prior_sigma2=c(), susieL=0,
                       prune_thresh=0.95, N=0, corg=c(), verbose=F, write.res=T) {
  check_run_inputs(methods, G, X, Y, x_betas, x_stderrs, y_betas, y_stderrs, N)

  # LD prune variables
  vars = ld_prune(G, G2, x_betas, x_stderrs, y_betas, y_stderrs, prior_weights,
                  prune_thresh, verbose)
  G = vars$G
  G2 = vars$G2

  # precompute univariate regressions if not provided
  if(length(y_betas) == 0) {
    if(verbose) { print('Performing genotype-phenotype regressions...') }
    ss = get_betas_stderrs(G, X, Y, G2, x_betas = vars$x_betas, x_stderrs = vars$x_stderrs)
  } else {
    ss = vars
  }

  # set some variables
  K = dim(ss$x_betas)[2]
  if(susieL == 0) { susieL = dim(ss$x_betas)[2] }
  if(length(G) == 0 && length(corg) == 0) {
    corg = diag(dim(ss$x_betas)[1])
  } else if(length(G) != 0) {
    corg = cor(G)
    if(N == 0) {
      N = dim(X)[1]
    }
  }

  # precompute X predicted by G
  if(length(G)  == 0) {
    x_pred = c()
  } else {
    if(verbose) { print('Performing first-stage multiple regressions...') }
    if(length(G2) == 0) {  # one sample setting
      x_pred = as.matrix(lm(X ~ G)$fitted.values)
    } else {  # two/three sample setting
      theta_gx_hat = as.matrix(as.matrix(lm(X ~ G)$coefficients)[2:(dim(G)[2]+1),])
      x_pred = G2 %*% theta_gx_hat
      # X = X2
    }
  }

  # precompute first-stage estimates of confounder Ztilde, if used
  fs = run_first_stage(methods, G, G2, X, ss$x_betas, ss$x_stderrs, N,
                       Z, Ztilde, theta_gy, verbose)

  # Run MR methods and write results
  if(verbose) { print('Running 2nd stage methods...') }
  all_res = list()
  for(i in 1:length(methods)) {
    if(verbose) { print(paste0('Running ', methods[i], '...')) }
    met = methods[i]
    splits = unlist(strsplit(met, '_'))
    second_stage = splits[1]
    first_stage = splits[2]

    # load appropriate first-stage estimate for confounder Ztilde and/or factors
    this.fs = load_first_stage(first_stage, fs, ss$x_betas, ss$x_stderrs)

    # run second stage and write results
    res = run_second_stage(second_stage, G, G2, X, X2, Y, x_pred, this.fs, ss,
                           corg, susieL, N, K, vars$prior_weights, prior_sigma2)
    if(write.res && length(res) > 0) {
      write.table(t(res), file = paste0(out, met, '.txt'),
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    all_res[[met]] = res
  }
  return(all_res)
}
