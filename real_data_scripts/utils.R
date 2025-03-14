# Provides various utility functions used by multiple scripts
library(data.table)
library(dplyr)
library(flashr)
library(GFA)
library(rpca)
library(elasticnet)


# helper function to subset summary statistics based on provided indices
subset_sumstats = function(ss, idx) {
  ss$betas = ss$betas[idx, , drop=F]
  ss$stderrs = ss$stderrs[idx, , drop=F]
  ss$weights = ss$weights[idx, , drop=F]
  ss$Z = ss$Z[idx, , drop=F]
  if(is.null(nrow(ss$pos))) {
    ss$pos = ss$pos[idx]
  } else {
    ss$pos = ss$pos[idx, , drop=F]
  }
  if('zscores_y' %in% names(ss)) {
    ss$zscores_y = ss$zscores_y[idx]
    ss$betas_y = ss$betas_y[idx]
    ss$stderrs_y = ss$stderrs_y[idx]
  }
  if('locus_idx' %in% names(ss)) {
    ss$locus_idx = ss$locus_idx[idx]
  }
  if('annih_y' %in% names(ss)) {
    ss$annih_y = ss$annih_y[idx]
    ss$orig_betas_y = ss$orig_betas_y[idx]
    ss$orig_zscores_y = ss$orig_zscores_y[idx]
  }
  return(ss)
}


# merge sumstats for a locus onto overall sumstats
merge_sumstats = function(all_sumstats, locus_sumstats, by_row=T) {
  if(length(all_sumstats) == 0 || !('betas' %in% names(all_sumstats))) {
    return(locus_sumstats)
  }
  if(by_row) {
    all_sumstats$Z = rbind(all_sumstats$Z, locus_sumstats$Z)
    all_sumstats$betas = rbind(all_sumstats$betas, locus_sumstats$betas)
    all_sumstats$stderrs = rbind(all_sumstats$stderrs, locus_sumstats$stderrs)
    all_sumstats$weights = rbind(all_sumstats$weights, locus_sumstats$weights)
  } else {
    all_sumstats$Z = cbind(all_sumstats$Z, locus_sumstats$Z)
    all_sumstats$betas = cbind(all_sumstats$betas, locus_sumstats$betas)
    all_sumstats$stderrs = cbind(all_sumstats$stderrs, locus_sumstats$stderrs)
    all_sumstats$weights = cbind(all_sumstats$weights, locus_sumstats$weights)
  }
  if(is.null(nrow(all_sumstats$pos)) && by_row) {
    all_sumstats$pos = c(all_sumstats$pos, locus_sumstats$pos)
  } else if(by_row) {
    all_sumstats$pos = rbind(all_sumstats$pos, locus_sumstats$pos)
  }
  if('zscores_y' %in% names(all_sumstats)) {
    all_sumstats$zscores_y = c(all_sumstats$zscores_y, locus_sumstats$zscores_y)
    all_sumstats$betas_y = c(all_sumstats$betas_y, locus_sumstats$betas_y)
    all_sumstats$stderrs_y = c(all_sumstats$stderrs_y, locus_sumstats$stderrs_y)
  }
  if('locus_idx' %in% names(all_sumstats)) {
    all_sumstats$locus_idx = c(all_sumstats$locus_idx, locus_sumstats$locus_idx)
  }
  if('lambda' %in% names(all_sumstats) && by_row) {
    all_sumstats$lambda = rbind(all_sumstats$lambda, locus_sumstats$lambda)
  }
  return(all_sumstats)
}


# read in VCF file, clean data, and preprocess for use in real data analysis
read_vcf = function(fname) {
  # read gwas file, filter for NA, MAF, multi-character allele, duplicates, etc
  gwas = fread(paste0("zgrep -v '^##' ", fname))
  setnames(gwas, '#CHROM', 'CHR')

  # Filter NAs, non-passing QC SNPs, reverse-complements, MHC region, and multi-allelic SNPs
  gwas = gwas %>% na.omit() %>%
    filter(FILTER == 'PASS',
           nchar(REF) == 1, nchar(ALT) == 1,
           !(INFO == "ReverseComplementedAlleles"),
           !(CHR == 6 & POS>=25000000 & POS<=36000000))


  # extract allele frequency (AF) from INFO field, filter at MAF < 0.01
  if(grepl('AF', gwas[['INFO']][1])) {
    AF = as.numeric(sapply(gwas$INFO, function(x) unlist(strsplit(x, '='))[2]))
    gwas = gwas[AF > 0.01 & AF < 0.99, ]  # MAF filter at 0.01
  }

  # gather the betas and stderrs from the last column using the FORMAT column
  #   to tell us which columns tell us those
  format_split = unlist(strsplit(gwas[1,]$FORMAT, ':'))
  beta_col = which(format_split == 'ES')
  se_col = which(format_split == 'SE')
  logp_col = which(format_split == 'LP')
  ss_col = gwas[[names(gwas)[length(gwas)]]]
  ss_split = strsplit(ss_col, ':')
  gwas$BETA = as.numeric(lapply(ss_split, function(x) x[beta_col]))
  gwas$SE = as.numeric(lapply(ss_split, function(x) x[se_col]))
  gwas$logp = as.numeric(lapply(ss_split, function(x) x[logp_col]))
  gwas$Z = gwas$BETA / gwas$SE
  gwas = gwas %>% na.omit()
  return(gwas)
}


# helper function for getting p-values from a wald test
getp = function(val, df) {
  return ((1-pt(val, df)) + pt(-val, df))
}

# run truncated svd to infer latent factors
get_svd_z = function(mat, trunc=5) {
  if(dim(mat)[2] == 1) {  # do not run SVD on univariate matrix
    return(mat)
  }
  s = svd(mat)
  if(trunc == 0){
    svd_z = s$u %*% diag(s$d) %*% t(s$v)
    # svd_z = mat %*% s$v
  } else {
    trunc_u = s$u[,1:trunc]
    trunc_d = diag(s$d[1:trunc])
    trunc_v = t(s$v)[1:trunc,]
    svd_z = trunc_u %*% trunc_d %*% trunc_v
    svd_z = svd_z[,1:trunc]
    # svd_z = mat %*% trunc_v
  }
  return(svd_z)
}

# helper function to return R^2 for all columns in mat1 explained by mat2
get_allr2 = function(mat1, mat2) {
  allr2=c()
  nvar1 = dim(mat1)[2]
  for(i in 1:nvar1) {
    r2 = summary(lm(mat1[,i] ~ mat2 + 0))$r.squared
    allr2 = c(allr2, r2)
  }
  return(allr2)
}


# helper function to heuristically prune excess factors for GFA: prune factors
#   that don't explain at least "thresh" variance of at least "min_count" exposures
gfa_factor_prune = function(x_betas, gfa_factors, thresh = 0.1, min_count = 2) {
    r2_matrix = cor(x_betas, gfa_factors)^2
    keep_indices = which(apply(r2_matrix, 2, function(x) sum(x > thresh)) >= min_count)
    kept_factors = gfa_factors[, keep_indices, drop = FALSE]
    return(kept_factors)
}


# equivalent to gfa_factor_prune but for full GFA object
gfa_factor_prune_full = function(x_betas, gfares, thresh = 0.1, min_count = 2) {
  r2_matrix = cor(x_betas, gfares$L_hat)^2
  keep_indices = which(apply(r2_matrix, 2, function(x) sum(x > thresh)) >= min_count)
  gfares$L_hat = gfares$L_hat[, keep_indices]
  gfares$F_hat = gfares$F_hat[, keep_indices]
  gfares$gfa_pve$pve = gfares$gfa_pve$pve[, keep_indices]
  return(gfares)
}


# run gfa to estimate confounders
run_gfa = function(x_betas, x_stderrs, N) {
  gfa_factors = tryCatch({
    gfares = gfa_fit(B_hat = x_betas, S = x_stderrs)
    return(gfa_factor_prune(x_betas, gfares$L_hat))
  }, error = function(e) {
    print(paste0('Error in GFA (may simply be no factors identified): ', e))
    return(c())
  })
  return(gfa_factors)
}


# run gfa to estimate confounders, and return the full GFA object
run_gfa_full = function(x_betas, x_stderrs, N) {
  gfa_factors = tryCatch({
    gfares = gfa_fit(B_hat = x_betas, S = x_stderrs)
    return(gfa_factor_prune_full(x_betas, gfares))
  }, error = function(e) {
    print(paste0('Error in GFA (may simply be no factors identified): ', e))
    return(c())
  })
  return(gfa_factors)
}


# run cate to estimate confounders
run_cate = function(myG, myX) {
  df = data.frame(X=myX,G=myG)
  trunc = est.confounder.num(X ~ G + 0, df, X, method='bcv', rmax = 30, nRepeat = 20)$r
  if(trunc == 0) {
    return(c())
  }
  cateres = cate(X ~ G + 0, df, X, r = trunc, adj.method='rr')
  cate_factors = t(cateres$alpha)
  return(cate_factors)
}


# run rpca to estimate confounders
run_rpca = function(x_betas, x_stderrs) {
  rpcares = rpca(x_betas)
  flash_xb = flash(x_betas)
  nfactors = flash_xb$nfactors
  if(nfactors == 0) {
    return(c())
  } else {
    svd_v = svd(t(rpcares$L))$v
    if(nfactors > dim(svd_v)[2]) {
      rpca_factors = svd_v
    } else {
      rpca_factors = svd_v[,1:nfactors]
    }
  }
  return(rpca_factors)
}

# run rpca to estimate confounders, and include both factors and loadings
run_rpca_full = function(x_betas, x_stderrs) {
  rpcares = rpca(x_betas)
  flash_xb = flash(x_betas)
  nfactors = flash_xb$nfactors
  if(nfactors == 0) {
    return(c())
  } else {
    svd_res = svd(t(rpcares$L))
    if(nfactors <= dim(svd_res$v)[2]) {
      svd_res$v = svd_res$v[,1:nfactors]
      svd_res$u = svd_res$u[,1:nfactors]
    }
  }
  return(list('L' = svd_res$v, 'F' = svd_res$u))
}


# run sparse pca
run_spca = function(x_betas, x_stderrs, full=FALSE) {
  # first use flash to estimate number of factors
  flash_xb = flash(x_betas)
  nfactors = flash_xb$nfactors
  if(nfactors == 0) {
    return(c())
  }

  # set penalty parameter, then run
  para = c(0.06,0.16,0.1, rep(0.5, 999))
  spca_res = spca(x_betas, K=nfactors, para=para[1:nfactors])
  factors = spca_res$loadings[,1:nfactors]
  factors = factors[, complete.cases(t(factors))]  # remove NA cols
  loadings = x_betas %*% factors
  if(!full) {
    return(loadings)
  } else {
    return(list('L' = loadings, 'F' = factors))
  }
}


# generates summary statistics by regressing X and Y on each variable in G
get_betas_stderrs = function(G, X, Y, G2=c(), x_betas = c(), x_stderrs = c()) {
  num_snps = dim(G)[2]
  num_exposures = dim(X)[2]
  # regress to get betas and stderrs for Y
  y_betas = integer(num_snps)
  y_stderrs = integer(num_snps)
  xb_precomp = length(x_betas) != 0  # check if x_betas precomputed
  if(!xb_precomp) {
    x_betas = matrix(0, num_exposures, num_snps)
    x_stderrs = matrix(0, num_exposures, num_snps)
  }
  for (i in 1:num_snps) {
    Gi = G[, i]
    Gi = (Gi - mean(Gi)) / sd(Gi)  # standardize genotypes
    if(length(G2) > 0) {  # two-sample mode -- run Y regressions on different data
      G2i = G2[, i]
      G2i = (G2i - mean(G2i)) / sd(G2i)  # standardize genotypes
      res_gy = summary(lm(Y ~ G2i))
    } else {
      res_gy = summary(lm(Y ~ Gi))
    }
    y_betas[i] = res_gy$coefficients[2,1]
    y_stderrs[i] = as.numeric(coef(res_gy)[, "Std. Error"][2])
    if(!xb_precomp) {
      for (j in 1:num_exposures) {
        res_gx = summary(lm(X[, j] ~ Gi))
        x_betas[j,i] = res_gx$coefficients[2,1]
        x_stderrs[j,i] = as.numeric(coef(res_gx)[, "Std. Error"][2])
      }
    }
  }
  if(!xb_precomp) {
    x_betas = t(x_betas)
    x_stderrs = t(x_stderrs)
  }
  return(list('x_betas' = x_betas, 'x_stderrs' = x_stderrs, 'y_betas' = y_betas, 'y_stderrs' = y_stderrs))
}
