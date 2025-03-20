# runs various methods on simulated MR data
source('utils.R')
library(cate)
library(susieR)
library(MendelianRandomization)
library(varbvs)
library(mr.ash.alpha)
library(VEB.Boost)
library(brms)
library(GRAPPLE)
source('mr-bma.R')


# calculate local false sign rate; used by mr.ash
calc_lfsr = function(beta, full.post) {
  prob0 = full.post$phi[,1]
  gridlen = dim(full.post$phi)[2]
  neg = beta < 0
  prob_fs = pnorm(0, mean = full.post$m[,2:gridlen], sd = sqrt(full.post$s2[,2:gridlen]))
  prob_fs[neg, ] = 1 - prob_fs[neg, ]
  prob_opp = rowSums(prob_fs * full.post$phi[,2:gridlen])
  lfsr = as.numeric(prob0 + prob_opp)
  return(lfsr)
}


# Run 2SLS, optionally with covariates to include
run_2sls = function(G, X, Y, covariates=c(), G2=c(), X2=c(), x_pred=c()) {
  ind = dim(X)[2] + 1
  # two stages of regression
  if(length(x_pred) == 0) {
    if(length(G2) == 0) {  # one sample setting
      x_pred = as.matrix(lm(X ~ G)$fitted.values)
    } else {  # two/three sample setting
      theta_gx_hat = as.matrix(lm(X ~ G)$coefficients[2:(dim(G)[2]+1),])
      x_pred = G2 %*% theta_gx_hat
      X = X2
    }
  }
  if(length(covariates) == 0) {
    tsls_reg = lm(Y ~ x_pred)
  } else {
    tsls_reg = lm(Y ~ x_pred + covariates)
  }

  # 2SLS stderr must be computed based on X, not predicted X. we recompute below.
  betas = as.numeric(tsls_reg$coefficients[2:ind])
  if(length(covariates) == 0) {
    Yhat = X %*% betas + tsls_reg$coefficients[1]
  } else {
    cov_ind = ind + dim(covariates)[2]
    cov_betas = as.numeric(tsls_reg$coefficients[(ind+1):cov_ind])
    Yhat = X %*% betas + covariates %*% cov_betas + tsls_reg$coefficients[1]
  }
  resid = Y - Yhat
  rss = sum(resid^2)
  sigma = sqrt(rss / tsls_reg$df.residual)
  stderrs = as.numeric(sqrt(diag(chol2inv(tsls_reg$qr$qr) * sigma^2)))[2:ind]
  stats = abs(betas / stderrs)
  pvals = lapply(stats, getp, df=tsls_reg$df.residual)

  # write results
  res = as.numeric(c(betas, stderrs, pvals))
  return(res)
}


# Run allele score MR approach -- essentially impute the predicted exposure
#   using the allele score built w/ the given weights, then run 2sls using those
run_allele_score = function(G, X, Y, weights, covariates=c(), G2=c(), X2=c()) {
  if(length(weights) == 0) {  # default to unweighted allele score
    weights = rep(1, ncol(G))
  }
  if(length(G2) > 0) {
    al_score = G2 %*% weights
    al_pred = sapply(1:ncol(X2), function(idx) as.matrix(lm(X2[,idx] ~ al_score[,idx])$fitted.values))
  } else {
    al_score = G %*% weights
    al_pred = sapply(1:ncol(X), function(idx) as.matrix(lm(X[,idx] ~ al_score[,idx])$fitted.values))
  }
  return(run_2sls(G, X, Y, covariates, G2, X2, x_pred=al_pred))
}


run_mr = function(method, x_betas, x_stderrs, y_betas, y_stderrs, corg, indices=c()) {
  num_exp = dim(x_betas)[2]
  if(length(indices) == 0) {  # specifies which indices to return values for, i.e. exposures (not covariates)
    indices = 1:num_exp  # if none specified, assume no covariates
  }

  # run MR analysis and write results
  if (method == 'egger') {
    if(num_exp == 1) {
      mr_res = mr_egger(mr_input(bx = as.numeric(x_betas), bxse = as.numeric(x_stderrs),
                                 by = as.numeric(y_betas), byse = as.numeric(y_stderrs), correlation=corg),
                        distribution = 't-dist', correl=TRUE)
    } else {
      mr_res = mr_mvegger(mr_mvinput(bx = x_betas, bxse = x_stderrs,
                                     by = y_betas, byse = y_stderrs, correlation=corg),
                          distribution = 't-dist', correl=TRUE)
    }
  } else if (method == 'median') {
    if(num_exp == 1) {
      mr_res = mr_median(mr_input(bx = as.numeric(x_betas), bxse = as.numeric(x_stderrs),
                                  by = as.numeric(y_betas), byse = as.numeric(y_stderrs), correlation=corg),
                         distribution = 't-dist')
    } else {
      mr_res = mr_mvmedian(mr_mvinput(bx = x_betas, bxse = x_stderrs,
                                      by = y_betas, byse = y_stderrs, correlation=corg),
                           distribution = 't-dist')
    }
  } else if (method == 'lasso') {
    if(num_exp == 1) {
      mr_res = mr_lasso(mr_input(bx = as.numeric(x_betas), bxse = as.numeric(x_stderrs),
                                 by = as.numeric(y_betas), byse = as.numeric(y_stderrs), correlation=corg),
                        distribution = 't-dist')
    } else {
      mr_res = mr_mvlasso(mr_mvinput(bx = x_betas, bxse = x_stderrs,
                                     by = y_betas, byse = y_stderrs, correlation=corg),
                          distribution = 't-dist')
    }
  } else if (method == 'robust') {
    if(num_exp == 1) {
      mr_res = mr_ivw(mr_input(bx = as.numeric(x_betas), bxse = as.numeric(x_stderrs),
                               by = as.numeric(y_betas), byse = as.numeric(y_stderrs)),
                      distribution = 't-dist', robust=TRUE)
    } else {
      mr_res = mr_mvivw(mr_mvinput(bx = x_betas, bxse = x_stderrs,
                                   by = y_betas, byse = y_stderrs),
                        distribution = 't-dist', robust=TRUE)
    }
  } else {  # IVW mode
    if(num_exp == 1) {
      mr_res = mr_ivw(mr_input(bx = as.numeric(x_betas), bxse = as.numeric(x_stderrs),
                               by = as.numeric(y_betas), byse = as.numeric(y_stderrs), correlation=corg),
                      distribution = 't-dist', robust=FALSE, correl=TRUE)
    } else {
      mr_res = mr_mvivw(mr_mvinput(bx = x_betas, bxse = x_stderrs,
                                   by = y_betas, byse = y_stderrs, correlation=corg),
                        distribution = 't-dist', robust=FALSE, correl=TRUE)
    }
  }
  if (method == 'egger') {  # egger has slightly different output value names
    mr_vals = as.numeric(c(mr_res$Estimate[indices], mr_res$StdError.Est[indices], mr_res$Pvalue.Est[indices]))
  } else{
    mr_vals = as.numeric(c(mr_res$Estimate[indices], mr_res$StdError[indices], mr_res$Pvalue[indices]))
  }
  return(mr_vals)
}

run_grapple = function(G, X, Y, indices, x_betas=c(), x_stderrs=c(), y_betas=c(), y_stderrs=c()) {
  if(length(x_betas) == 0) {  # regression results not provided
    regr_res = get_betas_stderrs(G, X, Y)
    x_betas = regr_res$x_betas
    x_stderrs = regr_res$x_stderrs
    y_betas = regr_res$y_betas
    y_stderrs = regr_res$y_stderrs
  }
  num_snps = dim(x_betas)[1]
  num_exp = dim(x_betas)[2]
  snp_names = 1:num_snps  # fake SNP names (these don't matter for our purposes)
  sel_p = rep(0.01, num_snps)  # flat prior on SNP selection for now

  # build dataframe of regression results needed by GRAPPLE
  df = data.frame('SNP' = snp_names, 'effect_allele' = snp_names, 'other_allele' = snp_names,
                  'gamma_out1' = y_betas, 'se_out1' = y_stderrs, 'selection_pvals' = sel_p)
  for(i in 1:num_exp) {
    name1 = paste0('gamma_exp', i)
    name2 = paste0('se_exp', i)
    df[,paste0(name1)] = as.numeric(x_betas[,i])
    df[,paste0(name2)] = as.numeric(x_stderrs[,i])
  }

  # run GRAPPLE, then collect and return results
  gres = grappleRobustEst(df, p.thres = 0.05)
  res = as.numeric(c(gres$beta.hat[indices], sqrt(diag(gres$beta.var))[indices], gres$beta.p.value[indices]))
  return(res)
}

run_susie_mrash = function(X, Y, G = c(), G2=c(), ztildehat=c(), x_pred=c(),
                           method='susie', L=0, indices=c(),
                           prior_weights=c(), prior_sigma2=c()) {
  # set up regression variables depending on what is provided
  if(length(indices) == 0) {
    indices = 1:dim(X)[2]
  }
  if(L == 0) {
    L = length(indices)
  }
  if(length(G) == 0) {  # summary statistics mode
    x_pred = X
  } else if(length(x_pred) == 0) {
    if(length(G2) == 0) {  # one sample setting
      x_pred = as.matrix(lm(X ~ G)$fitted.values)
    } else {  # two/three sample setting
      theta_gx_hat = as.matrix(lm(X ~ G)$coefficients[2:(dim(G)[2]+1),])
      x_pred = G2 %*% theta_gx_hat
      G = G2
    }
  }
  aug = cbind(x_pred, ztildehat, G)  # unaffected if ztildehat or G not provided

  # now run susie or mr.ash
  if(method == 'susie') {
    if(length(prior_weights) == 0) {
      res = susieR::susie(aug, Y, L=L)
    } else {
      if(length(prior_sigma2) == 0 || length(prior_sigma2) > length(prior_weights)) {
        # second condition because susie has some bug if L > num predictors
        res = susieR::susie(aug, Y, L=L, prior_weights=prior_weights)
      } else {
        scaled_sigma2 = prior_sigma2 / var(Y)
        res = susieR::susie(aug, Y, L=L, prior_weights=prior_weights,
                            scaled_prior_variance=scaled_sigma2,
                            estimate_prior_variance=F)
      }
    }
    pips = as.numeric(res$pip)[indices]
    post_means = as.numeric(colSums(res$alpha * res$mu))[indices]
  } else if(method == 'mrash' || method == 'mr.ash') {  # run mr.ash
    res.mr.ash = mr.ash(aug, Y)
    full.post = get.full.posterior(res.mr.ash)
    # pips = as.numeric(1 - full.post$phi[,1])
    post_means = as.numeric(res.mr.ash$beta)[indices]
    pips = (1 - calc_lfsr(res.mr.ash$beta, full.post))[indices]
  } else {  # run varbvs
    if(length(prior_weights) == 0) {
      fit.varbvs = varbvs(X=aug, y=Y, Z=NULL, family="gaussian")
    } else {
      logodds = log10(prior_weights / (1 - prior_weights))
      logodds_l = matrix(logodds, length(logodds), L)
      if(length(prior_sigma2) == 0) {
        prior_sigma2 = rep(0.01, L)  # decent default value
      }
      fit.varbvs = varbvs(X=aug, y=Y, Z=NULL, family="gaussian",
                          logodds = logodds_l, sa = prior_sigma2)
    }
    pips = as.numeric(fit.varbvs$pip)[indices]
    post_means = as.numeric(fit.varbvs$beta)[indices]
  }
  post_stderrs = post_means * 0 + 1
  res_vec = c(post_means, post_stderrs, pips)
  return(res_vec)
}


run_vebboost = function(G, X, Y, ztildehat=c(), G2=c(), x_pred=c()) {
  num_exposures = dim(X)[2]
  if(length(x_pred) == 0) {
    if(length(G2) == 0) {  # one sample setting
      x_pred = as.matrix(lm(X ~ G)$fitted.values)
    } else {  # two/three sample setting
      theta_gx_hat = as.matrix(lm(X ~ G)$coefficients[2:(dim(G)[2]+1),])
      x_pred = G2 %*% theta_gx_hat
      G = G2
    }
  }
  learner1 = makeMrAshLearner(x_pred, growMode="NA")
  learner2 = makeMrAshLearner(G, growMode="NA")
  if(length(ztildehat) == 0) {
    vfit = veb_boost(list(learner1,learner2), as.numeric(Y))
    xy.post = vfit$leaves[[1]]$learner$currentFit$mr.ash.post
  } else {  # fit effects of ztildehat first --> X-Y effect is 2nd leaf
    learner0 = makeMrAshLearner(ztildehat,growMode="NA")
    vfit = veb_boost(list(learner0,learner1,learner2), as.numeric(Y))
    xy.post = vfit$leaves[[2]]$learner$currentFit$mr.ash.post
  }
  # pips = as.numeric(1 - xy.post$phi[,1])[1:num_exposures]
  betas = as.numeric(rowSums(xy.post$phi * xy.post$m))
  pips = (1 - calc_lfsr(betas, xy.post))[1:num_exposures]
  betas = betas[1:num_exposures]
  stderrs = betas * 0 + 1
  res_vec = c(betas, stderrs, pips)
  return(res_vec)
}

run_brms = function(x_betas, y_betas, indices, prior='') {
  df = data.frame(yb = y_betas, xb = x_betas)
  varnames = names(df)[2:length(names(df))]
  formstr = paste0('yb ~ ', paste(varnames, collapse=' + '))
  if(prior == 'horseshoe') {
    fit = brm(formula = formstr, data=df, family=student(),
              prior=c(set_prior("horseshoe()", class = "b")),
              warmup=1000, iter=2000, chains=4)
  } else {
    fit = brm(formula = formstr, data=df, family=student(),
              prior=c(set_prior("normal(0, 0.25)", class = "b")),
              warmup=1000, iter=2000, chains=4)
  }

  effects = as.numeric(sapply(varnames, function(x)
    hypothesis(fit, paste0(x," > 0"))$hypothesis$Estimate))[indices]
  stderrs = as.numeric(sapply(varnames,
    function(x) hypothesis(fit, paste0(x," > 0"))$hypothesis$Est.Error))[indices]
  post.probs = as.numeric(sapply(varnames, function(x)
    max(hypothesis(fit, paste0(x," > 0"))$hypothesis$Post.Prob,
        hypothesis(fit, paste0(x," < 0"))$hypothesis$Post.Prob)))[indices]
  return(c(effects, stderrs, post.probs))
}

run_cml = function(x_betas, x_stderrs, y_betas, y_stderrs, rho_mat, N, indices, dp=F) {
  Sig_inv_l = invcov_mvmr(se_bx=x_stderrs, se_by=y_stderrs, rho_mat = rho_mat)
  
  if(!dp) {
    # regular (non-DP) MVMR-cML
    MVcML_res = MVmr_cML(b_exp=x_betas,
                         b_out=as.matrix(y_betas),
                         se_bx=x_stderrs,
                         Sig_inv_l=Sig_inv_l, n = N, 
                         K_vec = 0:20
    )
    MVcML_BIC_SE = MVcML_SdTheta(b_exp=x_betas,
                                 b_out=as.matrix(y_betas),
                                 Sig_inv_l=Sig_inv_l,
                                 theta=MVcML_res$BIC_theta,
                                 zero_ind = setdiff(1:length(y_betas),MVcML_res$BIC_invalid))
    
    effects = as.numeric(MVcML_res$BIC_theta)
    stderrs = MVcML_BIC_SE
    
    pvals = MVcMLBIC_pval = pnorm(-abs(MVcML_res$BIC_theta/MVcML_BIC_SE))*2
  } else {
    # DP version
    MVcML_res = MVmr_cML_DP(b_exp=x_betas,
                            b_out=as.matrix(y_betas),
                            se_bx=x_stderrs,
                            Sig_inv_l=Sig_inv_l, n = N, num_pert = 100,
                            K_vec = 0:20
    )
    effects = as.numeric(MVcML_res$BIC_DP_theta)
    stderrs = as.numeric(MVcML_res$BIC_DP_se)
    pvals = as.numeric(pnorm(-abs(MVcML_res$BIC_DP_theta/MVcML_res$BIC_DP_se))*2)
  }

  return(c(effects[indices], stderrs[indices], pvals[indices]))
}

run_bma = function(G, X, Y, indices, x_betas=c(), x_stderrs=c(), y_betas=c(), y_stderrs=c()) {
  if(length(x_betas) == 0) {  # regression results not provided
    regr_res = get_betas_stderrs(G, X, Y)
    x_betas = regr_res$x_betas
    x_stderrs = regr_res$x_stderrs
    y_betas = regr_res$y_betas
    y_stderrs = regr_res$y_stderrs
  }
  num_exposures = dim(x_betas)[2]
  exposures = rep('rs', num_exposures)  # set up variable "names" so we know which were picked
  for (i in 1:length(exposures)) { exposures[i] = paste0(exposures[i], i) }
  
  # run Zuber et al method and report best model and PIPs
  zuber_in = new('mvMRInput', betaX = x_betas, betaY = as.matrix(y_betas), 
                 exposure=exposures, outcome='y')
  BMA_output = summarymvMR_SSS(zuber_in, kmin=1, kmax=6, 
                               prior_prob=1/num_exposures, max_iter=10000)
  pips = as.numeric(BMA_output@pp_marginal)
  betas = as.numeric(BMA_output@BMAve_Estimate)
  stderrs = betas * 0 + 1
  res_vec = c(betas[indices], stderrs[indices], pips[indices])
  return(res_vec)
}

