# simulations for MVMR paper with sparse (point-normal) effects
library("MASS")

# Simulate a covariance matrix
sim_sigma = function(num_vars, rho_max=0.5, rho_min=-0.3) {
  if(num_vars == 0) {
    return(diag(0))
  }
  Sigma = matrix(0, num_vars, num_vars)
  for (i in 1:num_vars) {
    for (j in 1:i) {
      if (i==j) {
        Sigma[i,j] = 1
      } else {
        Sigma[i,j] = Sigma[j,i] = runif(1, rho_min, rho_max)
      }
    }
  }
  Sigma = t(Sigma) %*% Sigma  # ensure positive semidefinite
  Sigma = cov2cor(Sigma)
  return(Sigma)
}


# Simulate point-(multivariate)normal, i.e. spike-and-slab prior
sim_point_mvn = function(num, phi, mu, Sigma) {
  res = mvrnorm(num, mu, Sigma)
  if (num == 1) {
    for (i in 1:length(res)) {
      bern = rbinom(1, 1, phi)
      res[i] = res[i] * bern
    }
  } else {
    for (i in 1:dim(res)[1]) {
      for (j in 1:dim(res)[2]) {
        bern = rbinom(1, 1, phi)
        res[i,j] = res[i,j] * bern
      }
    }
  }
  return(res)
}


# Standardize columns in a matrix
standardize = function(mat) {
  if(length(mat) <= 1) {
    return(0)
  }
  cols = dim(mat)[2]
  if (is.null(cols) || cols == 1) {
    if(var(mat) == 0) {
      mat = mat - mean(mat)
    } else {
      mat = (mat - mean(mat)) / sd(mat)
    }
  } else {
    for (i in 1:cols) {
      if(var(mat[,i]) == 0) {
        mat[,i] = (mat[,i] - mean(mat[,i]))
      } else {
        mat[,i] = (mat[,i] - mean(mat[,i])) / sd(mat[,i])
      }
    }
  }
  return(mat)
}


# Simulate SNPs, ensuring they aren't correlated with noise term of Z
sim_genos = function(N, M, epsilon_z, freqs=c(), freqs2=c(), tol=0.005) {
  G = c()
  snps_done = 0
  if(length(freqs) == 0) {  # if allele freqs not given, sim uniform(0.05, 0.5)
    freqs = runif(M, 0.05, 0.5)
  }
  while(snps_done < M) {
    # newsnp = rnorm(N)
    newsnp = rbinom(N, 2, freqs[snps_done + 1])
    if(length(freqs2) > 0) {
      newsnp2 = rbinom(N, 2, freqs2[snps_done + 1])
      newsnp = c(newsnp, newsnp2)
    }
    maxcor = max(abs(cor(newsnp, epsilon_z)))
    while(maxcor > tol) {
      # newsnp = rnorm(N)
      newsnp = rbinom(N, 2, freqs[snps_done + 1])
      if(length(freqs2) > 0) {
        newsnp2 = rbinom(N, 2, freqs2[snps_done + 1])
        newsnp = c(newsnp, newsnp2)
      }
      maxcor = max(abs(cor(newsnp, epsilon_z)))
    }
    G = cbind(G, newsnp)
    snps_done = snps_done + 1
  }
  return(G)
}


# simulate SNPs with population stratification (balding-nichols model)
# wrapper around sim_G which provides allele frequencies
sim_genos_balding_nichols = function(N, M, epsilon_z, fst, tol=0.005) {
  # simulate allele frequencies from balding-nichols model
  mult = (1 - fst) / fst
  avg_freqs = runif(M, 0.05, 0.5)
  pop1_freqs = rbeta(M, avg_freqs * mult, (1 - avg_freqs) * mult)
  pop2_freqs = rbeta(M, avg_freqs * mult, (1 - avg_freqs) * mult)
  # now draw the actual genotypes, using the sim_G method
  half = floor(N/2)
  G = sim_genos(half, M, epsilon_z, freqs=pop1_freqs, freqs2=pop2_freqs, tol=tol)
  # record population labels and standardize to control R^2 on phenotypes
  labels = c(rep(1, half), rep(2, N - half))
  labels = (labels - mean(labels)) / sd(labels)  # standardize
  return(list('G' = G, 'labels' = labels, 'p1frq' = pop1_freqs, 'p2frq' = pop2_freqs))
}


# wrapper for simulating genotype matrix G
sim_g = function(G, N, M, epsilon, covar, fst, corr_g, tol=0.005) {
  if(length(G) == 0) {  # if we are not using a pre-built genotype matrix
    if(corr_g) {  # outprioritizes balding-nichols, which cannot ensure correlations
      G = mvrnorm(N, rep(0, M), covar$g)
      labels = rep(0, N)
    } else if(fst != 0.0) {  # simulate with population structure, ignored if corr_g
      sim_res = sim_genos_balding_nichols(N, M, epsilon$z, fst, tol=tol)
      G = sim_res$G
      labels = sim_res$labels  # population labels
    }
    else {
      G = sim_genos(N, M, epsilon$z, tol=tol)
      labels = rep(0, N)
    }
  } else {
    labels = rep(0, N)
  }
  G = standardize(G)
  return(list('G' = G, 'labels' = labels))
}


# assess whether R^2 between variables is as intended
assess_r2 = function(G, X, Y, Z, K, J) {
  allvars = c()
  allr2 = c()
  if(K > 1) {
    r2gx = summary(lm(X[,1]~G))$r.squared
    r2zx = summary(lm(X[,1]~Z))$r.squared
    varx = var(X[,1])
  } else {
    r2gx = summary(lm(X~G))$r.squared
    r2zx = summary(lm(X~Z))$r.squared
    varx = var(X)
  }
  if(J > 1) {
    r2gz = summary(lm(Z[,1]~G))$r.squared
    varz = var(Z[,1])
  } else {
    r2gz = summary(lm(Z~G))$r.squared
    varz = var(Z)
  }
  r2gy = summary(lm(Y~G))$r.squared
  r2zy = summary(lm(Y~Z))$r.squared
  r2xy = summary(lm(Y~X))$r.squared
  vars = c(varz, varx, var(Y))
  #print(vars)
  allvars = rbind(allvars, vars)
  r2s = c(r2gx, r2gy, r2gz, r2zx, r2zy, r2xy)
  allr2 = rbind(allr2, r2s)

  print('Average variances for Z, X, Y:')
  print(c(mean(allvars[,1]), mean(allvars[,2]), mean(allvars[,3])))
  print('Variance in variances for Z, X, Y:')
  print(c(var(allvars[,1]), var(allvars[,2]), var(allvars[,3])))
  print('Average R^2 for GX, GY, GZ, ZX, ZY, XY:')
  print(c(mean(allr2[,1]), mean(allr2[,2]), mean(allr2[,3]), mean(allr2[,4]), mean(allr2[,5]), mean(allr2[,6])))
  print('Variance in R^2 for GX, GY, GZ, ZX, ZY, XY:')
  print(c(var(allr2[,1]), var(allr2[,2]), var(allr2[,3]), var(allr2[,4]), var(allr2[,5]), var(allr2[,6])))
  cat('\n\n')
}


# check that input arguments are valid
check_inputs = function(N, M, K, J, H, L, beta, gamma, phi) {
  if(min(N, M, K, J, H, L) < 1) {
    stop("N, M, K, J, H, and L must be at least 1.")
  }
  if(min(unlist(gamma), unlist(phi)) < 0) {
    stop("All gamma and phi must be between 0 and 1.")
  }
  if(max(unlist(gamma), unlist(phi)) > 1) {
    stop("All gamma and phi must be between 0 and 1.")
  }
  if(gamma$gz + gamma$x1z > 1) {
    stop("gamma_gz + gamma_x1z must be less than or equal to 1.")
  }
  if(gamma$gx2 + gamma$zx + gamma$wx + gamma$yx2 > 1) {
    stop("gamma_gx2 + gamma_zx + gamma_wx + gamma_yx2 must be less than or equal to 1.")
  }
  if(gamma$gy + gamma$zy + gamma$uy + sum(beta^2) > 1) {
    stop("gamma_gy + gamma_zy + gamma_uy + sum(beta^2) must be less than or equal to 1.")
  }
  if(length(beta) > K) {
    stop("Cannot specify more betas than number of exposures K.")
  }
}


# set values of certain variables based on settings & arguments from user
set_vars = function(gamma_gx, gamma_gy, gamma_gz, gamma_zx,
                    gamma_zy, gamma_gu, gamma_uy, gamma_gw, gamma_wx,
                    phi_gx, phi_gy, phi_gz, phi_zx, phi_gu, phi_gw, phi_wx,
                    gamma_gx2, gamma_x1z, gamma_yx2, phi_gx2, phi_x1z, phi_yx2,
                    gamma_psx1, gamma_psx2, gamma_psy, gamma_psz,
                    beta, K, upstream_pct, reverse_cause, X1, X2, Z) {
  # store gamma and phi variables in lists for easier manipulation
  gamma = list('gx' = gamma_gx, 'gy' = gamma_gy, 'gz' = gamma_gz, 'gu' = gamma_gu,
               'gw' = gamma_gw, 'zx' = gamma_zx, 'zy' = gamma_zy, 'uy' = gamma_uy,
               'wx' = gamma_wx, 'gx2' = gamma_gx2, 'x1z' = gamma_x1z, 'yx2' = gamma_yx2,
               'psx1' = gamma_psx1, 'psx2' = gamma_psx2, 'psy' = gamma_psy, 'psz' = gamma_psz)
  phi = list('gx' = phi_gx, 'gy' = phi_gy, 'gz' = phi_gz, 'gu' = phi_gu,
             'gw' = phi_gw, 'zx' = phi_zx, 'wx' = phi_wx,
             'gx2' = phi_gx2, 'x1z' = phi_x1z, 'yx2' = phi_yx2)

  # if left to defaults, all X are in X2, so default theta/phi for x2 to vals for x1
  if(upstream_pct == 0.0 && gamma$gx2 == 0.0) {
    gamma$gx2 = gamma$gx
    phi$gx2 = phi$gx
  }

  # determine how many exposures will be upstream of factors
  K1 = floor(upstream_pct * K)
  K2 = K - K1
  if(K < length(beta)) {
    stop('Error: length of beta greater than K.')
  }
  if(K1 < length(beta) && reverse_cause) {
    stop(paste0("Error: More causal exposures than number of exposures ",
                "upstream of Y. Need upstream_pct * K >= length(beta)."))
  }
  if(reverse_cause) {  # fill in unspecified betas with 0
    beta = append(beta, integer(K1 - length(beta)))
  } else {
    beta = append(beta, integer(K - length(beta)))
  }

  # if X1/X2/Z pre-specified, set genetic effects to 0
  gamma$gx = ifelse(length(X1)==0, gamma$gx, 0)
  gamma$gx2 = ifelse(length(X2)==0, gamma$gx2, 0)
  gamma$gz = ifelse(length(Z)==0, gamma$gz, 0)

  return(list('gamma' = gamma, 'phi' = phi, 'beta' =  beta, 'K1' = K1, 'K2' = K2))
}


# simulate covariance variances for MVN  and point-MVN random variables
sim_covar = function(G, N, M, K, J, H, L, K1, K2, corr_g, corr_eff) {
  covar = list()
  if (length(G) == 0) {  # if we are not using a pre-built genotype matrix
    if (!corr_g) {
      covar$g = diag(M)
    } else{
      covar$g = sim_sigma(M, 0.5)
    }
  }
  if (corr_eff) {
    covar$theta_gx = sim_sigma(K1, rho_max=0.5)
    covar$theta_gz = sim_sigma(J, rho_max=0.5)
    covar$theta_gu = sim_sigma(H, rho_max=0.5)
    covar$theta_gw = sim_sigma(L, rho_max=0.5)
    covar$theta_zx = sim_sigma(K2, rho_max=0.5)
    covar$theta_wx = sim_sigma(K2, rho_max=0.5)
    covar$theta_gx2 = sim_sigma(K2, rho_max=0.5)
    covar$theta_yx2 = sim_sigma(K2, rho_max=0.5)
    covar$theta_x1z = sim_sigma(J, rho_max=0.5)
  } else {
    covar$theta_gx = diag(K1)
    covar$theta_gz = diag(J)
    covar$theta_gu = diag(H)
    covar$theta_gw = diag(L)
    covar$theta_zx = diag(K2)
    covar$theta_wx = diag(K2)
    covar$theta_gx2 = diag(K2)
    covar$theta_yx2 = diag(K2)
    covar$theta_x1z = diag(J)
  }
  covar$theta_gy = 1
  covar$theta_zy = 1
  covar$theta_uy = 1
  return(covar)
}


# compute per-SNP effect sizes, adjusted for sparsity and num of effects.
# for instance, for G-Z effects, we want the SNPs to collectively account for
#   gamma_gz percent of the variance of Z. We divide by M so this effect is
#   divided evenly among the M SNPs that might affect Z. Then we divide by
#   the density parameter phi$gz -- for example, if only 25% of the M SNPs are
#   expected to have an effect, we multiply their effect by 4. Finally,
#   we take the square root because this will multiply a covariance matrix.
adj_gamma = function(gamma, phi, M, K, J, H, L, K1) {
  adj = list()
  adj$gamma_gx = sqrt(gamma$gx / M / phi$gx)
  adj$gamma_gy = sqrt(gamma$gy / M / phi$gy)
  adj$gamma_gz = sqrt(gamma$gz / M / phi$gz)
  adj$gamma_gu = sqrt(gamma$gu / M / phi$gu)
  adj$gamma_gw = sqrt(gamma$gw / M / phi$gw)
  adj$gamma_zx = sqrt(gamma$zx / J / phi$zx)
  adj$gamma_wx = sqrt(gamma$wx / L / phi$wx)
  adj$gamma_zy = sqrt(gamma$zy / J)
  adj$gamma_uy = sqrt(gamma$uy / H)
  adj$gamma_gx2 = sqrt(gamma$gx2 / M / phi$gx2)
  adj$gamma_yx2 = sqrt(gamma$yx2 / phi$yx2)
  if(K1 > 0) {
    adj$gamma_x1z = sqrt(gamma$x1z / K1 / phi$x1z)
  } else {
    adj$gamma_x1z = 0
  }
  return(adj)
}


# simulate effect matrices
sim_effects = function(adj, phi, covar, mu, M, K, J, H, L, K1, K2) {
  theta = list()
  theta$gy = sim_point_mvn(M, phi$gy, mu, covar$theta_gy) * adj$gamma_gy
  theta$gz = sim_point_mvn(M, phi$gz, rep(mu, J), covar$theta_gz) * adj$gamma_gz
  theta$gu = sim_point_mvn(M, phi$gu, rep(mu, H), covar$theta_gu) * adj$gamma_gu
  theta$gw = sim_point_mvn(M, phi$gw, rep(mu, L), covar$theta_gw) * adj$gamma_gw
  theta$zy = mvrnorm(J, mu, covar$theta_zy) * adj$gamma_zy
  theta$uy = mvrnorm(H, mu, covar$theta_uy) * adj$gamma_uy
  if(K1 > 0) {
    theta$x1z = sim_point_mvn(K1, phi$x1z, rep(mu, J), covar$theta_x1z) * adj$gamma_x1z
    theta$gx = sim_point_mvn(M, phi$gx, rep(mu, K1), covar$theta_gx) * adj$gamma_gx
  } else {
    theta$x1z = NULL
    theta$gx = -999  # prevents "theta$gx" from grabbing theta$gx2
  }
  if(K2 > 0) {
    theta$zx = sim_point_mvn(J, phi$zx, rep(mu, K2), covar$theta_zx) * adj$gamma_zx
    theta$wx = sim_point_mvn(L, phi$wx, rep(mu, K2), covar$theta_wx) * adj$gamma_wx
    theta$gx2 = sim_point_mvn(M, phi$gx2, rep(mu, K2), covar$theta_gx2) * adj$gamma_gx2
    theta$yx2 = sim_point_mvn(1, phi$yx2, rep(mu, K2), covar$theta_yx2) * adj$gamma_yx2
  } else {
    theta$zx = NULL
    theta$wx = NULL
    theta$gx2 = NULL
    theta$yx2 = NULL
  }
  return(theta)
}


# simulate noise terms (i.e. those not from genetic/phenotypic effects)
sim_noise = function(gamma, beta, X1, X2, Z, N, M, K, J, H, L, K1, K2) {
  epsilon = list()
  # If X1/X2/Z pre-specified, then the existing values are the "noise matrices"
  #   since they are not simulated by this script
  noise_var_x1 = sqrt(1 - gamma$gx)
  noise_var_x2 = sqrt(1 - gamma$gx2 - gamma$zx - gamma$wx - gamma$yx2)
  noise_var_y = sqrt(1 - sum(beta^2) - gamma$gy - gamma$zy - gamma$uy)
  noise_var_z = sqrt(1 - gamma$gz - gamma$x1z)
  noise_var_u = sqrt(1 - gamma$gu)
  noise_var_w = sqrt(1 - gamma$gw)
  if(K1 > 0) {
    if(length(X1) == 0) {
      epsilon$x1 = mvrnorm(N, rep(0, K1), diag(K1)) * noise_var_x1
    } else {
      epsilon$x1 = standardize(X1) * noise_var_x1
    }
  } else {
    epsilon$x1 = NULL
  }
  if(K2 > 0) {
    if(length(X2) == 0) {
      epsilon$x2 = mvrnorm(N, rep(0, K2), diag(K2)) * noise_var_x2
    } else {
      epsilon$x2 = standardize(X2) * noise_var_x2
    }
  } else {
    epsilon$x2 = NULL
  }
  epsilon$y = matrix(rnorm(N, 0, 1) * noise_var_y, N, 1)
  if(length(Z) == 0) {
    epsilon$z = mvrnorm(N, rep(0, J), diag(J)) * noise_var_z
  } else {
    epsilon$z = standardize(Z) * noise_var_z
  }
  epsilon$u = mvrnorm(N, rep(0, H), diag(H)) * noise_var_u
  epsilon$w = mvrnorm(N, rep(0, L), diag(L)) * noise_var_w
  return(epsilon)
}


# simulate phenotypes from the structural equation model
sim_phenos = function(G, labels, theta, gamma, epsilon, beta, K1, K2, reverse_cause) {
  phe = list()
  if(K1 > 0) {
    X1tilde = G %*% theta$gx + labels * sqrt(gamma$psx1)
    X1 = G %*% theta$gx + labels * sqrt(gamma$psx1) + epsilon$x1
    Ztilde = G %*% theta$gz + X1tilde %*% theta$x1z + labels * sqrt(gamma$psz)
    Z = G %*% theta$gz + X1 %*% theta$x1z + labels * sqrt(gamma$psz) + epsilon$z
  } else {
    X1tilde = NULL
    X1 = NULL
    Ztilde = G %*% theta$gz + labels * sqrt(gamma$psz)
    Z = G %*% theta$gz + labels * sqrt(gamma$psz) + epsilon$z
  }
  Utilde = G %*% theta$gu
  U = G %*% theta$gu + epsilon$u
  Wtilde = G %*% theta$gw
  W = G %*% theta$gw + epsilon$w
  if(reverse_cause) {
    if(K1 > 0) {
      Ytilde = X1tilde %*% beta + G %*% theta$gy + Ztilde %*% theta$zy +
        Utilde %*% theta$uy + labels * sqrt(gamma$psy)
      Y = X1 %*% beta + G %*% theta$gy + Z %*% theta$zy +
        U %*% theta$uy + labels * sqrt(gamma$psy) + epsilon$y
    } else {
      Ytilde = G %*% theta$gy + Ztilde %*% theta$zy +
        Utilde %*% theta$uy + labels * sqrt(gamma$psy)
      Y = G %*% theta$gy + Z %*% theta$zy +
        U %*% theta$uy + labels * sqrt(gamma$psy) + epsilon$y
    }
    if(K2 > 0) {
      X2tilde = G %*% theta$gx2 + Ztilde %*% theta$zx + Wtilde %*% theta$wx +
        Ytilde %*% theta$yx2 + labels * sqrt(gamma$psx2)
      X2 = G %*% theta$gx2 + Z %*% theta$zx + W %*% theta$wx +
        Y %*% theta$yx2 + labels * sqrt(gamma$psx2) + epsilon$x2
    } else {
      X2tilde = NULL
      X2 = NULL
    }
    X = cbind(X1, X2)
    Xtilde = cbind(X1tilde, X2tilde)
  } else {
    if(K2 > 0) {
      X2tilde = G %*% theta$gx2 + Ztilde %*% theta$zx +
        Wtilde %*% theta$wx + labels * sqrt(gamma$psx2)
      X2 = G %*% theta$gx2 + Z %*% theta$zx +
        W %*% theta$wx + labels * sqrt(gamma$psx2) + epsilon$x2
    } else {
      X2tilde = NULL
      X2 = NULL
    }
    X = cbind(X1, X2)
    Xtilde = cbind(X1tilde, X2tilde)
    # # shuffle cols; lets causal Xs be from either X1 or X2
    # order = sample(ncol(X))
    # X = X[,order]
    # Xtilde = Xtilde[,order]
    Ytilde = Xtilde %*% beta + G %*% theta$gy + Ztilde %*% theta$zy +
      Utilde %*% theta$uy + labels * sqrt(gamma$psy)
    Y = X %*% beta + G %*% theta$gy + Z %*% theta$zy +
      U %*% theta$uy + labels * sqrt(gamma$psy) + epsilon$y
  }
  return(list('X' = X, 'X1' = X1, 'X2' = X2, 'Z' = Z, 'Y' = Y, 'W' = W, 'U' = U,
              'Xtilde' = Xtilde, 'X1tilde' = X1tilde, 'X2tilde' = X2tilde,
              'Ztilde' = Ztilde, 'Ytilde' = Ytilde, 'Wtilde' = Wtilde, 'Utilde' = Utilde))
}


# main user-facing wrapper method for simulating phenotypes
run_sim = function(N=10000, M=100, K=1, J=1, H=1, L=1, beta=0.2,
                   gamma_gx=0.5, gamma_gy=0.2, gamma_gz=0.2, gamma_gu=0.0, gamma_gw=0.0,
                   gamma_zx=0.2, gamma_zy=0.2, gamma_uy=0.0, gamma_wx=0.0,
                   phi_gx=0.1, phi_gy=0.1, phi_gz=0.1, phi_gu=0.1, phi_gw=0.1,
                   phi_zx=0.5, phi_wx=0.5,
                   gamma_gx2=0.0, gamma_x1z=0.0, gamma_yx2=0.0,
                   phi_gx2=0.1, phi_x1z=1.0, phi_yx2=1.0,
                   upstream_pct=0, reverse_cause=FALSE,
                   fst=0.0, gamma_psx1=0.1, gamma_psx2=0.1, gamma_psy=0.1, gamma_psz=0.1,
                   G=c(), X1=c(), X2=c(), Z=c(),
                   corr_g=FALSE, corr_eff=FALSE, mu=0, tol=0.005) {

  # set variable values based on arguments & check that they are valid
  vars = set_vars(gamma_gx, gamma_gy, gamma_gz, gamma_zx,
                  gamma_zy, gamma_gu, gamma_uy, gamma_gw, gamma_wx,
                  phi_gx, phi_gy, phi_gz, phi_zx, phi_gu, phi_gw, phi_wx,
                  gamma_gx2, gamma_x1z, gamma_yx2, phi_gx2, phi_x1z, phi_yx2,
                  gamma_psx1, gamma_psx2, gamma_psy, gamma_psz,
                  beta, K, upstream_pct, reverse_cause, X1, X2, Z)
  gamma = vars$gamma
  phi = vars$phi
  beta = vars$beta
  K1 = vars$K1
  K2 = vars$K2
  check_inputs(N, M, K, J, H, L, beta, gamma, phi)

  # simulate covariance matrices, sparsity-adjusted per-snp effect sizes,
  #   effect matrices, and noise matrices
  covar = sim_covar(G, N, M, K, J, H, L, K1, K2, corr_g, corr_eff)
  adj = adj_gamma(gamma, phi, M, K, J, H, L, K1)
  theta = sim_effects(adj, phi, covar, mu, M, K, J, H, L, K1, K2)
  theta$xy = beta
  epsilon = sim_noise(gamma, beta, X1, X2, Z, N, M, K, J, H, L, K1, K2)

  # Simulate genotypes (if appropriate) and phenotypes
  genos = sim_g(G, N, M, epsilon, covar, fst, corr_g, tol=tol)
  G = genos$G
  labels = genos$labels
  phenos = sim_phenos(G, labels, theta, gamma, epsilon, beta, K1, K2, reverse_cause)

  return(list('G' = G, 'phe' = phenos, 'theta' = theta, 'epsilon' = epsilon,
              'covar' = covar, 'adj' = adj, 'labels' = labels))
}
