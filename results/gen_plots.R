# Code for plots of simulation results
library(ggplot2)
library(ggpubr)
library(cowplot)
library(qvalue)
library(pROC)

# colorblind-friendly palette
cbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#000000")

# check if this is a pip method
pip_checker = function(method) {
  return(grepl('susie', method) || grepl('mrash', method) || grepl('zuber', method) || grepl('brms', method) || grepl('varbvs', method) || grepl('vebboost', method) || grepl('ctwas', method) || grepl('famr', method))
}

# change internal method name to displayed method name
display_name = function(method) {
  if(method == 'famr_gfa')  return('FAMR-Susie')
  if(method == 'famr_none')  return('FAMR-Susie_noFA')
  if(method == 'famr_prune0')  return('FAMR-Susie_clump0')
  if(method == 'ivw')  return('IVW')
  if(method == 'ivw_gfa')  return('IVW-GFA')
  if(method == 'median')  return('Median')
  if(method == 'median_gfa')  return('Median-GFA')
  if(method == 'robust')  return('Robust')
  if(method == 'robust_gfa')  return('Robust-GFA')
  if(method == 'grapple')  return('GRAPPLE')
  if(method == 'cml')  return('cML')
  if(method == 'bma')  return('BMA')
  return(method)  # if not one of the above, return the original method name
}

# Plot per-variable bar showing true/false positive rate
group_barchart = function(dirname, prefix='all_res_', n_x_show = 8, thresh = 0.05,
                          pip_thresh = 0.95, methods=c(), oldpip=F) {
  method_labels = c()
  beta_labels = c()
  all_pos_rates = c()

  for(m in methods) {
    fname = paste0(dirname, prefix, m, '.txt')
    res = na.omit(read.table(fname))
    num_cols = dim(res)[2]
    ss = strsplit(m, '_')[[1]][1]  # isolate the second-stage method name
    if(pip_checker(ss)) {
      num_betas = num_cols
      if(!oldpip) {
        res = res[, (num_betas/3*2+1):ncol(res)]
      }
    } else {
      num_betas = num_cols / 3
    }

    # compute rates of (true or false) positives
    for (i in 1:n_x_show) {
      if(pip_checker(ss)) {
        posrate = sum(res[,i] >= pip_thresh) / length(res[,i])
      } else {
        posrate = sum(res[,(i+num_betas*2)] <= thresh) / length(res[,(i+num_betas*2)])
      }
      method_labels = c(method_labels, m)
      beta_labels = c(beta_labels, as.character(i))
      all_pos_rates = c(all_pos_rates, posrate)
    }
  }

  data = data.frame(method_labels, beta_labels, all_pos_rates)
  ggplot(data, aes(fill=method_labels, y=all_pos_rates, x=beta_labels)) +
    geom_bar(position="dodge", stat="identity") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
    scale_fill_manual(values=cbPalette) + theme(text = element_text(size = 20)) +
    xlab('Exposures') + ylab('Positive Rate') + labs(fill = "Method") + ylim(0, 1)
}


# stacked bar plot showing true and false positives, similar to cTWAS paper
stacked_barplot = function(dirname, truevars, prefix='all_res_', thresh = 0.05,
                           pip_thresh = 0.95, methods=c(), oldpip=F, use_qval=F,
                           labels=c()) {
  method_labels = c()
  counts = c()
  tf_labels = c()  # whether true or false
  truevars = as.numeric(unlist(strsplit(truevars, ",")))

  for(m in methods) {
    if(length(labels) > 0) {
      lab = labels[which(methods == m)]
    } else {
      lab = m
    }
    fname = paste0(dirname, prefix, m, '.txt')
    res = na.omit(read.table(fname))
    ss = strsplit(m, '_')[[1]][1]  # isolate the second-stage method name

    # get columns with pips or pvalues
    if(pip_checker(ss) && oldpip) {
      colsUsed = 1:ncol(res)
    } else {
      start = ncol(res) / 3 * 2 + 1  # last third of cols are pvals
      colsUsed = start:ncol(res)
    }
    subset_res = as.matrix(res[, colsUsed])

    if(use_qval) {
      all_qvals = c()
      all_pvals = subset_res
      for(i in 1:nrow(res)) {
        qvals = qvalue(all_pvals[i,], lambda=0)$qvalues
        all_qvals = rbind(all_qvals, qvals)
      }
      subset_res = all_qvals
    }

    # count true and false positives
    true_vals = subset_res[, truevars]
    false_vals = subset_res[, -truevars]
    if(pip_checker(ss)) {
      tp = sum(true_vals >= pip_thresh)
      fp = sum(false_vals >= pip_thresh)
    } else {
      tp = sum(true_vals <= thresh)
      fp = sum(false_vals <= thresh)
    }

    method_labels = c(method_labels, lab, lab)
    counts = c(counts, tp, fp)
    tf_labels = c(tf_labels, 'True Pos.', 'False Pos.')
  }

  df = data.frame(method_labels, counts, tf_labels)
  ggplot(df, aes(fill=tf_labels, y=counts, x=method_labels)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=c("#888888", "#44AA99")) +
    theme(text = element_text(size = 20), axis.title.x=element_blank(),
          axis.title.y=element_blank(), legend.title=element_blank())
}


# bar plot showing power at a fixed FPR threshold
power_fpr_thresh_plot = function(dirname, truevars, prefix='all_res_', n_x_show = 4,
                                 thresh = 0.05, methods=c()) {
  method_labels = c()
  beta_labels = c()
  all_pos_rates = c()
  truevars = as.numeric(unlist(strsplit(truevars, ",")))

  for(m in methods) {
    fname = paste0(dirname, prefix, m, '.txt')
    res = na.omit(read.table(fname))
    num_sims = dim(res)[1]
    num_cols = dim(res)[2]
    ss = strsplit(m, '_')[[1]][1]  # isolate the second-stage method name
    if(pip_checker(ss)) {
      num_betas = num_cols
      if(!oldpip) {
        res = res[, (num_betas/3*2+1):ncol(res)]
      }
    } else {
      num_betas = num_cols / 3
    }
    num_falsevars = num_betas - length(truevars)
    tot_false = num_sims * num_falsevars
    falselim = tot_false * thresh

    # determine power of each true variable at specified FPR threshold
    sorted_pvals = c()
    for(row in 1:num_sims) {
      for(col in 1:num_betas) {
        if(pip_checker(ss)) {
          entry = res[row,col]
        } else {
          entry = res[row,(2*num_betas+col)]
        }
        sorted_pvals = rbind(sorted_pvals, c(entry, col))
      }
    }
    if(pip_checker(ss)) {
      sorted_pvals = sorted_pvals[order(sorted_pvals[,1],decreasing=TRUE),]
    } else {
      sorted_pvals = sorted_pvals[order(sorted_pvals[,1],decreasing=FALSE),]
    }

    truepos = rep(0, num_betas)
    falsepos = 0
    for(row in 1:dim(sorted_pvals)[1]) {
      var = sorted_pvals[row,2]
      if(var %in% truevars) {
        truepos[var] = truepos[var] + 1
      } else{
        falsepos = falsepos + 1
      }
      if(falsepos >= falselim) {
        break
      }
    }
    truepos = (truepos / num_sims)[1:n_x_show]
    for(i in 1:length(truepos)) {
      method_labels = c(method_labels, m)
      beta_labels = c(beta_labels, as.character(i))
      all_pos_rates = c(all_pos_rates, truepos[i])
    }
  }

  data = data.frame(method_labels, beta_labels, all_pos_rates)
  ggplot(data, aes(fill=method_labels, y=all_pos_rates, x=beta_labels)) +
    geom_bar(position="dodge", stat="identity") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
    scale_fill_manual(values=cbPalette) + theme(text = element_text(size = 20)) +
    xlab('Exposures') + ylab(paste0('Power @ FPR = ', thresh)) + labs(fill = "Method") + ylim(0, 1)
}


# line plots showing FPR/Power trending as a function of confounding/true signal strength
param_lineplot = function(dirvec, param_name, param_vals, truepos,
                          prefix='all_res_', thresh = 0.05, pip_thresh = 0.95,
                          methods=c(), type='both', use_qval=T, oldpip=F,
                          saveloc='NONE') {
  all_method_labels = c()
  all_param_vals = c()
  all_fpr = c()
  all_power = c()

  for(i in 1:length(dirvec)) {
    dirname = dirvec[i]
    thisval = param_vals[i]
    for(m in methods) {
      fpr = 0
      power = 0
      fname = paste0(dirname, prefix, m, '.txt')
      res = na.omit(read.table(fname))
      num_cols = dim(res)[2]
      if(pip_checker(m) && !oldpip) {
        res = res[, (num_cols/3*2+1):ncol(res)]
      }
      num_betas = num_cols / 3
      if(oldpip)  num_betas = num_cols

      # convert to q-values if requested for non-pip (p value-based) method
      # then take 1-qval and treat them like pips
      if(use_qval && !pip_checker(m)) {
        all_pvals = res[,(num_betas*2+1):num_cols]
        all_qvals = t(apply(all_pvals, 1, function(x) qvalue(x, lambda=0)$qvalues))
        res = 1 - all_qvals
      }

      # compute rates of (true or false) positives, store results
      tp = sum(res[,truepos] >= pip_thresh)
      fp = sum(res[,-truepos] >= pip_thresh)
      fn = sum(res[,truepos] < pip_thresh)
      fpr = fp / (tp + fp)
      power = tp / (tp + fn)
      all_method_labels = c(all_method_labels, display_name(m))
      all_param_vals = c(all_param_vals, thisval)
      all_power = c(all_power, power)
      all_fpr = c(all_fpr, fpr)
    }
  }

  data_fpr = data.frame(all_method_labels, all_param_vals, all_fpr)
  data_power = data.frame(all_method_labels, all_param_vals, all_power)
  fpr_plot = ggplot(data_fpr,
                    aes(group=all_method_labels, color=all_method_labels, y=all_fpr, x=all_param_vals)) +
    geom_line(linewidth=2) + geom_point(size=4) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
    scale_fill_manual(values=cbPalette) + theme(text = element_text(size = 16)) +
    xlab(param_name) + ylab('False Positive Rate') + labs(fill = "Method") + ylim(0, 1)
  power_plot = ggplot(data_power,
                     aes(group=all_method_labels, color=all_method_labels, y=all_power, x=all_param_vals)) +
   geom_line(linewidth=2) + geom_point(size=4) +
   scale_fill_manual(values=cbPalette) + theme(text = element_text(size = 16)) +
   xlab(param_name) + ylab('Power') + labs(fill = "Method") + ylim(0, 1)
  fpr_plot$labels$colour <- " "
  power_plot$labels$colour <- " "
  if(type == 'fpr') {
    ggarrange(fpr_plot, nrow=1, ncol=1, common.legend = TRUE, legend="bottom")
  } else if(type == 'power') {
    ggarrange(power_plot, nrow=1, ncol=1, common.legend = TRUE, legend="bottom")
  } else {
    ggarrange(fpr_plot, power_plot, nrow=2, ncol=1, common.legend = TRUE, legend="bottom")
  }
  
  if(saveloc != 'NONE')  ggsave(saveloc, dpi=300, bg='white')
}

# plot showing pip calibration for a single pip-based method
pip_calibration_plot = function(res, truecols, method_name, oldpip=F, saveloc='NONE') {
  if(oldpip) {
    pip_res = res
  } else {
    num_betas = ncol(res) / 3
    pip_res = res[, (num_betas*2+1):ncol(res)]
  }

  truevars = as.numeric(unlist(strsplit(truecols, ",")))
  true_pips = as.matrix(pip_res[, truevars])
  false_pips = as.matrix(pip_res[, -truevars])
  n_bins = 10
  bin_min = 0:(n_bins-1) / n_bins
  bin_max = 1:n_bins / n_bins
  bin_mid = (bin_min + bin_max) / 2
  # pct correct, stderr, lower & upper CI limits, and num. samples for each bin
  bin_pcts = rep(0, n_bins)
  bin_se = rep(0, n_bins)
  bin_lower = rep(0, n_bins)
  bin_upper = rep(1, n_bins)
  bin_n_samp = rep(0, n_bins)

  for(i in 1:n_bins) {
    rowfp = as.numeric(apply(false_pips, 1, function(x) sum(x >= bin_min[i] & x <= bin_max[i])))
    rowtp = as.numeric(apply(true_pips, 1, function(x) sum(x >= bin_min[i] & x <= bin_max[i])))
    rowtotals = rowtp + rowfp
    rowpcts = rowtp[rowtotals > 0] / rowtotals[rowtotals > 0]
    if(length(rowpcts) > 0) {
      bin_pcts[i] = mean(rowpcts)
      bin_n_samp[i] = sum(rowtotals) #length(rowpcts)
      if(length(rowpcts) > 1) {
        bin_se[i] = sd(rowpcts) / sqrt(bin_n_samp[i])
        bin_lower[i] = bin_pcts[i] - bin_se[i]
        bin_upper[i] = bin_pcts[i] + bin_se[i]
      }
    }
  }

  df = data.frame(bin_mid = bin_mid, bin_pcts = bin_pcts, bin_n = bin_n_samp)
  ggplot(df, aes(x=bin_mid, y=bin_pcts)) +
    geom_errorbar(aes(ymin=bin_lower, ymax=bin_upper),
                  colour="black", size = 0.5, width=.01) +
    geom_point(size=1.5, shape=21, fill="#002b36") + # 21 is filled circle
    geom_text(aes(label = paste0('n = ',bin_n)), vjust = -0.5, hjust = 0.5, size = 3) +
    xlab("Expected") + ylab("Observed") +
    coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
    geom_abline(slope=1,intercept=0,colour='red', size=0.2) +
    ggtitle(paste0(method_name, ' calibration')) +
    expand_limits(y=0) + theme_cowplot() +
    theme(panel.grid.major = element_line(colour = "grey",size=0.2,linetype="dashed"),
          plot.title = element_text(size=20)) +
    guides(fill = guide_legend(title = "Sample Size")) +
    scale_x_continuous(breaks = df$bin_mid)
  
  if(saveloc != 'NONE')  ggsave(saveloc, dpi=300, bg='white')
}


# stacked bar plot showing true and false positives, similar to cTWAS paper
auc_roc_plot = function(dirname, truevars, prefix='all_res_', thresh = 0.05,
                           pip_thresh = 0.95, methods=c(), oldpip=F, use_qval=F,
                           labels=c()) {
  method_labels = c()
  counts = c()
  tf_labels = c()  # whether true or false
  all_roc = list()
  truevars = as.numeric(unlist(strsplit(truevars, ",")))
  
  for(m in methods) {
    if(length(labels) > 0) {
      lab = labels[which(methods == m)]
    } else {
      lab = m
    }
    fname = paste0(dirname, prefix, m, '.txt')
    res = na.omit(read.table(fname))
    ss = strsplit(m, '_')[[1]][1]  # isolate the second-stage method name
    
    # get columns with pips or pvalues
    if(pip_checker(ss) && oldpip) {
      colsUsed = 1:ncol(res)
    } else {
      start = ncol(res) / 3 * 2 + 1  # last third of cols are pvals
      colsUsed = start:ncol(res)
    }
    subset_res = as.matrix(res[, colsUsed])
    
    if(use_qval) {
      all_qvals = c()
      all_pvals = subset_res
      for(i in 1:nrow(res)) {
        qvals = qvalue(all_pvals[i,], lambda=0)$qvalues
        all_qvals = rbind(all_qvals, qvals)
      }
      subset_res = all_qvals
    }
    
    true_labels = subset_res * 0
    true_labels[, truevars] = 1
    true_labels = as.numeric(true_labels)
    subset_res_vec = as.numeric(subset_res)
    roc_curve <- roc(true_labels, subset_res_vec)
    all_roc[[m]] = roc_curve
  }
  
  color_vec = c('blue', 'red', 'green', 'gold', 'purple', 'brown', 'orange', 'black')
  for(idx in 1:length(methods)) {
    m = methods[[idx]]
    plt = plot(all_roc[[m]], col = color_vec[idx], main = "ROC Curve", 
               print.auc = TRUE, print.auc.y = 0.05*(10-idx), add=(idx>1))
  }
  legend(x = "top",inset = 0, legend = methods, 
         col=color_vec, lwd=7, cex=.7, horiz = TRUE)
}

# -----------------------

# Basic FPR (scaling theta_zy) and power (scaling beta) plots

this_dirvec = c('results/results_lineplot_comprehensive/array_zy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.025_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.075_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.1_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability mediated \n by confounders',
               c(0.0, 0.025, 0.05, 0.075, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_zy.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.0,0.0,0.0,0.0_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.1,0.1,0.1,0.1_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.2,0.2,0.2,0.2_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.3,0.3,0.3,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, 'Effect size of causal exposures',
               c(0.0, 0.1, 0.2, 0.3), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_beta.png')

# Plots showing effects of uncorrelated pleiotropy on both FPR and Power

this_dirvec = c('results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.05_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.1_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability due to \n direct effects',
               c(0.0, 0.05, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_gy.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.05_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.1_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability due to \n direct effects',
               c(0.0, 0.05, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_gy.png')

# Plots showing effect of number of exposures on FPR and power

this_dirvec = c('results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_10_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_60_prune_max/')
param_lineplot(this_dirvec, 'Number of exposures',
               c(10, 30, 60), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_nexpo.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_10_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_60_prune_max/')
param_lineplot(this_dirvec, 'Number of exposures',
               c(10, 30, 60), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_nexpo.png')

# Plots showing effect of prune percentage on performance

this_dirvec = c('results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.0/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.1/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.3/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.5/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.75/')
param_lineplot(this_dirvec, 'LD Clumping threshold',
               c(0.0, 0.1, 0.3, 0.5, 0.75), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_prune.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.0/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.1/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.3/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.5/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.75/')
param_lineplot(this_dirvec, 'LD Clumping threshold',
               c(0.0, 0.1, 0.3, 0.5, 0.75), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('famr_gfa', 'ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_prune.png')

# Same as above plots, but withour FAMR-Susie

this_dirvec = c('results/results_lineplot_comprehensive/array_zy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.025_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.075_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.1_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability mediated \n by confounders',
               c(0.0, 0.025, 0.05, 0.075, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_zy_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.0,0.0,0.0,0.0_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.1,0.1,0.1,0.1_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.2,0.2,0.2,0.2_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_beta_psiy_0.05_beta_0.3,0.3,0.3,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, 'Effect size of causal exposures',
               c(0.0, 0.1, 0.2, 0.3), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_beta_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.05_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.1_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability due to \n direct effects',
               c(0.0, 0.05, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_gy_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.05_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.1_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability due to \n direct effects',
               c(0.0, 0.05, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_gy_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_10_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_60_prune_max/')
param_lineplot(this_dirvec, 'Number of exposures',
               c(10, 30, 60), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_nexpo_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_10_prune_max/',
                'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_lineplot_comprehensive/array_nexpo_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_60_prune_max/')
param_lineplot(this_dirvec, 'Number of exposures',
               c(10, 30, 60), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_nexpo_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.0/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.1/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.3/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.5/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.75/')
param_lineplot(this_dirvec, 'LD Clumping threshold',
               c(0.0, 0.1, 0.3, 0.5, 0.75), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_fpr_array_prune_no_famr.png')

this_dirvec = c('results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.0/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.1/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.3/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.5/',
                'results/results_lineplot_comprehensive/array_prune_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_0.75/')
param_lineplot(this_dirvec, 'LD Clumping threshold',
               c(0.0, 0.1, 0.3, 0.5, 0.75), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('ivw', 'ivw_gfa', 'median', 'median_gfa'),
               saveloc='plot_images/lineplot_power_array_prune_no_famr.png')


# PIP calibration plots for FAMR-Susie under correlated and uncorrelated pleiotropy

res = na.omit(read.table(
  'results/results_lineplot_comprehensive/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/all_res_famr_gfa.txt'))
pip_calibration_plot(res, '1,2,3,4', 'FAMR-Susie', saveloc='plot_images/pip_calibration_famr_susie_cor_plei.png')

res = na.omit(read.table(
  'results/results_lineplot_comprehensive/array_gy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.05_nexpo_30_prune_max/all_res_famr_gfa.txt'))
pip_calibration_plot(res, '1,2,3,4', 'FAMR-Susie', saveloc='plot_images/pip_calibration_famr_susie_uncor_plei.png')


# Plots comparing FAMR versions with different parts missing (max prune, no FA)
# Prune is a misnomer, LD clumping is what is actually performed

this_dirvec = c('results/results_compare_famr_versions/array_zy_psiy_0.0_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_zy_psiy_0.025_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_zy_psiy_0.05_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_zy_psiy_0.075_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_zy_psiy_0.1_beta_0.05,0.1,0.2,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, '% heritability mediated \n by confounders',
               c(0.0, 0.025, 0.05, 0.075, 0.1), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='fpr',
               methods=c('famr_susie', 'famr_prune0', 'famr_none'),
               saveloc='plot_images/lineplot_fpr_array_zy_famr_versions.png')

this_dirvec = c('results/results_compare_famr_versions/array_beta_psiy_0.05_beta_0.0,0.0,0.0,0.0_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_beta_psiy_0.05_beta_0.1,0.1,0.1,0.1_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_beta_psiy_0.05_beta_0.2,0.2,0.2,0.2_gy_0.0_nexpo_30_prune_max/',
                'results/results_compare_famr_versions/array_beta_psiy_0.05_beta_0.3,0.3,0.3,0.3_gy_0.0_nexpo_30_prune_max/')
param_lineplot(this_dirvec, 'Effect size of causal exposures',
               c(0.0, 0.1, 0.2, 0.3), c(1,2,3,4),
               prefix='all_res_', thresh = 0.05, pip_thresh=0.95, type='power',
               methods=c('famr_susie', 'famr_prune0', 'famr_none'),
               saveloc='plot_images/lineplot_power_array_zy_famr_versions.png')
