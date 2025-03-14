# evaluate all simulation results in a given directory
library(optparse)
library(qvalue)

### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--indir"), type="character", default='./',
              help="Input directory with simulation results", metavar="character"),
  make_option(c("--out"), type="character", default="DEFAULT",
              help="Output results files base name", metavar="character"),
  make_option(c("--truecols"), type="character", default="",
              help="Which columns are true effect vars", metavar="character"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Write more verbose performance results", metavar="logical")
);


eval_sim_nonpip = function(res, outfile, truecols='', verbose=FALSE) {
  num_rows = dim(res)[1]
  num_cols = dim(res)[2]
  num_betas = num_cols / 3
  thresh = 0.05
  qval_thresh = 0.05

  # compute qvalues
  all_qvals = c()
  all_pvals = res[,(num_betas*2+1):num_cols]
  for(i in 1:num_rows) {
    qvals = qvalue(all_pvals[i,], lambda=0)$qvalues
    all_qvals = rbind(all_qvals, qvals)
  }

  # compute average betas, average std errors, and positive rates
  avg_betas = integer(num_betas)
  sd_betas = integer(num_betas)
  avg_se = integer(num_betas)
  positive_rates = integer(num_betas)
  qval_pos_rates = integer(num_betas)
  for (i in 1:num_betas) {
    avg_betas[i] = mean(res[,i])
    sd_betas[i] = sd(res[,i])
    avg_se[i] = mean(res[,(i+num_betas)])
    positive_rates[i] = sum(all_pvals[,i] <= thresh) / num_rows
    qval_pos_rates[i] = sum(all_qvals[,i] <= qval_thresh) / num_rows
  }

  # print results
  write('Average Estimates for Betas: ', file=outfile)
  write(avg_betas, file=outfile, append=TRUE)
  write('Positive rate for exposures: ', file=outfile, append=TRUE)
  write(positive_rates, file=outfile, append=TRUE)
  write('Standard Deviations of Betas: ', file=outfile, append=TRUE)
  write(sd_betas, file=outfile, append=TRUE)
  write('Average StdErrors for Betas: ', file=outfile, append=TRUE)
  write(avg_se, file=outfile, append=TRUE)
  write('Q-value positive rate for exposures: ', file=outfile, append=TRUE)
  write(qval_pos_rates, file=outfile, append=TRUE)

  # if given, calculate precision/recall/f1 based on provided true positive cols
  if(truecols != '') {
    true_cols = as.numeric(unlist(strsplit(truecols, ',')))
    pval_tp = sum(all_pvals[,true_cols] <= thresh)
    pval_fp = sum(all_pvals[,-true_cols] <= thresh)
    pval_fn = sum(all_pvals[,true_cols] > thresh)
    pval_prec = pval_tp / (pval_tp + pval_fp)
    pval_rec = pval_tp / (pval_tp + pval_fn)
    pval_f1 = (2 * pval_prec * pval_rec) / (pval_prec + pval_rec)

    qval_tp = sum(all_qvals[,true_cols] <= qval_thresh)
    qval_fp = sum(all_qvals[,-true_cols] <= qval_thresh)
    qval_fn = sum(all_qvals[,true_cols] > qval_thresh)
    qval_prec = qval_tp / (qval_tp + qval_fp)
    qval_rec = qval_tp / (qval_tp + qval_fn)
    qval_f1 = (2 * qval_prec * qval_rec) / (qval_prec + qval_rec)

    write('Precision/Recall/F1 based on p-values:', file=outfile, append=TRUE)
    write(c(pval_prec, pval_rec, pval_f1), file=outfile, append=TRUE)
    write('Precision/Recall/F1 based on Q-values:', file=outfile, append=TRUE)
    write(c(qval_prec, qval_rec, qval_f1), file=outfile, append=TRUE)
  }
}


eval_sim_pip = function(res, outfile, truecols, verbose=FALSE) {
  num_betas = ncol(res) / 3
  pip_res = res[, (num_betas*2+1):ncol(res)]

  # compute average betas, average std errors, and positive rates
  avg_betas = integer(num_betas)
  sd_betas = integer(num_betas)
  avg_se = integer(num_betas)

  # compute average betas, average std errors, and positive rates
  avg_betas = integer(num_betas)
  sd_betas = integer(num_betas)
  avg_se = integer(num_betas)
  for (i in 1:num_betas) {
    avg_betas[i] = mean(res[,i])
    sd_betas[i] = sd(res[,i])
    avg_se[i] = mean(res[,(i+num_betas)])
  }

  # print results
  write('Average Estimates for Betas: ', file=outfile)
  write(avg_betas, file=outfile, append=TRUE)
  write('Standard Deviations of Betas: ', file=outfile, append=TRUE)
  write(sd_betas, file=outfile, append=TRUE)
  write('Average StdErrors for Betas: ', file=outfile, append=TRUE)
  write(avg_se, file=outfile, append=TRUE)

  # check calibration: true pos rate for each PIP decile
  write('Calibration results: ', file=outfile, append=TRUE)
  truevars = as.numeric(unlist(strsplit(truecols, ",")))
  true_pips = pip_res[, truevars]
  false_pips = pip_res[, -truevars]
  for(i in 0:9) {
    bin_min = i / 10
    bin_max = (i + 1) / 10
    tp = sum(true_pips >= bin_min & true_pips <= bin_max)
    fp = sum(false_pips >= bin_min & false_pips <= bin_max)
    if(tp + fp == 0) {
      write(paste0('No pips in bin from ', bin_min, ' to ', bin_max),
            file=outfile, append=TRUE)
    } else {
      pct = tp / (tp + fp)
      n_samp = tp + fp
      write(paste0('Pct. correct for pips in bin from ', bin_min, ' to ',
                   bin_max,': ', pct, '; n=', n_samp), file=outfile, append=TRUE)
    }
  }

  # compute average pips and percent of time pips meet certain threshold
  num_pips = ncol(pip_res)
  thresh1 <- 0.1
  thresh2 <- 0.5
  thresh3 <- 0.95
  avg_pips = integer(num_pips)
  posrate1 = integer(num_pips)
  posrate2 = integer(num_pips)
  posrate3 = integer(num_pips)
  for (i in 1:num_pips) {
    avg_pips[i] = mean(pip_res[,i])
    posrate1[i] = sum(pip_res[,i] >= thresh1) / length(pip_res[,i])
    posrate2[i] = sum(pip_res[,i] >= thresh2) / length(pip_res[,i])
    posrate3[i] = sum(pip_res[,i] >= thresh3) / length(pip_res[,i])
  }

  # print results
  write('Average PIPs: ', file=outfile, append=TRUE)
  write(avg_pips, file=outfile, append=TRUE)
  write(paste0('Percent of time PIP >= ', thresh1, ': '), file=outfile, append=TRUE)
  write(posrate1, file=outfile, append=TRUE)
  write(paste0('Percent of time PIP >= ', thresh2, ': '), file=outfile, append=TRUE)
  write(posrate2, file=outfile, append=TRUE)
  write(paste0('Percent of time PIP >= ', thresh3, ': '), file=outfile, append=TRUE)
  write(posrate3, file=outfile, append=TRUE)
}


eval_sim_pow_fpr_thresh = function(res, outfile, truecols, is_pip, verbose=FALSE) {
  thresh = 0.05
  truevars = as.numeric(unlist(strsplit(truecols, ",")))
  num_sims = dim(res)[1]
  num_cols = dim(res)[2]
  num_betas = num_cols / 3
  num_falsevars = num_betas - length(truevars)
  tot_false = num_sims * num_falsevars
  falselim = tot_false * thresh

  sorted_pvals = c()
  for(row in 1:num_sims) {
    for(col in 1:num_betas) {
      entry = res[row,(2*num_betas+col)]
      sorted_pvals = rbind(sorted_pvals, c(entry, col))
    }
  }
  sorted_pvals = sorted_pvals[order(sorted_pvals[,1],decreasing=is_pip),]

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
  truepos = truepos / num_sims

  write(paste0('Power at FPR threshold of ', thresh, ': '), file=outfile, append=TRUE)
  write(truepos[truevars], file=outfile, append=TRUE)
}


# given a list of files, return the list of methods that were run
get_method_list = function(file_list) {
  method_list = c()
  for(fname in file_list) {
    fname = unlist(strsplit(fname, '.txt'))
    if(grepl('all_res', fname)) {
      method_name = paste(unlist(strsplit(fname, "_"))[3:4], collapse='_')
    } else {
      method_name = paste(unlist(strsplit(fname, "_"))[4:5], collapse='_')
    }
    method_list = c(method_list, method_name)
  }
  method_list = gsub('_NA', '', unique(method_list))
  return(method_list)
}


# read in results from many single-simulation files into a table
read_results = function(fnames) {
  all_res = c()
  for(fname in fnames) {
    res = read.table(fname)
    all_res = rbind(all_res, res)
  }
  return(na.omit(all_res))
}


# print summary of famr priors and susie results
eval_famr_prior_calibration = function(famr_res_fnames, truecols, n_expo, outfile) {
  all_priors = list('pi' = list(), 'sigma2' = list())
  all_susie_res = list('avg_set_size' = 0, 'total_set_size' = 0,
                       'set_precision' = 0, 'set_recall' = 0, 'coverage' = 0)
  # gather variable names and initialize list
  res1 = readRDS(famr_res_fnames[1])
  varnames = names(res1$priors$pi)
  varnames = varnames[varnames != 'n_vars']  # not an actual prior
  for(vn in varnames) {
    all_priors$pi[[vn]] = c()
    all_priors$sigma2[[vn]] = c()
  }

  # gather prior values and susie results from all famr rds files
  num_nonzero_cs = 0
  truevars = as.numeric(unlist(strsplit(truecols, ",")))
  for(fname in famr_res_fnames) {
    res = readRDS(fname)
    for(vn in varnames) {
      all_priors$pi[[vn]] = c(all_priors$pi[[vn]], res$priors$pi[[vn]])
      all_priors$sigma2[[vn]] = c(all_priors$sigma2[[vn]], res$priors$sigma2[[vn]])
    }

    # factors are true effect variables as well
    if(res$priors$n_vars$exposures > n_expo) {
      this_truevars = c(truevars, (n_expo+1):res$priors$n_vars$exposures)
      this_n_expo = res$priors$n_vars$exposures
    } else if('factors' %in% names(res$priors$n_vars)) {
      this_truevars = c(truevars, (n_expo+1):(n_expo+res$priors$n_vars$factors))
      this_n_expo = n_expo + res$priors$n_vars$factors
    } else {
      this_truevars = truevars
      this_n_expo = n_expo
    }
    susie_res = assess_susie_sets(res$susieres, this_truevars, 1:this_n_expo)
    if(susie_res$avg_set_size > 0) {  # if any sets were returned
      num_nonzero_cs = num_nonzero_cs + 1
      all_susie_res$avg_set_size = all_susie_res$avg_set_size + susie_res$avg_set_size
      all_susie_res$total_set_size = all_susie_res$total_set_size + susie_res$total_set_size
      all_susie_res$set_precision = all_susie_res$set_precision + susie_res$set_precision
      all_susie_res$coverage = all_susie_res$coverage + susie_res$coverage
    }
    all_susie_res$set_recall = all_susie_res$set_recall + susie_res$set_recall
  }

  # summarize and write prior results
  write(paste0('Summary of famr priors: '), file=outfile, append=TRUE)
  for(vn in varnames) {
    vn_pi_mean = mean(all_priors$pi[[vn]])
    vn_pi_sd = sd(all_priors$pi[[vn]])
    vn_s2_mean = mean(all_priors$sigma2[[vn]])
    vn_s2_sd = sd(all_priors$sigma2[[vn]])
    write(paste0('Mean prior PIP for ', vn, ': ', vn_pi_mean), file=outfile, append=TRUE)
    write(paste0('SD of prior PIP for ', vn, ': ', vn_pi_sd), file=outfile, append=TRUE)
    write(paste0('Mean prior sigma^2 for ', vn, ': ', vn_s2_mean), file=outfile, append=TRUE)
    write(paste0('SD of prior sigma^2 for ', vn, ': ', vn_s2_sd), file=outfile, append=TRUE)
  }

  # summarize and write susie results
  ovr_avg_set_size = all_susie_res$avg_set_size / num_nonzero_cs
  ovr_tot_set_size = all_susie_res$total_set_size / num_nonzero_cs
  ovr_avg_precision = all_susie_res$set_precision / num_nonzero_cs
  ovr_avg_recall = all_susie_res$set_recall / length(famr_res_fnames)
  ovr_avg_coverage = all_susie_res$coverage / num_nonzero_cs
  write(paste0('Summary of susie sets: '), file=outfile, append=TRUE)
  write(paste0('Avg. set size: ', ovr_avg_set_size), file=outfile, append=TRUE)
  write(paste0('Avg. exposures in sets: ', ovr_tot_set_size), file=outfile, append=TRUE)
  write(paste0('Avg. precision: ', ovr_avg_precision), file=outfile, append=TRUE)
  write(paste0('Avg. recall: ', ovr_avg_recall), file=outfile, append=TRUE)
  write(paste0('Avg. coverage: ', ovr_avg_coverage), file=outfile, append=TRUE)
}


# assess susie sets for coverage, recall, average precision, etc
assess_susie_sets = function(res, truecols, include) {
  # check if no credible sets
  all_cs = as.numeric(unique(unlist(res$sets$cs)))
  if(length(all_cs) == 0) {
    return(list('avg_set_size' = 0, 'total_set_size' = 0,
                'set_precision' = 0, 'set_recall' = 0, 'coverage' = 0))
  }

  # filter out variables that are not exposures
  res$sets$cs = lapply(res$sets$cs, function(x) x[which(x %in% include)])
  res$sets$cs = res$sets$cs[sapply(res$sets$cs, function(x) length(x) > 0)]

  # compute set lengths, true/false positives, sets with true positives, etc
  all_cs = as.numeric(unique(unlist(res$sets$cs)))
  if(length(all_cs) == 0) {  # check again after filtering non-exposures
    return(list('avg_set_size' = 0, 'total_set_size' = 0,
                'set_precision' = 0, 'set_recall' = 0, 'coverage' = 0))
  }
  tot_cs_len = length(as.numeric(unlist(res$sets$cs)))
  num_cs = length(res$sets$cs)
  true_pos = length(which(truecols %in% all_cs))
  false_pos = length(all_cs) - true_pos
  false_neg = length(truecols) - true_pos
  sets_with_tp = sum(unlist(lapply(res$sets$cs, function(x) length(which(truecols %in% x)) > 0)))

  # summarize into avg set size, precision/recall, and coverage
  avg_set_size = tot_cs_len / num_cs
  set_precision = true_pos / (true_pos + false_pos)
  set_recall = true_pos / (true_pos + false_neg)
  coverage = sets_with_tp / num_cs
  return(list('avg_set_size' = avg_set_size, 'total_set_size' = tot_cs_len,
              'set_precision' = set_precision, 'set_recall' = set_recall, 'coverage' = coverage))
}


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(opt$out == 'DEFAULT') {
  opt$out = opt$indir
}
num_expo = 0
file_list = list.files(path = opt$indir, pattern = '.*res.*.txt')
method_list = get_method_list(file_list)
for(method in method_list) {
  if(grepl('TEMP', method))  next
  print(method)
  m2 = gsub('\\.', '\\\\.', method)  # escape '.' for string matching
  method_files = file_list[grepl(paste0(m2, '.txt'), file_list)]
  # check if all_res file already exists; if so, don't rewrite it
  allres_idx = which(grepl('all_res', method_files))
  if(length(allres_idx) > 0) {  # just read in all_res file
    all_resname = paste0(opt$out, method_files[allres_idx])
    if (file.size(all_resname) == 0)  next
    all_res = na.omit(read.table(all_resname))
  } else {  # all_res file = concatenated results of individual job files
    all_res = read_results(paste0(opt$indir, method_files))
    all_resname = paste0(opt$out, 'all_res_', method, '.txt')
    write.table(all_res, all_resname, row.names=F, col.names=F, quote=F)
  }
  if(length(all_res) == 0)  next
  num_expo = ncol(all_res) / 3
  performance_resname = paste0(opt$out, 'performance_', method, '.txt')

  # now run the evalutions, e.g. compute and print the metrics
  if(grepl('susie', method) || grepl('mrash', method) || grepl('varbvs', method)
     || grepl('brms', method) || grepl('vebboost', method) || grepl('bma', method)
     || grepl('ctwas', method) || grepl('famr', method)) {
    eval_sim_pip(all_res, performance_resname, opt$truecols, verbose=opt$verbose)
    eval_sim_pow_fpr_thresh(all_res, performance_resname, opt$truecols,
                            TRUE, verbose=opt$verbose)
  } else{
    eval_sim_nonpip(all_res, performance_resname,
                    truecols=opt$truecols, verbose=opt$verbose)
    eval_sim_pow_fpr_thresh(all_res, performance_resname, opt$truecols,
                            FALSE, verbose=opt$verbose)
  }
  # if(grepl('famr', method) || grepl('ctwas', method)) {
  #   # additionally check if famr priors are calibrated
  #   famr_res_files = list.files(path = opt$indir, pattern = paste0('job.*', method, '.rds'), full.names = T)
  #   if(length(famr_res_files) > 0) {
  #     eval_famr_prior_calibration(famr_res_files, opt$truecols, num_expo, performance_resname)
  #   }
  # }
}
