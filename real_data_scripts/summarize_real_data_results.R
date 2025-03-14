# Summarize real data FAMR, MVMR, and UVMR results into a handy table

### IMPORTS

library("optparse")
library(data.table)
library(dplyr)
library(qvalue)


### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--in_dir"), type="character", default='./',
              help="Directory of FAMR and MVMR results files", metavar="character"),
  make_option(c("--uvmr_dir"), type="character", default='./',
              help="Directory of univariate MR results files", metavar="character"),
  make_option(c("--trait_names"), type="character", default='NONE',
              help="File mapping trait names to trait IDs", metavar="character"),
  make_option(c("--extension"), type="character", default='_loci',
              help="Extension to trait names that should be removed", metavar="character"),
  make_option(c("--out"), type="character", default='results.rds',
              help="Output results file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# helper function to extract outcome (and possibly exposure) trait name(s) from a file name
extract_trait_name = function(fname, trait_names, delim='_', get_expo=F) {
  splits = strsplit(strsplit(basename(fname), '\\.')[[1]][1], delim)[[1]]
  outcome_trait_id = splits[[length(splits)]]
  outcome_name = trait_names[trait_names[['Trait ID']] == outcome_trait_id, ][['Trait Name']]
  if(get_expo) {
    exposure_trait_id = strsplit(splits[[1]][1], 'res_')[[1]][2]
    exposure_name = trait_names[trait_names[['Trait ID']] == exposure_trait_id, ][['Trait Name']]
    return(list('exp' = exposure_name, 'outc' = outcome_name))
  } else {
    return(outcome_name)
  }
}


# gather univariate MR results across each exposure for an outcome, compute q-values,
#   and summarize those q-values in a dataframe
process_uvmr = function(uvmr_fnames, trait_names) {
  # collect the univariate MR p-values for each exposure, grouped by outcome
  per_outc_res = list()
  for(fname in uvmr_fnames) {
    traits = extract_trait_name(fname, trait_names, delim='---', get_expo=T)
    pval = as.numeric(fread(fname)$V3)
    if(traits$outc %in% names(per_outc_res)) {
      rows = nrow(per_outc_res[[traits$outc]])
      per_outc_res[[traits$outc]][rows+1, ] = c(traits$exp, pval)
    } else {
      per_outc_res[[traits$outc]] = data.frame('Exposure' = traits$exp, 'pvalue' = pval)
    }
  }

  # for each outcome, compute the q-values and add to the final results dataframe
  all_df = data.frame('Exposure' = c(), 'pvalue' = c(), 'ivw.univar' = c(), 'Outcome' = c())
  for(outcome in names(per_outc_res)) {
    qvals = qvalue(as.numeric(per_outc_res[[outcome]]$pvalue), lambda=0)$qvalues
    per_outc_res[[outcome]]$ivw.univar = 1 - qvals
    per_outc_res[[outcome]]$Outcome = outcome
    all_df = rbind(all_df, per_outc_res[[outcome]])
  }
  all_df[['pvalue']] = NULL
  return(all_df)
}


# extract MVMR results
process_mvmr = function(mvmr_fnames, trait_names, extension) {
  all_df = data.frame()
  for(fname in mvmr_fnames) {
    # read in data
    res = readRDS(fname)

    # extract outcome and exposure names and add results to the dataframe
    expo_ids = sapply(res[[names(res)[1]]]$name, function(x) strsplit(x, extension)[[1]][1])
    df = data.frame(Exposure = sapply(expo_ids, function(x)
      trait_names[trait_names[['Trait ID']] == x, ][['Trait Name']]))
    df$Outcome = extract_trait_name(fname, trait_names, delim='_', get_expo=F)

    # add results for each MVMR method to the dataframe, then append to overall df
    for(method in names(res)) {
      df[[method]] = 1 - res[[method]]$qvalue
    }
    if(length(all_df) == 0) {
      all_df = df
    } else {
      all_df = rbind(all_df, df)
    }
  }
  return(all_df)
}


# extract FAMR results
process_famr = function(famr_fnames, trait_names, extension) {
  all_df = data.frame()
  for(fname in famr_fnames) {
    # read in data
    res = readRDS(fname)
    df = data.frame(res$pips$exposures)

    # extract outcome and exposure names and add results to the dataframe
    outcome_name = extract_trait_name(fname, trait_names, delim='_', get_expo=F)
    df$Outcome = outcome_name
    expo_ids = sapply(rownames(df), function(x) strsplit(x, extension)[[1]][1])
    df$Exposure = sapply(expo_ids, function(x)
      as.character(trait_names[trait_names[['Trait ID']] == x, ][['Trait Name']]))
    rownames(df) = NULL

    # rename the column to reflect the present method
    fa_method = strsplit(basename(fname), '_')[[1]][3]
    setnames(df, 'res.pips.exposures', paste0('famr_', fa_method))
    df$Exposure = as.character(df$Exposure)

    # add to overall results dataframe
    if(length(all_df) == 0) {
      all_df = df
    } else {
      all_df = rbind(all_df, df)
    }
  }
  return(all_df)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # read trait names and list files to process
  all_results = data.frame()
  trait_names = fread(opt$trait_names, sep='\t', header=T)
  famr_none_fnames = Sys.glob(paste0(opt$in_dir, '/results_famr_none*.rds'))
  famr_gfa_fnames = Sys.glob(paste0(opt$in_dir, '/results_famr_gfa*.rds'))
  mvmr_fnames = Sys.glob(paste0(opt$in_dir, '/results_mvmr*.rds'))
  uvmr_fnames = Sys.glob(paste0(opt$uvmr_dir, '*.txt'))

  # read in results and join them together in a table, then write that table
  uvmr_results = process_uvmr(uvmr_fnames, trait_names)
  mvmr_results = process_mvmr(mvmr_fnames, trait_names, opt$extension)
  famr_gfa_results = process_famr(famr_gfa_fnames, trait_names, opt$extension)
  all_results = full_join(famr_gfa_results, mvmr_results, by=c('Exposure', 'Outcome'))
  if(length(famr_none_fnames) > 0) {
    famr_none_results = process_famr(famr_none_fnames, trait_names, opt$extension)
    all_results = full_join(all_results, famr_none_results, by=c('Exposure', 'Outcome'))
  }
  all_results = full_join(all_results, uvmr_results, by=c('Exposure', 'Outcome'))
  all_results = all_results %>% relocate(c('Exposure', 'Outcome'))
  all_results = na.omit(all_results)  # removes factors and other unknown exposures
  write.table(all_results, file = opt$out, quote=F, sep='\t', row.names=F)
}
