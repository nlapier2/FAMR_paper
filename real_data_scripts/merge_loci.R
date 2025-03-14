# Merge per-trait loci across all traits, then split into individual locus files
# In other words, convert the single-trait, all-regions files from extract_loci
#   to all-traits, single-region files

### IMPORTS

library(optparse)
library(dplyr)
library(data.table)

### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--loci_dir"), type="character", default='NONE',
              help="Directory with locus rds files from extract_loci.R", metavar="character"),
  make_option(c("--merged_regions"), type="character", default='NONE',
              help="File with merged region files", metavar="character"),
  make_option(c("--region"), type="numeric", default=0,
              help="Which region to process (default: 0 = all)", metavar="numeric"),
  make_option(c("--unmerged"), type="logical", action="store_true", default=FALSE,
              help="Use if running only on an unmerged region file"),
  make_option(c("--out"), type="character", default='res_',
              help="Prefix of output and temporary files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# take summary statistics from one file and inner join it with those from other files
append_sumstats = function(sumstats, dat, trait) {
  z.trait = paste0('z|', trait)
  b.trait = paste0('b|', trait)
  s.trait = paste0('s|', trait)
  if("T_STAT" %in% names(dat)) {
    dt = dat  %>% rename(!!z.trait := T_STAT, !!b.trait := BETA, !!s.trait := SE)
  } else {
    dt = dat  %>% rename(!!z.trait := Z, !!b.trait := BETA, !!s.trait := SE)
  }
  dt$CHR = as.numeric(dt$CHR)
  if(is.null(nrow(sumstats))){
    sumstats = dt
  } else{
    sumstats = inner_join(sumstats, dt, by=c('ID', 'CHR', 'POS', 'REF', 'ALT'))
  }
  sumstats$CHR = as.numeric(sumstats$CHR)
  return(sumstats)
}


# write the merged summary statistics file
write_sumstats = function(sumstats, outname) {
  pos = sumstats %>% select(CHR, POS, ID, REF, ALT)
  zscores  = sumstats %>% select(starts_with('z|')) %>% rename_with(~ gsub("^z\\|", "", .x))
  betas  = sumstats %>% select(starts_with('b|')) %>% rename_with(~ gsub("^b\\|", "", .x))
  stderrs  = sumstats %>% select(starts_with('s|')) %>% rename_with(~ gsub("^s\\|", "", .x))
  saveRDS(list(pos=pos, Z = zscores, betas = betas, stderrs = stderrs), outname)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # get file names and region, set output name
  all.fnames = list.files(path = opt$loci_dir, pattern = '*.rds')

  # read merged regions file
  merged_regions = read.table(opt$merged_regions)
  if(opt$unmerged) {  # slightly different columns for unmerged regions
    names(merged_regions) = c('CHR', 'MIDDLE', 'MAXZ', 'START', 'END')
  } else {
    names(merged_regions) = c('CHR', 'START', 'END', 'LENGTH', 'MAXZ')
  }
  if(opt$region > 0) {  # select region specified by user, if appropriate
    merged_regions = merged_regions[opt$region, , drop=F]
    chr = merged_regions[1,1]
    if(opt$unmerged) {
      start = merged_regions[1,4]
      end = merged_regions[1,5]
    } else {
      start = merged_regions[1,2]
      end = merged_regions[1,3]
    }
    outname = paste0(opt$out, 'sumstats_chr_', chr, '_', start, '-', end, '.rds')
  } else {
    outname = paste0(opt$out, 'sumstats_all_regions.rds')
  }
  if(file.exists(outname)) {
    stop(paste0(outname, ': file already exists; will not overwrite.'))
  }

  sumstats = list()
  for(fname in all.fnames) {
    full.fname = paste0(opt$loci_dir, fname)
    trait = strsplit(fname, '\\.')[[1]][1]
    dat = readRDS(full.fname)
    # remove NaN values, non-standard REF/ALT, and duplicate SNPs
    dat = na.omit(dat)
    dat = dat %>% filter(REF %in% c('A', 'C', 'G', 'T') & ALT %in% c('A', 'C', 'G', 'T'))
    dat = dat %>% filter (! duplicated(POS) & !duplicated(POS, fromLast = TRUE))
    # Create a combined condition for all regions
    combined_condition <- merged_regions %>%
      mutate(CHR = as.numeric(CHR), START = as.numeric(START), END = as.numeric(END)) %>%
      rowwise() %>%
      mutate(condition = list(dat$CHR == CHR & dat$POS >= START & dat$POS <= END)) %>%
      ungroup() %>%
      pull(condition) %>%
      Reduce(`|`, .)
    dat = dat %>% filter(combined_condition)
    # append these sumstats onto the sumstats over all regions
    sumstats = append_sumstats(sumstats, dat, trait)
    sumstats = sumstats %>% arrange(POS)
  }
  write_sumstats(sumstats, outname)
}
