# Given per-phenotype GWAS loci, merge into non-overlapping regions across
#   all phenotypes.

### IMPORTS

library(optparse)
library(dplyr)
library(data.table)

### COMMAND LINE OPTION PARSING

if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  option_list = list(
    make_option(c("--region_dir"), type="character", default='NONE',
                help="Directory with per-phenotype regions", metavar="character"),
    make_option(c("--trait_list"), type="character", default='NONE',
                help="File with list of traits to merge", metavar="character"),
    make_option(c("--old_col_order"), type="logical", action="store_true", default=FALSE,
                help="Use older region file column order (e.g. for INTERVAL)"),
    make_option(c("--out"), type="character", default='res_',
                help="Output file name", metavar="character")
  );

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}


merge_regions = function(init_regions) {
  proc_regions = list()
  for(trait in names(init_regions)){
    region = init_regions[[trait]]
    region = region %>% arrange(desc(logp))
    region_r = c()
    for(i in 1:22){
      region.chr = region %>% filter(CHR == i) %>% arrange(START)
      if(nrow(region.chr) == 0){
        next
      }
      tmp = region.chr %>% group_by(g = cumsum(cummax(lag(END, default = first(END))) < START)) %>%
        summarise(START = first(START), END = max(END), logp = max(logp),.groups = 'drop') %>%
        mutate(length = END - START) %>%
        mutate(CHR = i) %>% select(CHR, START, END, length, logp)
      region_r = rbind(region_r, tmp)
    }
    proc_regions[[trait]] = region_r
  }

  tb = bind_rows(proc_regions, .id = "column_label")
  res.final = c()
  for(i in 1:22){
    tb.chr = tb %>% filter(CHR == i) %>% arrange(START)
    if(nrow(tb.chr) == 0){
      next
    }
    tmp = tb.chr %>% group_by(g = cumsum(cummax(lag(END, default = first(END))) < START)) %>%
      summarise(START = first(START), END = max(END), logp = max(logp), .groups = 'drop') %>%
      mutate(length = END - START) %>%
      mutate(CHR = i) %>% select(CHR, START, END, length, logp)
    res.final = rbind(res.final, tmp)
  }
  return(res.final)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # file names of all traits with available regions in regions_dir
  all.fnames = list.files(path = opt$region_dir, pattern = '*', full.names = T)
  # which traits to use in this analysis
  if(opt$trait_list != 'NONE')  traits_to_use = fread(opt$trait_list, header = T)
  if('Trait ID' %in% names(traits_to_use))  traits_to_use[['gwas_id']] = traits_to_use[['Trait ID']]

  # store per-file regions, skipping those not in traits_to_use
  trait_regions = list()
  for(fname in all.fnames) {
    if(file.size(fname) == 0L)  next
    if(opt$trait_list != 'NONE') {
      trait_id = unlist(strsplit(basename(fname), '\\.'))[1]
      if(!(trait_id %in% traits_to_use$gwas_id))  next
    }

    dat = read.table(fname)
    names(dat) = c('CHR', 'START', 'END', 'LENGTH', 'logp')
    if(opt$old_col_order)  names(dat) = c('CHR', 'POS', 'logp', 'START', 'END')
    trait_regions[[fname]] = dat
  }

  # merge regions, then write the merged regions file
  merged_regions = merge_regions(trait_regions)
  write.table(merged_regions, file = opt$out, quote = F, row.names = F, col.names = F)
}
