# Given a directory containing VCF-format summary statistics files for multiple
#   phenotypes, extract regions (loci) surrounding genome-wide significant SNPs

### IMPORTS

library(optparse)
library(dplyr)
library(data.table)
source("merge_regions.R")
source("utils.R")

### COMMAND LINE OPTION PARSING

option_list = list(
  make_option(c("--vcf"), type="character", default='NONE',
              help="Input VCF summary statistics file", metavar="character"),
  make_option(c("--out"), type="character", default='res_',
              help="Output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# find regions with a p-value stronger than 5e-8 and center a window around each
#   such peak to create regions (loci)
extract_trait_regions = function(gwas) {
  res = c()
  while(max(gwas$logp, na.rm=T) > -log10(5e-8)){
    signal = gwas %>% filter(logp == max(logp, na.rm=T)) %>% top_n(1, POS) %>% select(CHR, POS, logp) %>%
      mutate(START = max(POS - 500000, 0), END = POS + 500000)
    res = rbind(res, signal)
    gwas[which(gwas$CHR == signal$CHR & gwas$POS <= signal$END & gwas$POS >= signal$START),] = NA
  }
  return(res)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  dat = read_vcf(opt$vcf)
  res = extract_trait_regions(dat)
  if(length(res) > 0)  res = merge_regions(list('regions' = res))
  write.table(res, file = opt$out, quote = F, row.names = F, col.names = F)
}
