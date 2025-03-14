# Utility script for converting INTERVAL-format TSV files to the VCF file format
#   used by IEU, which extract_loci_vcf.R and univar_real_data_analysis.R use.

### IMPORTS

library(optparse)
library(data.table)
library(dplyr)

### COMMAND LINE OPTION PARSING

if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  option_list = list(
    make_option(c("--tsv"), type="character", default='NONE',
                help="Input tsv file", metavar="character"),
    make_option(c("--out"), type="character", default='out.vcf',
                help="Output VCF file name", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  trait_id = unlist(strsplit(basename(opt$tsv), '\\.'))[1]
  temp_fname = paste0(opt$tsv, 'TEMP')
  system(paste0('zcat ', opt$tsv, ' > ', temp_fname))  # extract from .tar.gz
  dat = fread(temp_fname, sep='\t', skip=1)
  system(paste0('rm ', temp_fname))
  dat = dat[dat$V2 != 'CHR', ]  # removes extra header lines
  names(dat) = c('ID', '#CHROM', 'POS', 'GENPOS', 'REF', 'ALT', 'INFO', 'tmp', 'BETA', 'SE', 'P', 'Pother')
  dat$P = as.numeric(dat$P)
  
  # add a few fields needed for VCF
  dat$QUAL = '.'
  dat$FILTER = 'PASS'
  dat$FORMAT = 'ES:SE:LP:AF:ID'
  # now add the formatted field for the trait
  dat[[trait_id]] <- paste(dat$BETA, dat$SE, log10(dat$P), dat$INFO, dat$ID, sep = ":")
  # redorder and select fields
  dat = dat %>% select('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', all_of(trait_id))
  dat$INFO = paste0('AF=', dat$INFO)
  
  # write the data.frame as gzipped vcf file
  gz1 <- gzfile(opt$out, "w")
  write.table(dat, gz1, sep='\t', row.names=F, quote=F)
  close(gz1)
}