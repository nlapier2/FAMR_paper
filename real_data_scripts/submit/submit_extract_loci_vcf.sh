#!/bin/bash
#SBATCH --job-name=vcf_extract_loci
#SBATCH --output=job-vcf_extract_loci.out
#SBATCH --error=job-vcf_extract_loci.err
#SBATCH --time=01:59:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=64000
#SBATCH --account=pi-mstephens
#SBATCH --array=1-43 #33

. ~/.bash_profile

region_file=../interval_data/nightingale/merged_regions_no_lipo_subfrac.txt  #merged_regions/merged_regions_combined_newexpo_1Mb.txt
outdir=../interval_data/nightingale/trait_loci_no_lipo_subfrac  #trait_loci_combined_newexpo

# get file corresponding to this job, set output file name
job=${SLURM_ARRAY_TASK_ID}
files=(../interval_data/nightingale/vcfs_no_lipo_subfraction/*.vcf.gz)
#files=(raw_vcfs_combined/*.vcf.gz)  #(raw_vcfs_metc/*.vcf.gz)
this_vcf=${files[$((job-1))]}
trait_id=$(basename $this_vcf | cut -d '.' -f 1)
outname=${outdir}/${trait_id}_loci.rds

Rscript extract_loci_vcf.R --exposure_vcf $this_vcf --merged_regions $region_file --out $outname

