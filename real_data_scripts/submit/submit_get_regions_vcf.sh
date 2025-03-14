#!/bin/bash
#SBATCH --job-name=vcf_get_regions
#SBATCH --output=job-vcf_get_regions.out
#SBATCH --error=job-vcf_get_regions.err
#SBATCH --time=08:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=16000
#SBATCH --account=pi-mstephens
#SBATCH --array=1-33

. ~/.bash_profile

# get file corresponding to this job, set output file name
job=${SLURM_ARRAY_TASK_ID}
files=(raw_vcfs_combined/*.vcf.gz) #(raw_vcfs_metc/*.vcf.gz)
this_vcf=${files[$((job-1))]}
trait_id=$(basename $this_vcf | cut -d '.' -f 1)

Rscript get_regions_vcf.R --vcf $this_vcf --out regions_combined_newexpo/${trait_id}_regions.txt

