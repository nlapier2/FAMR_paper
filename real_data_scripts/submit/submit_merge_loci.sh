#!/bin/bash
#SBATCH --job-name=merge_loci
#SBATCH --output=job-merge_loci.out
#SBATCH --error=job-merge_loci.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=16000
#SBATCH --account=pi-mstephens
#SBATCH --array=1-166 #728

. ~/.bash_profile

Rscript merge_loci.R --loci_dir ../interval_data/nightingale/trait_loci_no_lipo_subfrac/ --merged_regions ../interval_data/nightingale/merged_regions_no_lipo_subfrac.txt --out ../interval_data/nightingale/merged_loci_interval_no_lipo_subfrac/ --region ${SLURM_ARRAY_TASK_ID}

#Rscript merge_loci.R --loci_dir trait_loci_combined_newexpo/ --merged_regions merged_regions/merged_regions_combined_newexpo_1Mb.txt --out merged_loci_combined_newexpo/ --region ${SLURM_ARRAY_TASK_ID}

