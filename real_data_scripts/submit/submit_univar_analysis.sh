#!/bin/bash
#SBATCH --job-name=univar_analysis
#SBATCH --output=job-univar_analysis.out
#SBATCH --error=job-univar_analysis.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=32000
#SBATCH --account=pi-mstephens
#SBATCH --array=706-1410 #1410 #705  #165  # 320

. ~/.bash_profile

outdir=results_mvmr_interval_nightingale_new_outc/univar  #results_mvmr_interval_nightingale/univar  #results_metc/univar

# select the exposure and outcome this job in the array will run on
job=${SLURM_ARRAY_TASK_ID}
nexpo=141  # 33
e_idx=$(( ($job - 1) % $nexpo  ))  # which exposure to run on in this job
o_idx=$(( ($job - 1) / $nexpo  ))  # which outcome  to run on in this job

# get exposure file and corresponding regions to draw instruments from
e_files=(../interval_data/nightingale/vcfs/*.vcf.gz)  # (raw_vcfs_combined/*.vcf.gz)  #(raw_vcfs_metc/*.vcf.gz)
this_e_infile=${e_files[$e_idx]}
e_trait_id=$(basename $this_e_infile | cut -d '.' -f 1)
region_fname=../interval_data/nightingale/loci_regions/${e_trait_id}.tar.gz.regions.txt  #regions_combined_newexpo/${e_trait_id}_regions.txt  # regions to draw instruments from for this exposure

# get outcome file
o_files=(raw_vcfs_interval_outcomes/*.vcf.gz)  #(raw_vcfs_morrison/outcomes/*.vcf.gz)
this_o_infile=${o_files[$o_idx]}
o_trait_id=$(basename $this_o_infile | cut -d '.' -f 1)

# set output file name
this_outfile=${outdir}/res_${e_trait_id}---${o_trait_id}.txt

Rscript univar_real_data_analysis_vcf.R --exposure_vcf $this_e_infile --outcome_vcf $this_o_infile --regions $region_fname --out $this_outfile --old_col_order

