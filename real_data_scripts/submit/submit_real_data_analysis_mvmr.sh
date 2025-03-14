#!/bin/bash
#SBATCH --job-name=real_data_analysis
#SBATCH --output=job-real_data_analysis.out
#SBATCH --error=job-real_data_analysis.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=32000
#SBATCH --account=pi-mstephens
#SBATCH --array=1-10

. ~/.bash_profile

# general options
outdir=results_mvmr_interval_nightingale_new_outc_no_lipo_subfrac/  #results_metc/mvmr_noannih_max/
methods=susie.ss,mrash.ss,susie.ss_gfa,mrash.ss_gfa,ivw,ivw_gfa,median,median_gfa  #robust,robust_gfa,brms.hs,brms.hs_gfa
N=37359  #115078  # 20000 
outcome_dir=raw_vcfs_interval_outcomes/  #raw_vcfs_morrison/outcomes
sumstat_dir=subset_sumstats_interval_nightingale_no_lipo_subfrac/  #subset_sumstats_metc/  #merged_loci_combined/
ld_dir=ukbb_ld_interval_nightingale_no_lipo_subfrac/  #ukbb_ld_metc/
sel_method=max #clumped
prune_thresh=0.1


# get outcome trait corresponding to this job, set output file name, create output directory
job=${SLURM_ARRAY_TASK_ID}
outc_files=(${outcome_dir}/*.vcf.gz)
outcome_vcf=${outc_files[$((job-1))]}
trait_id=$(basename $outcome_vcf | cut -d '.' -f 1)
outname=${outdir}/results_mvmr_${sel_method}_${trait_id}.rds
[ ! -d ${outdir} ] && mkdir ${outdir}
log_file=${outdir}/logfile_${sel_method}_${trait_id}


# run and log lines -- generally do not change
cmd="Rscript real_data_analysis_vcf.R --sumstat_dir ${sumstat_dir} --ld_dir ${ld_dir} --outcome_vcf ${outcome_vcf} --sel_method ${sel_method} --N ${N} --methods ${methods} --prune_thresh ${prune_thresh} --out ${outname}"

echo 'command: ' > ${log_file}
echo ${cmd} >> ${log_file}
echo 'start: ' >> ${log_file}
date >> ${log_file}
time ${cmd} >> ${log_file} 2>&1  # here is where we actually execute the analysis command
echo 'end: ' >> ${log_file}
date >> ${log_file}
#

