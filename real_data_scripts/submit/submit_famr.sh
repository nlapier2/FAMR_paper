#!/bin/bash
#SBATCH --job-name=famr_real_data_analysis
#SBATCH --output=job-famr_real_data_analysis.out
#SBATCH --error=job-famr_real_data_analysis.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=32000
#SBATCH --account=pi-mstephens
#SBATCH --array=1-1 #1-10

. ~/.bash_profile

# general options
outdir=results_famr_susie_interval_nightingale_new_outc/ #_no_lipo_subfrac/  #results_metc/mvmr/
outcome_dir=raw_vcfs_interval_outcomes/  #raw_vcfs_morrison/outcomes/
in_dir=subset_sumstats_interval_nightingale/ #_no_lipo_subfrac/
ld_dir=ukbb_ld_interval_nightingale/ #_no_lipo_subfrac/
fa_method="gfa" #"none"
N=37359 #20000  #115078
susieL=30
prune_thresh=0.5
fa_prune_thresh=0.1
final_prune_thresh=0.1
n_iter=30

# get outcome trait corresponding to this job, set output file name, create output directory
job=${SLURM_ARRAY_TASK_ID}
outc_files=(${outcome_dir}/*.vcf.gz)
outcome_vcf=${outc_files[$((job-1))]}
trait_id=$(basename $outcome_vcf | cut -d '.' -f 1)
outname=${outdir}/results_famr_${fa_method}_noannih_${trait_id}.rds
[ ! -d ${outdir} ] && mkdir ${outdir}
log_file=${outdir}/logfile_famr_${fa_method}_noannih_${trait_id}


# run and log lines -- generally do not change
cmd="Rscript famr.R --in_dir ${in_dir} --ld_dir ${ld_dir} --y_gwas_file ${outcome_vcf} --fa_method ${fa_method} -N ${N} --susieL ${susieL} --prune_thresh ${prune_thresh} --fa_prune_thresh ${fa_prune_thresh} --final_prune_thresh ${final_prune_thresh} --n_iter ${n_iter} --out ${outname} --idcol 1 --betacol 2 --secol 3"  # --annihilate_factors"

echo 'command: ' > ${log_file}
echo ${cmd} >> ${log_file}
echo 'start: ' >> ${log_file}
date >> ${log_file}
time ${cmd} >> ${log_file} 2>&1  # here is where we actually execute the analysis command
echo 'end: ' >> ${log_file}
date >> ${log_file}
#
