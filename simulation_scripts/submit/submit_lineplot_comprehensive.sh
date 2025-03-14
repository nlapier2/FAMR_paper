#!/bin/bash
#SBATCH --job-name=job-ukbb-sim-mvmr
#SBATCH --output=job-ukbb-sim-mvmr.out
#SBATCH --error=job-ukbb-sim-mvmr.err
#SBATCH --time=23:00:00
#SBATCH --partition=caslake  # general-purpose compute partition
#SBATCH --nodes=1  # 1 node requested
#SBATCH --ntasks-per-node=1  # 1 core requested per node
#SBATCH --mem-per-cpu=16000
#SBATCH --array=67-100  # array of jobs. access job id with $SLURM_ARRAY_TASK_ID.
#SBATCH --account=pi-mstephens

. ~/.bash_profile

# default parameter values
tld=results_lineplot_comprehensive   # top-level directory for results; do not add /
[ ! -d ${tld} ] && mkdir ${tld}
default_other="-i /project2/mstephens/yuxin/ukb-bloodcells/genotypes/ --ld_dir /project2/mstephens/yuxin/ukb-bloodcells/LD/ --num_samples 10000 --num_loci 100 --snps_per_locus 100 --num_confounders 3 --max_causal 5 --phi_gy 0.05 --phi_gz 0.2 --phi_zx 0.5 --susieL 30 --methods=famr_gfa,ivw,ivw_gfa,median,median_gfa --fa_prune_thresh 0.1 --famr_n_iter 30 --final_prune_thresh 0.1 --two_sample"
default_gx=0.2
default_psix=0.1 
default_psiy=0.05
default_psiw=0.0
default_gz=0.2  
default_gw=0.4
default_zx=`echo "scale=8 ; $default_psix / $default_gz" | bc`
default_zy=`echo "scale=8 ; $default_psiy / $default_gz" | bc`
default_wx=`echo "scale=8 ; $default_psiw / $default_gw" | bc`
default_gy=0.0
default_mu=0.0
default_n_expo=30
default_prune=0.5
default_beta="0.05,0.1,0.2,0.3"
default_phi_gx=0.05

# parameter values in array other than defaults
psiy_vals=(0.0 0.025 0.05 0.075 0.1)
beta_vals=("0.0,0.0,0.0,0.0" "0.1,0.1,0.1,0.1" "0.2,0.2,0.2,0.2" "0.3,0.3,0.3,0.3")
gy_vals=(0.0 0.05 0.1)
n_expo_vals=(10 60)
prune_vals=(0.1 0.3 0.5 0.75)
phi_gx_vals=(0.01 0.025 0.1 0.2)

args=${default_other}" --gamma_gx $default_gx --gamma_zx $default_zx --mu $default_mu --gamma_gw $default_gw --gamma_wx $default_wx --gamma_gz $default_gz"

sims_per_job=1
job=$((${SLURM_ARRAY_TASK_ID}-1))
base=$((${job}*${sims_per_job}))

### Vary psi_y parameter

for this_psiy in "${psiy_vals[@]}"
do
        this_zy=`echo "scale=8 ; $this_psiy / $default_gz" | bc`
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy $this_zy --gamma_gy $default_gy --beta=$default_beta --num_exposures $default_n_expo --famr_prune $default_prune --criterion max --phi_gx $default_phi_gx"
        outdir=${tld}/array_zy_psiy_${this_psiy}_beta_${default_beta}_gy_${default_gy}_nexpo_${default_n_expo}_prune_max/ 
	[ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done

for this_beta in "${beta_vals[@]}"
do
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy $default_zy --gamma_gy $default_gy --beta=$this_beta --num_exposures $default_n_expo --famr_prune $default_prune --criterion max --phi_gx $default_phi_gx"
        outdir=${tld}/array_beta_psiy_${default_psiy}_beta_${this_beta}_gy_${default_gy}_nexpo_${default_n_expo}_prune_max/
        [ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done

for this_gy in "${gy_vals[@]}"
do
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy 0.0 --gamma_gy $this_gy --beta=$default_beta --num_exposures $default_n_expo --famr_prune $default_prune --criterion max --phi_gx $default_phi_gx"
        outdir=${tld}/array_gy_psiy_0.0_beta_${default_beta}_gy_${this_gy}_nexpo_${default_n_expo}_prune_max/
        [ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done

for this_n_expo in "${n_expo_vals[@]}"
do
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy $default_zy --gamma_gy $default_gy --beta=$default_beta --num_exposures $this_n_expo --famr_prune $default_prune --criterion max --phi_gx $default_phi_gx"
        outdir=${tld}/array_nexpo_psiy_${default_psiy}_beta_${default_beta}_gy_${default_gy}_nexpo_${this_n_expo}_prune_max/
        [ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done

for this_prune in "${prune_vals[@]}"
do
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy $default_zy --gamma_gy $default_gy --beta=$default_beta --num_exposures $default_n_expo --famr_prune $this_prune --criterion clumped --phi_gx $default_phi_gx"
        outdir=${tld}/array_prune_psiy_${default_psiy}_beta_${default_beta}_gy_${default_gy}_nexpo_${default_n_expo}_prune_${this_prune}/
        [ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done

for this_phi_gx in "${phi_gx_vals[@]}"
do
        cmd="Rscript ukbb_run_sim_and_methods.R "${args}" --gamma_zy $default_zy --gamma_gy $default_gy --beta=$default_beta --num_exposures $default_n_expo --famr_prune $default_prune --criterion max --phi_gx $this_phi_gx"
        outdir=${tld}/array_phigx_psiy_${default_psiy}_beta_${default_beta}_gy_${default_gy}_nexpo_${default_n_expo}_prune_max_phigx_${this_phi_gx}/
	[ ! -d ${outdir} ] && mkdir ${outdir}
        echo ${cmd} > ${outdir}info.txt
        ${cmd} --num_sims ${sims_per_job} --out ${outdir}job_${job}_res_
done
#
