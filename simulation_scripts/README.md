# Simulation scripts for the FAMR paper

This folder contains scripts used to perform simulations and generate results
for the FAMR paper.


### Simulating phenotypes

Phenotypes were simulated according to the DAG in the FAMR paper 
(see the paper linked on the top-level README).
The phenotype simulation is performed by `sim_mvmr_sparse.R`.


### Scripts to run the methods

FAMR is run by the [FAMR package](https://github.com/nlapier2/FAMR/).
MR-BMA code is in the `mr-bma.R` script, and code for all other methods is
in the `methods.R` script. `run_methods.R` is a wrapper script to run a list of
these methods, adding factors as extra exposures if specified by the user.

`ukbb_run_sim_and_methods.R` is the big wrapper script around all of this
that takes in input data, simulates phenotypes using it, and runs all methods.
The input data (genotypes) are sourced from the UK Biobank, so they cannot
be included here. We randomly downsampled genotypes from independent loci
selected for analysis in the
[mvSuSiE paper](https://doi.org/10.1101/2023.04.14.536893).
Please consult their paper for more details.

`utils.R` is a utility script that contains common functions used by many of the
scripts above.


### Submit and evaluate scripts

The results in the paper are averaged over many simulation replicates.
In practice we did this by running array jobs on our compute cluster, with one
job per simulation, and each simulation writing separate results files.
The scripts to do this, which also give the exact paramter settings used,
are in the `submit/` folder.
The main script for generating most of the simulation results is
`submit_lineplot_comprehensive.sh`; the others were used to generate the
supplementary figure comparing FAMR-Susie versions with different parts
disabled.

After all simulation replicates were run, `eval_all_sim_res.R` was used to
summarize simulation results over all replicates and compute metrics.

For more details on each of these, please peruse the code.
For examples of how to run them, please see the `submit/` folder.


##### Note

We note for the sake of clarity that some of these scripts contain a lot of 
extraneous code that is not used in the manuscript, such as code to simulate 
individual-level data, reverse causation, population stratification, etc.
This was used in internal testing but is now obsolete, and unfortunately it 
would be rather difficult to remove all of these parts of the code at this
stage. This is particularly true for `ukbb_run_sim_and_methods.R`, 
`sim_mvmr_sparse.R`, and `run_methods.R`.
