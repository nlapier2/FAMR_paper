# Simulation scripts for the FAMR paper

This folder contains scripts used to perform simulations and generate results
for the FAMR paper.

Phenotypes were simulated according to the DAG in the FAMR paper 
(see the paper linked on the top-level README).
The phenotype simulation is performed by `sim_mvmr_sparse.R`.

FAMR is run by the [FAMR package](https://github.com/nlapier2/FAMR/).
MR-BMA code is in the `mr-bma.R` script, and code for all other methods is
in the `methods.R` script. `run_methods.R` is a wrapper script to run a list of
these methods, adding factors as extra exposures if specified by the user.

`ukbb_run_sim_and_methods.R` is the big wrapper script around all of this
that takes in input data, simulates phenotypes using it, and runs all methods.

`utils.R` is a utility script that contains common functions used by many of the
scripts above.

The results in the paper are averaged over many simulation replicates.
In practice we did this by running array jobs on our compute cluster, with one
job per simulation, and each simulation writing separate results files.
The scripts to do this, which also give the exact paramter settings used,
are in the `submit/` folder.

After all simulation replicates were run, `eval_all_sim_res.R` was used to
summarize simulation results over all replicates and compute metrics.

For more details on each of these, please peruse the code.
For examples of how to run them, please see the `submit/` folder.

We note for the sake of clarity that these scripts contain a lot of extraneous 
code that is not used in the manuscript, such as code to simulate 
individual-level data, reverse causation, population stratification, etc.
We note this so users will not be confused by it, but leave it in for 
completeness and generality.
