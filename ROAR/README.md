The files within this (ROAR) directory are those associated with fitting the jsPGA model on Penn State's ROAR computing cluster.

**submit_setup.sh** is a SLURM script to run the R script **jspga_setup.R**

**jspga_setup.R** is an R script that initializes the jsPGA model within stan using the R package `cmdstanr`. It also reads in the data object (MWfish_standat.rds) created within the Local directory and saves it as a new data object (MWfish_model.rds). This is a relic from previous scripts where this file created changes to the data object, which is no longer the case.

**submit_model.sh** is a SLURM script to run the R script **jspga_model.R**. It runs the R script across an array (currently 100) to account for uncertainty in the thermal tolerance parameters as described in the paper. 

**jspga_model.R** is an R script that runs a single iteration of the jsPGA model. Each iteration is unique based on the `jobid` as assigned by SLURM (via **submit_model.sh**). The script reads in the R data object (MWfish_model.rds) and performs the final data manipulations to fit the jsPGA model. The script then saves the model output from each iteration (mod_`jobid`.rds).
