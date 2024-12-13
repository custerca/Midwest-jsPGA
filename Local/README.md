The files within this directory are the R scripts used on the local machine to prepare the data for model fitting, and the R scripts used to analyze the output from the completed model fitting process. Associated data objects for these scripts are too large and not available on Github.

Below are descriptions of the role of each R script found within this directory.

**MW_LakePrep.R** prepares the lake specific data for analysis. This includes all associated habitat data and predicted lake surface temperature data from the GCM models. The file produces two R data objects (.rds) that contain the relevant habitat and lake temperature data for each lake included in the analysis. The first data object (MWlakesall.rds) has the raw habitat data values and the second (MWlakes_z.rds) has the standardized (mean 0 and standard deviation 1) data values.

**MW_SBV** creates the spatial basis vectors for each lake within the study.  This script requires the output from **MW_LakePrep.R** (MWlakes_z.rds). This script produces two R data objects: the raw spatial basis vectors calculated directly from the locations (MWlakes_Psi_nonortho.qs) and the orthogonal (with respect to the habitat effects) basis vectors (MWlakes_Psi_ortho.qs).

**MW_FishPrep.R** prepares the fish and lake data for fitting the jsPGA model. This script requires the output from **MW_LakePrep.R** (MWlakes_z.rds) and **MW_SBV** (MWlakes_Psi_ortho.qs). This script produces an R data object (MWfish_standat.rds) that is formatted for fitting the jsPGA model using Stan within R. This script also includes code to create species specific thermal response curve figures (redundant with post-fit analysis scripts).

**MW_predictions.R** calculates the predictions across the entire study region from the model fitting output. This script requires each unique model fit object (i.e., each of the output files created from the numerical integration approach performed on the ROAR supercomputer) to crate a single posterior distribution for each of the parameters.  The script creates a number of different data objects and prediction summaries that are required for the analysis and summary scripts found within this directory.

**MW_Analysis.R** creates all the figures and tables that summarize the jsPGA parameter estimates. This script also produces the map of the raw values for each of the 6 habitat covariates and the mean catch values for each species for each sampled lake (where catch is greater than 0 for that species).

