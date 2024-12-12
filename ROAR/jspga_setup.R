library(cmdstanr,lib.loc="/storage/home/cac6877/R")

# Read in list object (simulated or case study)
dat <- readRDS("MWfish_standat.rds")


#identifying stan file path
stanfile <- "jspga.stan"

#creating stan model
mod <- cmdstan_model(stanfile)

saveRDS(dat,"MWfish_model.rds")

