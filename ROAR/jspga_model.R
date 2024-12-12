library(cmdstanr,lib.loc="/storage/home/cac6877/R")
library(tidyverse,lib.loc="/storage/home/cac6877/R")
# This code works by bringing in the jobid from the Roar supercomputer
# Each model realization is indexed by jobid to use a different randomly sampled value for CTmax and Topt
# Output from each model is saved and then merged into single posterior within analysis scripts

jobid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

# Read in list object (simulated or case study)
dat <- readRDS("MWfish_model.rds")

dat$CTmax_prior = tibble(species=dat$species) %>%
  left_join(dat$kcode,by="species") %>%
  left_join(tibble(ctmax=dat$ctmax[jobid,],
                   COMMON_NAME=names(dat$ctmax[jobid,])),
            by="COMMON_NAME") %>%
  pull(ctmax)

dat$Topt_prior = tibble(species=dat$species) %>%
  left_join(dat$kcode,by="species") %>%
  left_join(tibble(topt=dat$topt[jobid,],
                   COMMON_NAME=names(dat$topt[jobid,])),
            by="COMMON_NAME") %>%
  pull(topt)

dat$ctmax <- NULL
dat$topt <- NULL

# Initial values function for parameters that were problematic
init_fun <- list(#theta=matrix(rep(1/dat$J,dat$K * dat$J),nrow=dat$K),
  recip_phi=1, 
  beta = matrix(1,nrow=dat$P,ncol=dat$K),
  z_A = matrix(1,nrow=dat$K,ncol=dat$M),
  tau = rep(1,dat$K),
  L_Sigma = diag(rep(1,dat$K))
)

#identifying stan file path
stanfile <- "jspga.stan"
#creating stan model
mod <- cmdstan_model(stanfile)

fit <- mod$sample(data = dat, 
                  seed=10,
                  init = list(init_fun),
                  chains = 1,
                  thin=5,
                  max_treedepth = 11,
                  iter_warmup = 1000,
                  iter_sampling = 1000,
                  refresh = 100)

gc()
# converting cmdstanr output into stanfit object to utilize familiar functions for posterior merging and analysis
stancmd <- rstan::read_stan_csv(fit$output_files())

# saving stanfit object
saveRDS(stancmd,file=paste0("Output/MW/mod_",jobid,".rds"))



