library(rstan)
library(tidyverse)
library(data.table)


TRC = function(temp, CTmax, Topt, sigma){
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  return(trc)
}

# # vector of all file names within folder
# x <- list.files("output/MW/")
# simnum <- sapply(strsplit(str_remove(x,".rds"),"_"),"[[",2)
# 
# # reading in each of the RDS files
# out_list <- list()
# for(i in 1:length(x)){
#   name <- paste0("iter",simnum[i])
#   out_list[[name]] <- readRDS(paste0("output/MW/",x[i]))
#   print(i)
# }
# # merging all stan objects (each model realization posterior) into single object (single posterior for analysis)
# out <- sflist2stanfit(out_list)
# rm(out_list)

# saving single stan object containing all model posteriors
# qs::qsave(out,file="output/MWpost.qs")
# 
# out <- qs::qread("output/MWpost.qs")
# chains <- rstan::extract(out,permuted=TRUE,
#                          pars=c("beta",
#                                 "theta",
#                                 "phi",
#                                 "A",
#                                 "tau",
#                                 "Sigma",
#                                 "Pi"))
# 
# # beta
# dimnames(chains$beta)[[2]] <- colnames(dat$X)
# dimnames(chains$beta)[[3]] <- colnames(dat$Y)
# # theta
# dimnames(chains$theta)[[2]] <- colnames(dat$Y)
# dimnames(chains$theta)[[3]] <- colnames(dat$E)
# # A
# dimnames(chains$A)[[2]] <- paste0("A",1:dat$M)
# dimnames(chains$A)[[3]] <- colnames(dat$Y)
# #tau
# dimnames(chains$tau)[[2]] <- colnames(dat$Y)
# #Sigma
# dimnames(chains$Sigma)[[2]] <- colnames(dat$Y)
# dimnames(chains$Sigma)[[3]] <- colnames(dat$Y)
# #Pi
# dimnames(chains$Pi)[[2]] <- colnames(dat$Y)
# dimnames(chains$Pi)[[3]] <- colnames(dat$Y)
# 
# qs::qsave(chains,"output/MWchains.qs")

out <- qs::qread("output/MWpost.qs")
standat <- readRDS("data/MWfish_standat.rds")

sigmadf <- tibble(sigma = standat$sigma,
                  spp = names(standat$sigma)) %>%
  distinct()

chains <- qs::qread("output/MWchains.qs")

lakedat <- readRDS("data/MWlakes_z.rds")

fullPsi <- as_tibble(qs::qread("data/MWlakes_Psi_ortho.qs"),rownames="site_id")

site_order <- lakedat %>%
  select(site_id) %>%
  distinct() %>%
  arrange(site_id)

Psimat <- fullPsi %>%
  select(site_id,paste0("comp",1:standat$M)) %>%
  distinct() %>%
  right_join(site_order,by = join_by(site_id)) %>%
  column_to_rownames("site_id") %>%
  as.matrix()

# number of unique lakes
n_lakes <- nrow(Psimat)
# species
spp <- colnames(standat$Y)

betachain = chains$beta
dim(betachain)
#20000 x 7 x8
covs <- dimnames(chains$beta)[[2]]

Xmat <- lakedat %>%
  mutate(beta0=1) %>%
  select(site_id,all_of(covs)) %>%
  distinct() %>%
  right_join(site_order,by = join_by(site_id)) %>%
  column_to_rownames("site_id") %>%
  as.matrix()

# betadf <- tibble()
# for(i in 1:8){
#   sppi <- spp[i]
#   x <- as_tibble(betachain[,,i]) %>%
#     rownames_to_column("iter") %>%
#     mutate(iter=as.numeric(iter),
#            spp=spp[i]) %>%
#     pivot_longer(-c(spp,iter),names_to = "covs",values_to = "beta_est")
# 
#   betadf <- bind_rows(betadf,x)
# 
# }
# 
# g <- betadf %>%
#   filter(covs != "beta0") %>%
#   ggplot(aes(x=iter,y=beta_est)) +
#   geom_line() +
#   facet_grid(rows=vars(spp),cols=vars(covs)) +
#   theme_bw()
# 
# ggsave(filename = "Figures/betachains.jpeg", plot=g,
#        units="in",width = 10,height = 10,dpi = 600)
# 
# g0 <- betadf %>%
#   filter(covs == "beta0") %>%
#   ggplot(aes(x=iter,y=beta_est)) +
#   geom_line() +
#   facet_grid(rows=vars(spp),cols=vars(covs)) +
#   theme_bw()
# 
# ggsave(filename = "Figures/beta0.jpeg", plot=g0,
#        units="in",width = 10,height = 10,dpi = 600)


#number of iterations across all chains
n_samps = dim(betachain)[1]

# A - basis coefficient matrix
Achain <- chains$A
dim(Achain)
# 20000 x 32 x 8

n_chains <- length(out@sim$samples) # number of chains ran
c_ind <- sort(rep(1:n_chains, n_samps/n_chains))



t_mods <- unique(lakedat$GCM)
t_eras <- c("LC","MC","Current","Historic")
temp_list <- list()
temp_list[["LC"]] <- lapply(t_mods,function(x) 
  lakedat %>% filter(GCM==x) %>% mutate(temp = NLDAS2000 + LC) %>% arrange(site_id) %>% select(site_id,temp) %>% pull(temp,name=site_id)
)
names(temp_list[["LC"]]) <- t_mods

temp_list[["MC"]] <- lapply(t_mods,function(x) 
  lakedat %>% filter(GCM==x) %>% mutate(temp = NLDAS2000 + MC) %>% arrange(site_id) %>% select(site_id,temp) %>% pull(temp,name=site_id)
)
names(temp_list[["MC"]]) <- t_mods

temp_list[["Current"]] <- lapply(t_mods,function(x) 
  lakedat %>% filter(GCM==x) %>% mutate(temp = NLDAS2021) %>% arrange(site_id) %>% select(site_id,temp) %>% pull(temp,name=site_id)
)
names(temp_list[["Current"]]) <- t_mods

temp_list[["Historic"]] <- lapply(t_mods,function(x) 
  lakedat %>% filter(GCM==x) %>% mutate(temp = NLDAS2000) %>% arrange(site_id) %>% select(site_id,temp) %>% pull(temp,name=site_id)
)
names(temp_list[["Historic"]]) <- t_mods


x = Sys.time()

for(t in 4:4){
  teraloop <- t_eras[t]
  for(z in 1:6){
    lamdf <- array(0,dim=c(n_samps,n_lakes,standat$K))
    dimnames(lamdf)[[2]] <- rownames(Psimat)
    dimnames(lamdf)[[3]] <- spp
    
    scadf <- lamdf
    
    tmodloop <- t_mods[z]
    
    for(i in 1:n_samps){
      for(j in 1:length(spp)){
        spploop <- spp[j]
        scadf[i,,spploop] = TRC(temp_list[[teraloop]][[tmodloop]],
                                standat$ctmax[c_ind[i],spploop],
                                standat$topt[c_ind[i],spploop],
                                sigma=filter(sigmadf,spp==spploop)$sigma)
      }
      
      lamdf[i,,] = exp((Xmat %*% betachain[i,,]) + (Psimat %*% Achain[i,,])) * scadf[i,,]
      
      if(((i %% 100)==0)|i==1) print(paste0(teraloop,"-",tmodloop,"-",i))
      
    }
    qs::qsave(lamdf,file=paste0("Output/Predictions/Arrays/",teraloop,"preds_lambda_",tmodloop,".qs"),preset="fast")
    qs::qsave(scadf,file=paste0("Output/Predictions/Arrays/",teraloop,"preds_scalar_",tmodloop,".qs"),preset="fast")
    
    # Mean predicted relative abundance across for each lake
    lammean <- apply(lamdf,c(2,3),mean)
    qs::qsave(lammean,file=paste0("Output/Predictions/lammean_",teraloop,"_",tmodloop,".qs"))
    rm(lammean)
    rm(lamdf)
    
    # Average percent extinction for each lake
    pct_extinct <- apply(scadf,c(2,3),function(x) mean(x==0))
    qs::qsave(pct_extinct,file=paste0("Output/Predictions/pctextinct_",teraloop,"_",tmodloop,".qs"))
    rm(pct_extinct)
    rm(scadf)
    
    gc()
  }
}

Sys.time() - x





x = Sys.time()
for(z in 1:6){
  lamdf <- array(0,dim=c(n_samps,n_lakes,standat$K))
  dimnames(lamdf)[[2]] <- rownames(Psimat)
  dimnames(lamdf)[[3]] <- spp
  
  scadf <- lamdf
  
  tmodloop <- t_mods[z]
  
  for(i in 1:n_samps){
    for(j in 1:length(spp)){
      spploop <- spp[j]
      scadf[i,,spploop] = TRC(temp_list$Current[[tmodloop]],
                              standat$ctmax[c_ind[i],spploop],
                              standat$topt[c_ind[i],spploop],
                              sigma=filter(sigmadf,spp==spploop)$sigma)
    }
    
    lamdf[i,,] = exp((Xmat %*% betachain[i,,]) + (Psimat %*% Achain[i,,])) * scadf[i,,]
    
    if(((i %% 100)==0)|i==1) print(paste0(t_mods[z],"-",i))
    
  }
  qs::qsave(lamdf,file=paste0("Output/Predictions/Arrays/LCpreds_lambda_",tmodloop,".qs"),preset="fast")
  qs::qsave(scadf,file=paste0("Output/Predictions/Arrays/LCpreds_scalar_",tmodloop,".qs"),preset="fast")
  
  # Mean predicted relative abundance across for each lake
  lammean_LC <- apply(lamdf,c(2,3),mean)
  qs::qsave(lammean_LC,file=paste0("Output/Predictions/lammean_LC_",tmodloop,".qs"))
  rm(lammean_LC)
  
  # Average percent extinction for each lake
  pct_extinct_LC <- apply(scadf,c(2,3),function(x) mean(x==0))
  qs::qsave(pct_extinct_LC,file=paste0("Output/Predictions/pctextinct_LC_",tmodloop,".qs"))
  rm(pct_extinct_LC)
  
  gc()
}
Sys.time() - x





