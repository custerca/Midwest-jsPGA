# Spatial basis functions
library(tidyverse)
library(fields)
library(sf)
library(data.table)

# From MW_LakePrep.R
MWlakesall <- readRDS("data/MWlakes_z.rds")


locs <- MWlakesall %>%
  select(site_id,lon,lat) %>%
  arrange(site_id) %>%
  distinct() %>% 
  select(lon, lat) %>% 
  as.matrix()

dmat <- rdist.earth(locs, miles = F)
diag(dmat) <- 0

C = fields::Matern(dmat,range = max(dmat)/6, smoothness = 2.5)

gc() 

dcomp <- svd(C)

site_order <- MWlakesall %>%
  select(site_id,lon,lat) %>%
  arrange(site_id) %>%
  distinct() %>% 
  select(site_id)

decomp_df <- as_tibble(dcomp$u %*% diag(sqrt(dcomp$d)), 
                       .name_repair = ~paste0("comp", 1:nrow(locs))) %>% 
  mutate(lon = locs[,1], 
         lat = locs[,2],
         site_id = site_order$site_id)

saveRDS(decomp_df,"data/MWlakes_sbv.rds")

###############################################

gc()
MWlakesall <- readRDS("data/MWlakes_z.rds")
decomp_df <- readRDS("data/MWlakes_sbv.rds")

site_order <- select(decomp_df,site_id) %>%
  arrange(site_id)
  
psitemp <- left_join(site_order,decomp_df,by="site_id") %>%
  distinct() %>%
  column_to_rownames("site_id") %>%
  select(-lon,-lat) %>%
  as.matrix()

rm(decomp_df)

xmat <- left_join(site_order,MWlakesall %>%
                     select(-total.for.z) %>%
                     select(site_id,ends_with(".z")) %>%
                     distinct(),
                   by="site_id") %>%
  column_to_rownames("site_id") %>%
  as.matrix()

rm(MWlakesall)
gc()

# (I - P) where P=X * inv (X’X) * X’

Pmat = xmat %*% solve(t(xmat) %*% xmat) %*% t(xmat)
Imat = diag(rep(1,nrow(Pmat)))

fullPsi <- (Imat - Pmat) %*% psitemp

qs::qsave(psitemp,"data/MWlakes_Psi_nonortho.qs")
qs::qsave(fullPsi,"data/MWlakes_Psi_ortho.qs")









