library(tidyverse)
library(sf)
library(data.table)
library(mwlaxeref)
############################################### initial data cleaning ############################################################
# Setting working directory
# setwd()

# List of all lakes for which we have future predicted temperature
laketemp <- arrow::read_feather("data/lake_temperature_metrics_GLM_NLDAS.feather",
                                col_select = c("site_id","year","mean_surf_jul"))

# Raw fish sample data
# Available on ScienceBase DOI:
rawfish <- read.csv("data/all_state_cpue_6Feb24.csv") %>%
  mutate(state=as.character(factor(state,levels=state.name,labels=state.abb)),
         date=as.Date(date,"%Y-%m-%d"),
         nhdhr_id=ifelse(!is.na(nhdhr_id),str_remove(nhdhr_id,"nhdhr_"),NA)) %>%
  mutate(nhdhr_id=ifelse(!is.na(nhdhr_id),paste0("nhdhr_",nhdhr_id),NA))

justIN <- filter(rawfish,state=="IN")

INids <- read.csv("data/IN_Lake_IDs.csv")

INdates <- read.csv("data/Indy_dates.csv")

INfull <- select(justIN,-lat_unspec,-lon_unspec,-date,-nhdhr_id) %>%
  left_join(distinct(select(INids,state_id,nhdhr_id,lon_unspec=lon,lat_unspec=lat)), by=c("lake_id"="state_id")) %>%
  left_join(distinct(select(INdates,total_effort_ident,olddate=date,month,year) %>%
                       mutate(total_effort_ident=as.character(total_effort_ident),
                              date=as.Date(ifelse(is.na(olddate),
                                                  paste0(month,"/01/",substr(year,3,4)),
                                                  olddate),
                                           "%m/%d/%y")) %>%
                       select(-olddate,-month,-year)), by="total_effort_ident")



wis_vert <- rawfish %>% 
  filter(state == "WI" & sampling_method == "vertical_gill_net") %>% 
  local_to_nhdhr(from_colname = "lake_id", states = "wi") %>% 
  select(-nhdhr_id) %>% 
  rename(nhdhr_id = nhdhr.id) %>%
  mutate(nhdhr_id=paste0("nhdhr_",nhdhr_id))
  
mwfish <- bind_rows(filter(rawfish,state!="IN"),INfull)  %>%
  filter(!(state == "WI" & sampling_method == "vertical_gill_net")) %>%
  bind_rows(wis_vert) %>%
  mutate(chkr=(nhdhr_id %in% unique(laketemp$site_id))*1L) 


# checking how many lakes per state don't have an NHD ID within the lake temperature dataset
notemp <- select(mwfish,state,nhdhr_id,chkr) %>% 
  filter(chkr==0) %>%
  distinct()
# numbers look good enough
group_by(notemp,state) %>% summarise(n=n())

#
mwfish_fix <- mwfish %>%
  select(state,
         lake_id,
         surv_id=total_effort_ident,
         site_id=nhdhr_id,
         date,
         sampling_method,
         effort=total_effort_1,
         e_unit=total_effort_1_units,
         species=species_1,
         count,
         lat=lat_unspec,
         lon=lon_unspec)

# Spatial object of state borders
state_sf <-st_read("data/stateshape/cb_2018_us_state_500k.shp") %>%
  filter(STUSPS %in% unique(mwfish$state))
st_write(state_sf,"data/MWstates.gpkg",delete_dsn = TRUE)


lakemeta <- read.csv("data/lake_metadata.csv") %>%
  filter(state %in% unique(mwfish$state)) %>%
  select(site_id,lat=centroid_lat,lon=centroid_lon)

write.csv(lakemeta,"data/NHDlaketemps.csv")

# Use state_sf to create spatial object for use in QGIS to get all NHD IDs in MW states
fish_noNHD <- mwfish_fix %>%
  filter(!is.na(lat)) %>%
  filter(is.na(site_id)) %>%
  select(state,lake_id,surv_id,lat,lon) %>%
  distinct() %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269))

# Use QGIS to link these lakes with nhdids
st_write(fish_noNHD,dsn="data/fish_noNHD.gpkg",delete_dsn = TRUE)

# Reading in spatially joined NHD IDs (joined with spatial object of lakes within predicted lake temp dataset)
NHDjoins <- st_read("data/MWNHD_joinedID.gpkg") %>%
  st_drop_geometry() %>%
  mutate(lat=as.numeric(feature_y),
         lon=as.numeric(feature_x),
         site_id=paste0("nhdhr_",permanent_identifier)) %>%
  select(state,lake_id,surv_id,site_id,lat,lon) %>%
  distinct()

mwfish_join <- mwfish_fix %>%
  left_join(NHDjoins,by=c("state","lake_id","surv_id")) %>%
  rowwise() %>%
  mutate(site_id=ifelse(is.na(site_id.x),site_id.y,site_id.x)) %>%
  ungroup() %>%
  select(state,site_id,lake_id,surv_id,date,sampling_method,effort,e_unit,species,count)

fullfish <- filter(mwfish_join,!is.na(site_id)) %>%
  left_join(lakemeta,by="site_id") %>%
  filter(!is.na(lat)) %>%
  mutate(gears=paste0(state,"_",sampling_method))

fullfish$gears = fct_recode(fullfish$gears,
                            "WI_FN" = "WI_fyke_net",
                            "WI_BS" = "WI_boom_shocker",
                            "WI_VGN" = "WI_vertical_gill_net",
                            "MN_GN" = "MN_Standard gill net sets",
                            "MN_TN2f" = "MN_Standard 3/4-in mesh, double frame trap net sets",
                            "MN_GNsh" = "MN_Standard gill nets, set shallow in stratified assessment",
                            "MN_VGN" = "MN_Standard Vertical Gillnet",
                            "MI_GN" = "MI_inland_gill_net",
                            "MI_FN" = "MI_large_mesh_fyke_net",
                            "MI_TN" = "MI_trap_net",
                            "MI_BS" = "MI_boomshocking",
                            "IL_EF_AC" = "IL_Boat electrofishing AC (Day)",
                            "IL_GN125" = "IL_Gill net - 125 ft experimental",
                            "IL_TN" = "IL_Trap net",
                            "IL_GN250" = "IL_Gill net - 250 ft experimental",
                            "IL_EF_DC" = "IL_Boat electrofishing DC (Day)",
                            "IL_TN05" = "IL_Trap net 0.5 in bar mesh",
                            "SD_FN34" = "SD_frame net (std 3/4 in)",
                            "SD_GN_exp" = "SD_std exp gill net",
                            "SD_GN_std" = "SD_AFS std gill net",
                            "SD_FN38" = "SD_std frame net (3/8 inch)",
                            "SD_FN" = "SD_AFS std frame net",
                            "SD_FN_Lrg" = "SD_large frame net",
                            "SD_FN_ST" = "SD_single throat frame net",
                            "SD_BS_night" = "SD_boat shocker (night)",
                            "SD_BS_DCnight" = "SD_boat shocker (night, DC)",
                            "SD_BS_day" = "SD_boat shocker (day)",
                            "SD_BS_ACnight" = "SD_boat shocker (night, AC)",
                            "SD_EF_snLMB" = "SD_spring night EF-LMB",
                            "IA_FN" = "IA_Fyke Netting Unspecified",
                            "IA_BS_2n" = "IA_DC Boat Shocker Day 2 Netters",
                            "IA_BS" = "IA_DC Boat Shocker Day",
                            "IA_GN" = "IA_Gill Net Unspecified",
                            "IN_GN_exp" = "IN_experimental gill net",
                            "IN_EF" = "IN_dc nighttime electrofishing",
                            "IN_GN_std" = "IN_gill net std exp",
                            "IN_TN" = "IN_standard trap net")


fishdat <- fullfish %>%
  mutate(SURVEYDATE = date,
         YEAR = year(SURVEYDATE),
         samp_id = paste0(state,"_",surv_id)) %>%
  select(STATE=state,site_id,samp_id,SURVEYDATE,YEAR,GEAR=gears,EFFORT=effort,COMMON_NAME=species,count,lon,lat)

saveRDS(fishdat,"data/fishdat_manip.rds")
########################################################## merging fish and lake data #############################################################

fishdat <- readRDS("data/fishdat_manip.rds") %>%
  filter(!is.na(SURVEYDATE)) %>%
  mutate(GEAR=as.character(GEAR))

# print(count(fishdat,GEAR) %>% arrange(n),n=40)
# print(
#   fishdat %>%
#     group_by(GEAR) %>%
#     summarise(n=length(unique(samp_id))) %>%
#     arrange(n),
#   n=100
# )

  
ggplot(fishdat,aes(x=YEAR)) + geom_histogram(binwidth = 1) + facet_wrap(~STATE,scales="free",nrow=3)

# Merging MN fish data with NHD ID crosswalk
# summarizing into one sample per species - gear - lake - sample date
setDT(fishdat)
MWfish <- fishdat[YEAR >= 2000 & YEAR <= 2021,
                  .(TOTAL_CATCH=sum(count),EFFORT=sum(EFFORT)),
                  by=.(site_id,STATE,samp_id,YEAR,COMMON_NAME,GEAR)
                  ]
#[GEAR %in% gearmin]

# setnames(MWfish,"species","COMMON_NAME")
# setnames(MWfish,"nhdid","site_id")

# Creating data table to include column for sample year and previous 4 years to create 5 year rolling means for lake temperature and secchi (water clarity)
MWyears <- unique(MWfish[,.(site_id,STATE,samp_id,YEAR)])[,.(YEAR5=seq(YEAR-4,YEAR,by=1)),
                                                             by=.(site_id,STATE,samp_id)]

# reading in NLDAS lake temperature data, only retaining columns of interest
# Corson-Dosch, H.R., Mcaliley, W.A., Platt, L.R.C., Padilla, J.A., and Read, J.S., 2023, 
# Daily water column temperature predictions for thousands of Midwest U.S. lakes between 1979-2022 
# and under future climate scenarios: U.S. Geological Survey data release, 
# https://doi.org/10.5066/P9EQQER7.

laketemp <- arrow::read_feather("data/lake_temperature_metrics_GLM_NLDAS.feather",
                                col_select = c("site_id","year","mean_surf_jul"))
setDT(laketemp)
# renaming column for merging
setnames(laketemp,"year","YEAR5")

# merging 5-year mean DT with temperature data
MWyear_temp <- data.table::merge.data.table(MWyears,
                                            laketemp,
                                            all.x=TRUE,
                                            by=c("site_id","YEAR5"))

# checking to make sure every lake has all 5 years of data
chkr1 <- MWyear_temp[,.(chkr=sum(!is.na(mean_surf_jul))),by=c("site_id","STATE","samp_id")]
count(chkr1,chkr)
# # some only have 0 or 4
yrtemp5 <- unique(chkr1[chkr==5,site_id])
# Cleanest data.  All lakes/surveys with missing data along the way (i.e., no site ID or surface temp) were filtered out
MWfish_temp <- data.table::merge.data.table(MWfish,
                                            # merging with a filtered DT to only lakes with all 5 years of lake temperatures
                                            MWyear_temp[site_id %in% yrtemp5,
                                                        # Mean temperature across all 5 years
                                                        .(july5yr=mean(mean_surf_jul)),
                                                        by=.(site_id,samp_id)],
                                            by=c("site_id","samp_id"))[
                                              # Removing NAs
                                              !is.na(july5yr)
                                              ]

n_distinct(MWfish_temp$site_id)
count(distinct(select(MWfish_temp,site_id,STATE)),STATE)

# Environmental and spatial data
MWlakes_z <- readRDS("data/MWlakes_z.rds") 

gcm_mods <- unique(MWlakes_z$GCM)

MWlakesall <- MWlakes_z %>%
  select(-MC) %>%
  pivot_wider(names_from = GCM,values_from = LC) %>%
  distinct()


setDT(MWlakesall)
setnames(MWlakesall,"state","STATE")

MWlakefish <- data.table::merge.data.table(MWfish_temp,
                                           MWlakesall,
                                           by=c("site_id","STATE"))

gearmin <- MWlakefish %>%
  select(samp_id,GEAR) %>%
  distinct() %>%
  count(GEAR) %>%
  arrange(n) %>%
  filter(n>30) %>%
  pull(GEAR)


MWfish_gear <- MWlakefish[GEAR %in% gearmin]

fullPsi <- as_tibble(qs::qread("data/MWlakes_Psi_ortho.qs"),rownames="site_id")
setDT(fullPsi)

MWfish_z <- data.table::merge.data.table(MWfish_gear,fullPsi,by="site_id")

saveRDS(MWfish_z,file="data/MWfishz.rds")


# create spatial object of all (sampled + unsampled) lakes
MWlakesall %>%
  mutate(fish=(site_id %in% unique(MWfish_z$site_id))*1L) %>%
  select(site_id,lon,lat,fish) %>%
  distinct() %>%
  st_as_sf(coords=c("lon","lat"),crs=st_crs(4269)) %>%
  st_write("data/spatial/MWlakes_pred.gpkg",delete_dsn = TRUE)

gc()
############################## stan data prep ################################

MWfish <- readRDS("data/MWfishz.rds")

M=32

# All data except TRC values
fishstan <- MWfish %>%
  select(site_id, samp_id, YEAR, COMMON_NAME, GEAR, EFFORT, TOTAL_CATCH, july5yr,
         ends_with(".z"),all_of(paste0("comp",1:M))) %>%
  select(-total.for.z) 


spp <- unique(fishstan$COMMON_NAME)
gears <- unique(fishstan$GEAR)

fishstan <- fishstan %>%
  pivot_wider(names_from = COMMON_NAME,values_from = TOTAL_CATCH,values_fill = 0) %>%
  pivot_wider(names_from = GEAR, values_from = EFFORT,values_fill = 0)

# stan data list
dat <- list()

dat$N <- nrow(fishstan)

dat$J <- length(gears)

dat$K <- length(spp)

dat$P <- ncol(select(fishstan,ends_with(".z"))) + 1
dat$alpha <- rep(1, dat$J)
dat$M <- M

dat$Y <- as.matrix(select(fishstan,all_of(spp)))
dat$y_vec <- c(dat$Y)

dat$E <- as.matrix(select(fishstan,all_of(gears)))

dat$X <- cbind(beta0=1,as.matrix(select(fishstan,ends_with(".z"))))

dat$temp <- fishstan$july5yr

dat$Psi <- as.matrix(select(fishstan,starts_with("comp")))
rownames(dat$Psi) <- paste0(fishstan$site_id)


# TRC Literature values

# CTmax
ctmax <- list()
# # black bullhead
# ctmax[["black bullhead"]] <- c(38.1,37.5,35:39)

# black_crappie
ctmax[["black_crappie"]] <- c(34.9,36:40,38:40)

# bluegill
ctmax[["bluegill"]] <- c(35.8,39.1,36.6,37.5,37.9,37.5,41.4,
                         35.6,38.5,37.9,36.8,38.8,40.4,36.7,
                         38,40.9,36.3,37.4,40.9,37,39.6,41.4,38.3)

# # brown bullhead
# ctmax[["brown bullhead"]] <- c(37.1,37.8,37.5,38)

# cisco species
ctmax[["cisco"]] <- c(26.2)

# largemouth_bass
ctmax[["largemouth_bass"]] <- c(37.80,36.70,40.10,36.30,35.40,36.70,38.50)

# northern_pike
ctmax[["northern_pike"]] <- c(33.6,33.25,30.8)

# smallmouth_bass
ctmax[["smallmouth_bass"]] <- c(36.90,36.30)

# walleye
ctmax[["walleye"]] <- c(34.40,34.80,35.00,34.30)

# # white sucker
# ctmax[["white sucker"]] <- c(34.9,31.6,35.1,36.1,32.7)

# # yellow bullhead
# ctmax[["yellow bullhead"]] <- c(36.4,38,37.9,35)

# yellow_perch
ctmax[["yellow_perch"]] <- c(33.50,35.00,33.4,34)

ctmax_mean <- sapply(ctmax,mean)
ctmax_sd <- sapply(ctmax,sd)
ctmax_sd["cisco"] <- mean(ctmax_sd,na.rm=TRUE)
ctmax_N <- sapply(ctmax,length)



# Topt
topt <- list()
# # black bullhead
# topt[["black bullhead"]] <- c(18:29,23:24)

# black_crappie
topt[["black_crappie"]] <- c(22:25)

# bluegill
topt[["bluegill"]] <- c(30,30,30,31,30.1,29,30,31,24,27,27,31.2)

# # brown bullhead
# topt[["brown bullhead"]] <- c(32,28.2,29.9)

# cisco species
topt[["cisco"]] <- c(18.1, 13:18) # Ty's code, C. artedii

# largemouth_bass
topt[["largemouth_bass"]] <- c(30.69,25.00,26.00,27.00,28.00,27.00,
                               30.00,23.90,27.00,25:30)

# northern_pike
topt[["northern_pike"]] <- c(19,21,20.9,19.8,26,18:25,19:21,21,26)

# smallmouth_bass
topt[["smallmouth_bass"]] <- c(28.00,25.00,26.00,25.00,29.00,27.00)

# walleye
topt[["walleye"]] <- c(22.1,25.2,22,22.6,22)

# # white sucker
# topt[["white sucker"]] <- c(27,24,26.9,26,24,26.9,24,24,24)

# # yellow bullhead
# topt[["yellow bullhead"]] <- c(28.3,28.8,27.6)

# yellow_perch
topt[["yellow_perch"]] <- c(22.5,23,24.2,23,28,23:24,29,23,26:30,24.7)

topt_mean <- sapply(topt,mean)
topt_sd <- sapply(topt,sd)
topt_N <- sapply(topt,length)

# CTmin
ctmin <- list()
# # black bullhead
# ctmin[["black bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# black_crappie
ctmin[["black_crappie"]] <- c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin from bluegill
                              3,3,5,7,10,11,15,6,11 # LILT from bluegill
)

# bluegill
ctmin[["bluegill"]] <-  c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin
                          3,3,5,7,10,11,15,6,11 # LILT
)

# # brown bullhead
# ctmin[["brown bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT

# cisco species
ctmin[["cisco"]] <- c(0.3, # CTmin
                      0,0.5,3,4.7 # LILT
)

# largemouth_bass
ctmin[["largemouth_bass"]] <- c(3.2,7.3, 10.7, # CTmin
                                5,7,11,5.5,11.8,10 # LILT
)

# northern_pike
ctmin[["northern_pike"]] <- c(0.1,5,3) #LILT

# smallmouth_bass
ctmin[["smallmouth_bass"]] <- c(2,4,4,7,10,10,2,4,7,10,10.1,1.6) # LILT

# walleye
ctmin[["walleye"]] <- c(0.1,5,3) #LILT from NP

# 

# # white sucker
# ctmin[["white sucker"]] <- c(2,3,6,6,2.5,6.6,4.8,6.1,4.8) #LILT
# 

# # yellow bullhead
# ctmin[["yellow bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# yellow_perch
ctmin[["yellow_perch"]] <- c(1.1, #CTmin
                             4,1.1,3.7,6.8 #LILT
)


ctmin_mean <- sapply(ctmin,mean)
ctmin_sd <- sapply(ctmin,sd)
ctmin_N <- sapply(ctmin,length)

# sigma
sigma_trc <- (topt_mean - ctmin_mean)/4

trcdf <- tibble(COMMON_NAME=spp) %>%
  left_join(tibble(COMMON_NAME=names(ctmax_mean),
                   ctmax=ctmax_mean,
                   sd_ctmax=ctmax_sd),
            by="COMMON_NAME") %>%
  left_join(tibble(COMMON_NAME=names(topt_mean),
                   topt=topt_mean,
                   sd_topt=topt_sd),
            by="COMMON_NAME") %>%
  left_join(tibble(COMMON_NAME=names(sigma_trc),
                   sigmaval=sigma_trc),
            by="COMMON_NAME")

# Number of iterations
nsim <- 100
ctmax_mat <- matrix(NA,nrow=nsim,ncol=dat$K)
topt_mat <- matrix(NA,nrow=nsim,ncol=dat$K)

mxtmp <- MWfish %>%
  filter(TOTAL_CATCH>0) %>%
  group_by(COMMON_NAME) %>%
  summarise(mx=max(july5yr))

trcdf2 <- left_join(trcdf,mxtmp,by="COMMON_NAME")
colnames(ctmax_mat) <- trcdf2$COMMON_NAME
colnames(topt_mat) <- trcdf2$COMMON_NAME

# loop
j=1L
tot=1
set.seed(1)
while(j <= nsim){
  print(tot)
  tot=tot+1L
  ctmax_vec <- rnorm(dat$K, mean = trcdf2$ctmax, sd = trcdf2$sd_ctmax)
  topt_vec <- rnorm(dat$K, mean = trcdf2$topt, sd = trcdf2$sd_topt)
  
  if(any(ctmax_vec - topt_vec < 1)) next
  # Testing how it looks if I make sure the fixed CTmax is never less than the observed max
  if(any(ctmax_vec < trcdf2$mx)) next
  
  #if(sum(round(ctmax_vec) == round(topt_vec)) > 0) next
  
  ctmax_mat[j,] <- ctmax_vec
  topt_mat[j,] <- topt_vec
  j=j+1L
  
}

dat$ctmax <- ctmax_mat
dat$topt <- topt_mat
dat$sigma <- trcdf2$sigmaval

saveRDS(dat,file=paste0("data/MWfish_standat.rds"))



################################# TRC Curves ################################# 
# TRC Literature values

# CTmax
ctmax <- list()
# # black bullhead
# ctmax[["black bullhead"]] <- c(38.1,37.5,35:39)

# black_crappie
ctmax[["black_crappie"]] <- c(34.9,36:40,38:40)

# bluegill
ctmax[["bluegill"]] <- c(35.8,39.1,36.6,37.5,37.9,37.5,41.4,
                         35.6,38.5,37.9,36.8,38.8,40.4,36.7,
                         38,40.9,36.3,37.4,40.9,37,39.6,41.4,38.3)

# # brown bullhead
# ctmax[["brown bullhead"]] <- c(37.1,37.8,37.5,38)

# cisco species
ctmax[["cisco"]] <- c(26.2)

# largemouth_bass
ctmax[["largemouth_bass"]] <- c(37.80,36.70,40.10,36.30,35.40,36.70,38.50)

# northern_pike
ctmax[["northern_pike"]] <- c(33.6,33.25,30.8)

# smallmouth_bass
ctmax[["smallmouth_bass"]] <- c(36.90,36.30)

# walleye
ctmax[["walleye"]] <- c(34.40,34.80,35.00,34.30)

# # white sucker
# ctmax[["white sucker"]] <- c(34.9,31.6,35.1,36.1,32.7)

# # yellow bullhead
# ctmax[["yellow bullhead"]] <- c(36.4,38,37.9,35)

# yellow_perch
ctmax[["yellow_perch"]] <- c(33.50,35.00,33.4,34)

ctmax_mean <- sapply(ctmax,mean)
ctmax_sd <- sapply(ctmax,sd)
ctmax_sd["cisco"] <- mean(ctmax_sd,na.rm=TRUE)
ctmax_N <- sapply(ctmax,length)



# Topt
topt <- list()
# # black bullhead
# topt[["black bullhead"]] <- c(18:29,23:24)

# black_crappie
topt[["black_crappie"]] <- c(22:25)

# bluegill
topt[["bluegill"]] <- c(30,30,30,31,30.1,29,30,31,24,27,27,31.2)

# # brown bullhead
# topt[["brown bullhead"]] <- c(32,28.2,29.9)

# cisco species
topt[["cisco"]] <- c(18.1, 13:18) # Ty's code, C. artedii

# largemouth_bass
topt[["largemouth_bass"]] <- c(30.69,25.00,26.00,27.00,28.00,27.00,
                               30.00,23.90,27.00,25:30)

# northern_pike
topt[["northern_pike"]] <- c(19,21,20.9,19.8,26,18:25,19:21,21,26)

# smallmouth_bass
topt[["smallmouth_bass"]] <- c(28.00,25.00,26.00,25.00,29.00,27.00)

# walleye
topt[["walleye"]] <- c(22.1,25.2,22,22.6,22)

# # white sucker
# topt[["white sucker"]] <- c(27,24,26.9,26,24,26.9,24,24,24)

# # yellow bullhead
# topt[["yellow bullhead"]] <- c(28.3,28.8,27.6)

# yellow_perch
topt[["yellow_perch"]] <- c(22.5,23,24.2,23,28,23:24,29,23,26:30,24.7)

topt_mean <- sapply(topt,mean)
topt_sd <- sapply(topt,sd)
topt_N <- sapply(topt,length)

# CTmin
ctmin <- list()
# # black bullhead
# ctmin[["black bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# black_crappie
ctmin[["black_crappie"]] <- c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin from bluegill
                              3,3,5,7,10,11,15,6,11 # LILT from bluegill
)

# bluegill
ctmin[["bluegill"]] <-  c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin
                          3,3,5,7,10,11,15,6,11 # LILT
)

# # brown bullhead
# ctmin[["brown bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT

# cisco species
ctmin[["cisco"]] <- c(0.3, # CTmin
                      0,0.5,3,4.7 # LILT
)

# largemouth_bass
ctmin[["largemouth_bass"]] <- c(3.2,7.3, 10.7, # CTmin
                                5,7,11,5.5,11.8,10 # LILT
)

# northern_pike
ctmin[["northern_pike"]] <- c(0.1,5,3) #LILT

# smallmouth_bass
ctmin[["smallmouth_bass"]] <- c(2,4,4,7,10,10,2,4,7,10,10.1,1.6) # LILT

# walleye
ctmin[["walleye"]] <- c(0.1,5,3) #LILT from NP

# 

# # white sucker
# ctmin[["white sucker"]] <- c(2,3,6,6,2.5,6.6,4.8,6.1,4.8) #LILT
# 

# # yellow bullhead
# ctmin[["yellow bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# yellow_perch
ctmin[["yellow_perch"]] <- c(1.1, #CTmin
                             4,1.1,3.7,6.8 #LILT
)


ctmin_mean <- sapply(ctmin,mean)
ctmin_sd <- sapply(ctmin,sd)
ctmin_N <- sapply(ctmin,length)

# sigma
sigma_trc <- (topt_mean - ctmin_mean)/4

# Thermal performance function 
TRC = function(temp, CTmax, Topt, sigma){
  
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  
  return(trc)
  
}

dat <- readRDS("data/MWfish_standat.rds")


ctmax_df <- bind_rows(ctmax_mean) %>%
  bind_rows(as_tibble(dat$ctmax)) %>%
  mutate(iter=row_number())

ctmax_l <- ctmax_df %>%
  pivot_longer(-iter,
               names_to = "Spp",
               values_to = "ctmax")


topt_df <- bind_rows(topt_mean) %>%
  bind_rows(as_tibble(dat$topt)) %>%
  mutate(iter=row_number())

topt_l <- topt_df %>%
  pivot_longer(-iter,
               names_to = "Spp",
               values_to = "topt")

sigma_df <- tibble(Spp=names(sigma_trc),
                   sigma=sigma_trc)

trc_plot <- tibble(iter=1:100) %>%
  mutate(t0=list(seq(1,43,by=0.5)),
         Spp=list(names(sigma_trc))) %>%
  unnest(cols=c(t0)) %>%
  unnest(cols=Spp) %>%
  left_join(ctmax_l,by=c("iter","Spp")) %>%
  left_join(topt_l,by=c("iter","Spp")) %>%
  left_join(sigma_df,by="Spp") %>%
  rowwise() %>%
  mutate(trc=TRC(temp = t0,CTmax = ctmax,Topt = topt,sigma = sigma)) %>%
  ungroup() %>%
  mutate(iter=paste0("iter",iter))

trc_plot$Spp <- factor(trc_plot$Spp,
                       levels=c("cisco","northern_pike",
                                "walleye","black_crappie",
                                "yellow_perch","smallmouth_bass",
                                "largemouth_bass","bluegill"),
                       labels=c("Cisco","Northern pike",
                                "Walleye","Black crappie",
                                "Yellow perch","Smallmouth bass",
                                "Largemouth bass","Bluegill"))


mypal <- c("#fddbc7","#f4a582","#d6604d","#C41C0E","#2166ac", "#4393c3", "#92c5de", "#d1e5f0")

toptgg <- tibble(Spp=names(topt_mean),
                 topt=topt_mean) %>%
  mutate(Spp = factor(Spp,
                      levels=c("cisco","northern_pike",
                               "walleye","black_crappie",
                               "yellow_perch","smallmouth_bass",
                               "largemouth_bass","bluegill"),
                      labels=c("Cisco","Northern pike",
                               "Walleye","Black crappie",
                               "Yellow perch","Smallmouth bass",
                               "Largemouth bass","Bluegill")))

p <- {ggplot() +
    geom_line(data=trc_plot[trc_plot$iter != "iter1",],
              aes(x=t0,y=trc,
                  color=Spp,
                  group=interaction(iter,Spp)
              ),
              alpha=0.2,linewidth=0.4) +
    geom_line(data=trc_plot[trc_plot$iter == "iter1",],
              aes(x=t0,y=trc,
                  color=Spp,
                  group=interaction(iter,Spp)
              ),
              alpha=2,size=2) +
    geom_vline(data=toptgg,aes(xintercept = topt),linetype=2) +
    facet_wrap(~Spp,nrow=2) +
    scale_color_viridis_d(guide=guide_legend(override.aes = list(linewidth = 10))) +
    labs(x="Temperature (\u00B0C)",y="Thermal Performance Scalar") +
    theme_bw() +
    theme(axis.title = element_text(size=22),
          axis.text = element_text(size=20),
          axis.text.x = element_text(angle=45,hjust=1),
          panel.grid = element_blank(),
          strip.text = element_text(size=26),
          legend.position = "none")}

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
k <- 1

for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- mypal[k]
  k <- k+1
}
plot(g)
gg1 <- gridExtra::arrangeGrob(g)
ggsave("Figures/TRCplot.jpeg", plot=gg1,width = 13, height = 7, units = 'in', dpi = 600)

