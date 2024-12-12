# Spatial basis functions
library(tidyverse)
library(fields)
library(sf)
library(data.table)

#fishfull <- readRDS(file="data/MWfishfull.rds")

mwstate <- c("MN","MI","SD","WI","IA","IN","IL")

lakemeta <- read.csv("data/lake_metadata.csv") %>%
  filter(state %in% mwstate) %>%
  select(site_id, state, clarity,max_depth) %>%
  filter(!is.na(clarity),!is.na(max_depth))

setDT(lakemeta)

nldas_historic <- arrow::read_feather("Data/lake_temperature_metrics_GLM_NLDAS.feather",
                                col_select = c("site_id","year","mean_surf_jul")) %>%
  filter(year <= 2000) %>%
  group_by(site_id) %>%
  summarise(NLDAS2000=mean(mean_surf_jul)) %>%
  ungroup()

setDT(nldas_historic)

nldas_current <- arrow::read_feather("Data/lake_temperature_metrics_GLM_NLDAS.feather",
                                      col_select = c("site_id","year","mean_surf_jul")) %>%
  filter(year <= 2021 & year >= 2002) %>%
  group_by(site_id) %>%
  summarise(NLDAS2021=mean(mean_surf_jul)) %>%
  ungroup()

setDT(nldas_current)

temp_nldas <- bind_rows(nldas_current,nldas)

temp_gcm <- arrow::read_feather("Data/lake_temperature_metrics_GLM_GCM.feather",
                                col_select = c("site_id","GCM","year","mean_surf_jul")) %>%
  mutate(time.period = ifelse(year <= 2000,"historic",
                              ifelse(year >= 2080,"late-century","mid-century"))) %>%
  group_by(site_id,GCM,time.period) %>%
  summarise(mean.GCM=mean(mean_surf_jul)) %>%
  ungroup()

setDT(temp_gcm)

# Difference between historical temps and future temps
gcm_diff <- temp_gcm %>%
  pivot_wider(names_from = "time.period",values_from = "mean.GCM") %>%
  mutate(LC=`late-century`-`historic`,
         MC=`mid-century`-`historic`) %>%
  select(site_id,GCM,LC,MC)

MWlakes <- data.table::merge.data.table(data.table::merge.data.table(data.table::merge.data.table(lakemeta,
                                                                                                  nldas_historic,
                                                                                                  by=c("site_id")),
                                                                     nldas_current,
                                                                     by=c("site_id")),
                                        gcm_diff, by=c("site_id"))

gc()
# All LAGOS data obtained from 
# LAGOS-US Geo https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1136.3
# LAGOS-US Locus https://portal.edirepository.org/nis/mapbrowse?packageid=edi.854.1
# LAGOS land use data
lagos_lc <- read.csv("Data/lagos geo/zone_landuse.csv")
setDT(lagos_lc)
# watershed spatial division
lagos_ws <- lagos_lc[spatial_division=="ws"]
rm(lagos_lc)

lagos_info <- read.csv("Data/lagos locus/lake_information.csv")
setDT(lagos_info)

lagos_nhd <- lagos_info[,site_id:=paste0("nhdhr_",lake_nhdid)
                        # filtering to sites within MN fish data
                        ][site_id %in% unique(MWlakes$site_id)
                          # selecting columns and mutating data types
                          ][,.(lagoslakeid=as.character(lagoslakeid),
                               zoneid=as.character(ws_zoneid),
                               lake_elevation_m,
                               lon=lake_lon_decdeg,
                               lat=lake_lat_decdeg),
                            by=site_id]

lagos_char <- read.csv("Data/lagos locus/lake_characteristics.csv")
setDT(lagos_char)

# selecting columns and mutating data types
lagos_char_2 <- lagos_char[,.(lagoslakeid=as.character(lagoslakeid),lake_waterarea_ha)]

lagos_nhd_2 <- data.table::merge.data.table(lagos_nhd,
                                            lagos_char_2,
                                            by="lagoslakeid")


# merging all LAGOS data
lagos_landuse <- merge(lagos_nhd_2,
                       lagos_ws,
                       by="zoneid")[,
                                    .(lon,lat,
                                      elevation=lake_elevation_m,
                                      lakearea=lake_waterarea_ha,
                                      total_dev = sum(c(nlcd_devopen21_pct, nlcd_devlow22_pct,
                                                        nlcd_devmed23_pct, nlcd_devhi24_pct), na.rm=T),
                                      total_ag = sum(c(nlcd_past81_pct, nlcd_cultcrop82_pct), na.rm=T),
                                      total_for = sum(nlcd_fordec41_pct,nlcd_forcon42_pct,nlcd_formix43_pct,na.rm=T)),
                                    by=c("site_id","year")][year==2011]

setnames(lagos_landuse,"year","NLCDyear")


MWlakesall <- data.table::merge.data.table(MWlakes,
                                           lagos_landuse,
                                           by=c("site_id"))

MWlakesall$secchi <- 1.7/MWlakesall$clarity


saveRDS(MWlakesall,file="data/MWlakesall.rds")


# selecting environmental covariates then transforming and standardizing them
zmat <- unique(MWlakesall[,.(site_id,
                              secchi,
                              max_depth,
                              elevation,
                              lakearea,
                              total_dev,
                              total_ag,
                              total_for)])[,.(site_id,
                                              secchi.z=as.numeric(scale(secchi)),
                                              elevation.z=as.numeric(scale(elevation)),
                                              depth.z=as.numeric(scale(log(max_depth))),
                                              lakearea.z=as.numeric(scale(log(lakearea))),
                                              total.dev.z=as.numeric(scale(car::logit(total_dev/100,adjust = 0.01))),
                                              total.ag.z=as.numeric(scale(car::logit(total_ag/100,adjust = 0.01))),
                                              total.for.z=as.numeric(scale(car::logit(total_for/100,adjust = 0.01))))]

zmat %>%
  distinct() %>%
  pivot_longer(ends_with(".z"),names_to = "covs",values_to = "vals") %>%
  ggplot() +
  geom_density(aes(x=vals,color=covs))

MWlakes_z <- data.table::merge.data.table(MWlakesall[,.(site_id,
                                                        state,
                                                        NLDAS2000,
                                                        NLDAS2021,
                                                        GCM,
                                                        LC, 
                                                        MC,
                                                        lon,
                                                        lat)],
                                         zmat,
                                         by=c("site_id"))

saveRDS(MWlakes_z,file="data/MWlakes_z.rds")















