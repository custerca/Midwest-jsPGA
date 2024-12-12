library(tidyverse)
library(sf)
library(RColorBrewer)
############################# Merging predictions from all eras and GCM models
x <- list.files("Output/Predictions",pattern = ".*qs$")
eras <- unique(sapply(str_split(x,"_"),"[[",2))
gcms <- str_remove(unique(sapply(str_split(x,"_"),"[[",3)),".qs")

# lammean6 <- tibble()
# for(i in 1:length(eras)){
#   for(j in 1:length(gcms)){
#     y <- as_tibble(qs::qread(paste0("Output/Predictions/lammean_",eras[i],"_",gcms[j],".qs")),rownames="site_id")
#     lammean6 <- bind_rows(lammean6, y %>% mutate(era=eras[i],
#                                                GCM=gcms[j]))
#   }
# }
# rm(y)
# qs::qsave(lammean6,"Output/lammean6.qs")
# # 
# lammean <- lammean6 %>%
#   pivot_longer(-c(site_id,era,GCM),names_to = "spp", values_to = "lambda") %>%
#   group_by(site_id,era,spp) %>%
#   summarise(x = mean(lambda)) %>%
#   pivot_wider(names_from = era, values_from = x) %>%
#   ungroup()
# 
# qs::qsave(lammean,"Output/lammean.qs")
# 
# 
# 
# pctextinct6 <- tibble()
# for(i in 1:length(eras)){
#   for(j in 1:length(gcms)){
#     y <- as_tibble(qs::qread(paste0("Output/Predictions/pctextinct_",eras[i],"_",gcms[j],".qs")),rownames="site_id")
#     pctextinct6 <- bind_rows(pctextinct6, y %>% mutate(era=eras[i],
#                                                      GCM=gcms[j]))
#   }
# }
# rm(y)
# qs::qsave(pctextinct6,"Output/pctextinct6.qs")
# 
# pctextinct <- pctextinct6 %>%
#   pivot_longer(-c(site_id,era,GCM),names_to = "spp", values_to = "probex") %>%
#   group_by(site_id,era,spp) %>%
#   summarise(x = mean(probex)) %>%
#   pivot_wider(names_from = era, values_from = x) %>%
#   ungroup()
# 
# qs::qsave(pctextinct,"Output/pctextinct.qs")
# 
# rm(x)
############################# 

# simplifying to only fish data
# fishdat <- readRDS("data/MWfishz.rds")
# justfish <- select(fishdat,site_id,STATE,COMMON_NAME,GEAR,TOTAL_CATCH,EFFORT,lon,lat)
# qs::qsave(justfish,"data/justfish.qs")
rm(x)

fishdat <- qs::qread("data/justfish.qs")
lakedat <- readRDS("data/MWlakes_z.rds")

sppv <- unique(fishdat$COMMON_NAME)
spp_sites <- list() 
for(i in 1:length(sppv)){
  spp_sites[[sppv[i]]] <- fishdat %>%
    filter(COMMON_NAME==sppv[i],TOTAL_CATCH>0) %>%
    select(site_id) %>%
    distinct() %>%
    pull()
}


# lammean6 <- qs::qread("output/lammean6.qs")
# pctextinct6 <- qs::qread("output/pctextinct6.qs")

lammean <- qs::qread("output/lammean.qs")
pctextinct <- qs::qread("output/pctextinct.qs")

lakelocs <- lammean %>%
  rowwise() %>%
  mutate(smpl = (site_id %in% spp_sites[[spp]])*1L) %>%
  ungroup() %>%
  left_join(distinct(select(lakedat,site_id,lon,lat)),
            by="site_id") %>%
  select(site_id,spp,lon,lat,smpl,pred21=Current)

lakelocs %>%
  group_by(spp) %>%
  summarise(x=sum(smpl))

lammean %>%
  filter(!is.na(Current)) %>%
  select(site_id,spp,Current) %>%
  pivot_wider(names_from = spp,values_from = Current) %>%
  summary()


################################################################################### GCM summaries ###################################################################################
usa <-st_as_sf(maps::map("state",fill=TRUE,plot=FALSE))

gcm_sf <- st_as_sf(lakedat,coords=c("lon","lat"),crs=st_crs(4269)) %>%
  select(site_id,GCM,MC,LC) %>%
  pivot_longer(c(MC,LC),names_to="era",values_to = "tinc") %>%
  arrange(tinc) %>%
  filter(!is.na(tinc)) %>%
  mutate(era=factor(era,levels=c("MC","LC"),labels=c("Mid-century","Late-century")))

ggplot() +
    geom_sf(data=usa) +
    geom_sf(data=gcm_sf,aes(color=tinc)) +
    coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
    facet_grid(cols=vars(GCM), rows=vars(era)) +
    scale_color_viridis_c(na.value = "#5E5D61") +
    scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
    theme_bw() +
    labs(color="Increase \u00B0C") +
    theme(text=element_text(size=18),
          legend.text = element_text(size=18),
          strip.text = element_text(face = "bold"))

ggsave("Figures/gcmmap.jpeg",units="in",width=16,height=8,dpi=600)

MWmap_sf <- gcm_sf %>%
  filter(era=="Late-century") %>%
  group_by(site_id) %>%
  summarize(avgt=mean(tinc))


library(cowplot)
MWmap <- ggplot() +
  geom_sf(data=usa) +
  geom_sf(data=MWmap_sf,aes(color=avgt)) +
  coord_sf(xlim=c(-97.47815,-80),ylim=c(37.15983,50)) +
  scale_color_viridis_c(na.value = "#5E5D61") +
  scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
  theme_bw() +
  labs(color="Increase \u00B0C") +
  theme(text=element_text(size=18),
        legend.text = element_text(size=18),
        strip.text = element_text(face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.12),
        legend.background = element_rect(color="black"),
        panel.grid = element_blank())

ggdraw(MWmap) +
  draw_plot(
    {
      usa %>%
        mutate(fillr=ifelse(ID %in% c("iowa","indiana","illinois","michigan","minnesota","south dakota","wisconsin"),
                            "Y","N")) %>%
      ggplot() +
        geom_sf(aes(fill=fillr)) +
        scale_fill_manual(values = c("Y"="red","N"="white")) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank())
    },
    x=0.58,
    y=0.67,
    width=0.4,
    height=0.4
  )


ggsave("Figures/MWmap.jpeg",units="in",width=9,height=9,dpi=600)


st_drop_geometry(gcm_sf) %>%
  group_by(era,GCM) %>%
  summarise(avg=mean(tinc),
            min=min(tinc),
            max=max(tinc))

st_drop_geometry(gcm_sf) %>%
  group_by(era) %>%
  summarise(avg=mean(tinc),
            min=min(tinc),
            max=max(tinc))

################################################################################### Percent change in relative abundance ###################################################################################

lammean %>%
  mutate(Lpc=100*(LC-Current)/Current,
         Mpc=100*(MC-Current)/Current) %>%
  inner_join(lakelocs,by=c("site_id","spp")) %>%
  filter(!is.na(Lpc),pred21 > 0.01|smpl==1) %>%
  group_by(spp) %>%
  summarise(avg=mean(Lpc),
            med=median(Lpc),
            iqr25=quantile(Lpc,0.25),
            iqr75=quantile(Lpc,0.75),
            pd=sprintf("%.1f",mean(Lpc < 0)*100)) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")),
         tabtxt = paste0(sprintf("%.1f",med),
                         " (",
                         sprintf("%.1f",iqr25),
                         ", ",
                         sprintf("%.1f",iqr75),
                         ")")) %>%
  arrange(desc(avg)) %>%
  select(spp,avg,tabtxt,pd) %>%
  knitr::kable(format="latex",digits=2,booktabs=TRUE,linesep="")

lammean %>%
  mutate(Lpc=100*(LC-Current)/Current,
         Mpc=100*(MC-Current)/Current) %>%
  inner_join(lakelocs,by=c("site_id","spp")) %>%
  filter(!is.na(Lpc),pred21 > 0.01|smpl==1) %>%
  group_by(spp) %>%
  summarise(avg=mean(Mpc),
            med=median(Mpc),
            iqr25=quantile(Mpc,0.25),
            iqr75=quantile(Mpc,0.75),
            pd=sprintf("%.1f",mean(Mpc < 0)*100)) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")),
         tabtxt = paste0(sprintf("%.1f",med),
                         " (",
                         sprintf("%.1f",iqr25),
                         ", ",
                         sprintf("%.1f",iqr75),
                         ")")) %>%
  arrange(desc(avg)) %>%
  select(spp,avg,tabtxt,pd) %>%
  knitr::kable(format="latex",digits=2,booktabs=TRUE,linesep="")

lampc_LC <- lammean %>%
  mutate(pc=100*(LC-Current)/Current) %>%
  filter(!is.na(pc))

lampc_MC <- lammean %>%
  mutate(pc=100*(MC-Current)/Current) %>%
  filter(!is.na(pc))

lampc <- lammean %>%
  filter(!is.na(Current),!is.na(MC),!is.na(LC)) %>%
  mutate(LCpc=100*(LC-Current)/Current,
         MCpc=100*(MC-Current)/Current) %>%
  filter(!is.na(MCpc),!is.na(LCpc))

# lammean %>%
#   filter(spp %in% c("walleye","largemouth_bass"),site_id != "nhdhr_155417216") %>%
#   select(site_id,spp, Current) %>%
#   pivot_wider(names_from = "spp", values_from = "Current") %>%
#   ggplot() +
#   geom_point(aes(x=`largemouth_bass`,y=`walleye`)) +
#   theme_bw()
# 
# lammean %>%
#   filter(spp %in% c("walleye","largemouth_bass"),site_id != "nhdhr_155417216") %>%
#   select(site_id,spp, LC) %>%
#   pivot_wider(names_from = "spp", values_from = "LC") %>%
#   ggplot() +
#   geom_point(aes(x=`largemouth_bass`,y=`walleye`))

# expand.grid(unique(lampc_LC$spp),unique(lampc_LC$spp)) %>%
#   rename(spp=Var1,spp2=Var2) %>%
#   left_join(select(lampc_LC,site_id,spp,x=pc),by=c("spp"="spp"),relationship = "many-to-many")  %>%
#   left_join(select(lampc_LC,site_id,spp,y=pc),by=c("spp2"="spp","site_id"),relationship = "many-to-many") %>%
#   mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
#                                  "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
#                     labels=c("Cisco","Northern pike","Walleye","Black crappie",
#                              "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")),
#          spp2=factor(spp2,levels=c("cisco","northern_pike","walleye","black_crappie",
#                                    "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
#                      labels=c("Cisco","Northern pike","Walleye","Black crappie",
#                               "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
#   filter(site_id != "nhdhr_155417216") %>%
#   ggplot(aes(x=x,y=y)) +
#   geom_point() +
#   facet_grid(rows=vars(spp2),cols=vars(spp)) +
#   theme_bw()

lamera <- lampc %>%
  select(-c(Current,Historic,LC,MC)) %>%
  pivot_longer(-c(site_id,spp),names_to="era",values_to = "pc") %>%
  group_by(spp,era) %>%
  summarise(avg=mean(pc),
            l25=quantile(pc,0.25),
            u75=quantile(pc,0.75)) %>%
  ungroup() %>%
  mutate(era=factor(era,levels=c("MCpc","LCpc"),labels=c("Mid-century","Late-century")),
         spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")))


# ggplot(lamera, aes(x=era,y=avg,group=spp)) +
#   geom_point(aes(color=spp),position = position_dodge(width=0.1)) +
#   geom_line(aes(color=spp),position = position_dodge(width=0.1)) +
#   geom_errorbar(aes(ymin=l25,ymax=u75,color=spp),width = 0.25,position = position_dodge(width=0.1)) +
#   ggrepel::geom_text_repel(
#     data = lamera |> filter(era == "Mid-century"),
#     aes(
#       # label = paste0(spp,
#       #                  ": ",
#       #                  sprintf("%.1f",avg),
#       #                  "% (",
#       #                  sprintf("%.1f",l25),
#       #                  "%, ",
#       #                  sprintf("%.1f",u75),
#       #                  "%)")),
#       label = spp),
#     
#     segment.linetype=2,
#     segment.alpha=0.5,
#     size = 8 / .pt,
#     hjust = 1,
#     direction = "y",
#     nudge_x = -0.3
#   ) +
#   ggrepel::geom_text_repel(
#     data = lamera |> filter(era == "Late-century"),
#     aes(
#       # label = paste0(spp,
#       #                  ": ",
#       #                  sprintf("%.1f",avg),
#       #                  "% (",
#       #                  sprintf("%.1f",l25),
#       #                  "%, ",
#       #                  sprintf("%.1f",u75),
#       #                  "%)")),
#         label=spp),
#     segment.linetype=2,
#     segment.position=position_dodge(width=0.1),
#     segment.alpha=0.5,
#     size = 8 / .pt,
#     hjust = 0,
#     direction = "y",
#     nudge_x = 0.3
#   ) +
#   labs(x=NULL,y="% change",color=NULL) +
#   scale_color_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none",
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.y=element_blank())
# 
# 
# ggplot(lamera, aes(x=era,y=avg,group=spp)) +
#   geom_point(aes(color=spp),position = position_dodge(width=0.25)) +
#   geom_line(aes(color=spp),position = position_dodge(width=0.25)) +
#   geom_errorbar(aes(ymin=l25,ymax=u75,color=spp),
#                 width = 0.25,position = position_dodge(width=0.25)) +
#   labs(x=NULL,y="% change",color=NULL) +
#   scale_color_viridis_d() +
#   guides(color = guide_legend(reverse=T)) +
#   theme_bw() +
#   theme(panel.grid.major.y=element_blank(),
#         panel.grid.minor.y=element_blank())
# 
# 
# lampc_LC %>%
#   filter(spp %in% c("walleye","largemouth_bass"),site_id != "nhdhr_155417216") %>%
#   select(site_id,spp, pc) %>%
#   pivot_wider(names_from = "spp", values_from = "pc") %>%
#   ggplot() +
#   geom_point(aes(x=`largemouth_bass`,y=`walleye`)) +
#   labs(x="LMB percent change",y="Walleye percent change") +
#   theme_bw()
# 
# lampc_LC %>%
#   filter(spp %in% c("northern_pike","bluegill"),site_id != "nhdhr_155417216") %>%
#   select(site_id,spp, pc) %>%
#   pivot_wider(names_from = "spp", values_from = "pc") %>%
#   ggplot() +
#   geom_point(aes(y=`bluegill`,x=`northern_pike`)) +
#   labs(y="Bluegill percent change",x="Northern pike percent change") +
#   theme_bw()
# 
# lampc_LC %>%
#   filter(spp %in% c("yellow_perch","black_crappie"),site_id != "nhdhr_155417216") %>%
#   select(site_id,spp, pc) %>%
#   pivot_wider(names_from = "spp", values_from = "pc") %>%
#   ggplot() +
#   geom_point(aes(y=`black_crappie`,x=`yellow_perch`)) +
#   labs(y="Black crappie percent change",x="Yellow perch percent change") +
#   theme_bw()

lamtemp <- left_join(lammean,
                     select(lakedat,site_id,Currtemp=NLDAS2021,GCM,LCdiff=LC,MCdiff=MC),
                     by='site_id',relationship = "many-to-many") %>%
  mutate(LCtemp=Currtemp+LCdiff,
         MCtemp=Currtemp+MCdiff,
         LCpc=100*(LC-Current)/Current,
         MCpc=100*(MC-Current)/Current) %>%
  group_by(site_id,spp) %>%
  summarise(Current=mean(Current),
            MC=mean(MC),
            LC=mean(LC),
            MCpc=mean(MCpc),
            LCpc=mean(LCpc),
            MCdiff=mean(MCdiff),
            LCdiff=mean(LCdiff),
            Currtemp=mean(Currtemp),
            LCtemp=mean(LCtemp),
            MCtemp=mean(MCtemp)) %>%
  ungroup()

# lamtemp %>%
#   filter(site_id != "nhdhr_155417216") %>%
#   mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
#                                  "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
#                     labels=c("Cisco","Northern pike","Walleye","Black crappie",
#                              "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
#   ggplot(aes(x=LCdiff,y=LCpc,color=spp,group=spp)) +
#   geom_point(alpha=0.01) +
#   geom_smooth() +
#   scale_x_continuous(limits=c(3.2,5)) +
#   scale_y_continuous(limits=c(-101,100),transform = "pseudo_log") +
#   scale_color_viridis_d() +
#   theme_bw()

lamtempfull <- left_join(lammean %>%
                       pivot_longer(-c(site_id,spp),names_to = "era",values_to = "lambda") %>%
                         filter(era != "Historic"),
                     select(lakedat,site_id,Current=NLDAS2021,GCM,LCdiff=LC,MCdiff=MC) %>% 
                       mutate(LC=Current+LCdiff,
                              MC=Current+MCdiff) %>%
                       select(-c(LCdiff,MCdiff)) %>%
                       group_by(site_id) %>%
                       summarise(Current=mean(Current),
                                 MC=mean(MC),
                                 LC=mean(LC)) %>%
                       pivot_longer(-c(site_id),names_to = "era",values_to = "tmpc"),
                     by=c('site_id','era'),relationship = "many-to-many")


# lamtempfull %>%
#   filter(site_id != "nhdhr_155417216") %>%
#   mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
#                                  "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
#                     labels=c("Cisco","Northern pike","Walleye","Black crappie",
#                              "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
#   ggplot(aes(x=LCtemp,y=LC,color=spp,group=spp)) +
#   geom_point(alpha=0.01) +
#   geom_smooth() +
#   #scale_x_continuous(limits=c(3.2,5)) +
#   scale_y_continuous(transform = "pseudo_log") +
#   scale_color_viridis_d() +
#   theme_bw()


# lamtempfull %>%
#   mutate(tgrp = ifelse(tmpc <))
#   mutate(spp=factor(spp,
#                     levels=c("cisco","northern_pike","walleye","black_crappie",
#                                  "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
#                     labels=c("Cisco","Northern pike","Walleye","Black crappie",
#                              "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")),
#          era=factor(era,levels=c("Current","MC","LC"))) %>%
#   ggplot(aes(fill=spp, y=lambda, x=tmpc)) + 
#   geom_bin_2d(position="fill", stat="identity") +
#   facet_wrap(~era,nrow=1)
  

lammean_wide <- lammean %>%
  filter(site_id != "nhdhr_155417216") %>%
  select(site_id,spp,LC) %>%
  pivot_wider(names_from = spp,values_from = LC,values_fill = 0)

lamcols <- names(select(lammean_wide,-site_id))

lammax <- lammean %>%
  group_by(site_id) %>%
  reframe(Current=which.max(Current),
            MC=which.max(MC),
            LC=which.max(LC)) %>%
  rowwise() %>%
  mutate(Curspp=lamcols[Current],
         MCspp=lamcols[MC],
         LCspp=lamcols[LC],
         mid=ifelse(Current==MC,0,1),
         late=ifelse(MC==LC,0,1))
  
  
lampc_LC_sf <- inner_join(lakelocs,as_tibble(lampc_LC),by=c("site_id","spp")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(pc) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
  arrange(abs(pc))

lampc_LC_coord <- inner_join(lakelocs,as_tibble(lampc_LC),by=c("site_id","spp"))
lampc_LC_coord$latbin <- ntile(lampc_LC_coord$lat,n=30)

lampc_LC_coord %>%
  group_by(latbin,spp) %>%
  summarise(lat=mean(lat),
            avg=mean(pc),
            mdn=median(pc),
            sdp=sd(pc)) %>%
  ggplot() +
  geom_col(aes(x=lat,y=avg))
  


ggplot() +
  #geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
  geom_sf(data=filter(lampc_LC_sf, site_id != "nhdhr_155417216",pred21>0.01|smpl==1),
          aes(color=pc)) +
  coord_sf(xlim=c(-87.47815,-86.0606),y=c(45.5,46.5)) +
  facet_wrap(~spp,nrow=2) +
  scale_color_viridis_c(na.value = "#5E5D61") +
  scale_shape_manual(values=c("Y"=4,"N"=16),guide=NULL) +
  theme_bw() +
  labs(color="% change") +
  theme(text=element_text(size=18),
        strip.text = element_text(face = "bold"))
  


usa <-st_as_sf(maps::map("state",fill=TRUE,plot=FALSE))

mypal <- c("#fee090","#f46d43","#d73027","#a50026","#313695", "#4575b4","#74add1","#e6f598")
# mypal <- c("#fddbc7","#f4a582","#d6604d","#C41C0E","#2166ac", "#4393c3", "#92c5de", "#d1e5f0")
p <- {ggplot() +
    geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
    geom_sf(data=filter(lampc_LC_sf, site_id != "nhdhr_155417216",pred21>0.01|smpl==1),
            aes(color=pc)) +
    coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
    facet_wrap(~spp,nrow=2) +
    # scale_color_gradient2(low = "red",high = "blue",
    #                      mid = "gray",midpoint = 0,
    #                      breaks=c(-100,-50,0,50),
    #                      labels=function(x) paste0(sprintf("%.0f",x),"%")) +
    scale_color_viridis_c(na.value = "#5E5D61") +
    scale_shape_manual(values=c("Y"=4,"N"=16),guide=NULL) +
    scale_x_continuous(breaks=c(-95,-90,-85)) +
    scale_y_continuous(breaks=c(38,43,48)) +
    theme_bw() +
    labs(color="% change") +
    theme(text=element_text(size=18),
          strip.text = element_text(face = "bold"))}

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- mypal[k]
  k <- k+1
}
# jpeg("Figures/lampc_LC.jpeg", width = 13, height = 7, units = 'in', res = 600)
plot(g)
gg1 <- gridExtra::arrangeGrob(g)
ggsave("Figures/lampc_LC.jpeg", plot=gg1,width = 13, height = 7, units = 'in', dpi = 600)


lampc_MC_sf <- inner_join(lakelocs,as_tibble(lampc_MC),by=c("site_id","spp")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(pc) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
  arrange(abs(pc))

p2 <- {ggplot() +
    geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
    geom_sf(data=filter(lampc_MC_sf, site_id != "nhdhr_155417216",pred21>0.01|smpl==1),
            aes(color=pc)) +
    coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
    facet_wrap(~spp,nrow=2) +
    # scale_color_gradient2(low = "red",high = "blue",
    #                       mid = "gray",midpoint = 0,
    #                       breaks=c(-100,-50,0,50),
    #                       labels=function(x) paste0(sprintf("%.0f",x),"%")) +
    scale_color_viridis_c(na.value = "#5E5D61") +
    scale_shape_manual(values=c("Y"=4,"N"=16),guide=NULL) +
    scale_x_continuous(breaks=c(-95,-90,-85)) +
    scale_y_continuous(breaks=c(38,43,48)) +
    theme_bw() +
    labs(color="% change") +
    theme(text=element_text(size=18),
          strip.text = element_text(face = "bold"))}

g2 <- ggplot_gtable(ggplot_build(p2))
strip_both <- which(grepl('strip-', g2$layout$name))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- mypal[k]
  k <- k+1
}
plot(g2)
gg2 <- gridExtra::arrangeGrob(g2)
ggsave("Figures/lampc_MC.jpeg", plot=gg2,width = 13, height = 7, units = 'in', dpi = 600)



lammean %>%
  filter(Current > 0.01,site_id != "nhdhr_155417216",!is.na(LC)) %>%
  mutate(pLC=100*(LC-Current)/Current,
         pMC=100*(MC-Current)/Current) %>%
  select(spp,pLC,pMC) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill"))) %>%
  pivot_longer(-spp,names_to = "era",values_to = "pc") %>%
  mutate(era=factor(era,levels=c("pMC","pLC"),labels=c("Mid-century","Late-century"))) %>%
  ggplot() +
  ggridges::geom_density_ridges(aes(x=pc,y=spp,fill=spp),alpha=0.5) +
  facet_wrap(~era,ncol=2,scales="free_x") +
  scale_fill_manual(values=c("Yellow perch"="#fee090",
                             "Smallmouth bass"="#f46d43",
                             "Largemouth bass"="#d73027",
                             "Bluegill" ="#a50026",
                             "Cisco"="#313695", 
                             "Northern pike"="#4575b4",
                             "Walleye"="#74add1",
                             "Black crappie"="#e6f598")) +
  #scale_fill_viridis_d() +
  labs(y=NULL,fill=NULL,x="% change") +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size=20))

ggsave("Figures/pc_dens.jpeg", width = 13, height = 7, units = 'in', dpi = 600)


bg_pc <- lakelocs %>%
  filter(spp=="bluegill") %>%
  distinct() %>%
  left_join(filter(lammean,spp=="bluegill"),
            by=c("site_id","spp")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(LC)

bg_gg <- bg_pc %>%
  filter(pred21>0.01|smpl==1,site_id != "nhdhr_155417216") %>%
  select(site_id,MC,LC,Current) %>%
  mutate(pLC=(LC-Current)/Current,
         pMC=(MC-Current)/Current) %>%
  select(pLC,pMC) %>%
  pivot_longer(c(pMC,pLC),names_to = "era",values_to = "pc") %>%
  mutate(era=factor(era,levels=c("pMC","pLC"),labels=c("Mid-century","Late-century"))) %>%
  filter(!is.na(pc)) %>%
  arrange(abs(pc))

ggplot() +
  geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
  geom_sf(data=bg_gg,
          aes(fill=pc),
          size=3,shape=21) +
  facet_wrap(~era,nrow=1) +
  scale_fill_gradient2(low = "red",high = "blue",
                       mid = "white",midpoint = 0,
                       transform = "exp",
                       breaks=c(-.6,0,.3),
                       labels=function(x) paste0(sprintf("%.0f",x*100),"%")) +
  coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
  scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
  theme_bw() +
  labs(fill="% change") +
  theme(text=element_text(size=18),
        strip.text = element_text(face = "bold"))

################################################################################### Probability of extirpation ###################################################################################
probex_sf <- inner_join(lakelocs,as_tibble(pctextinct),by=c("site_id","spp")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(LC) %>%
  mutate(spp=factor(spp,levels=c("cisco","northern_pike","walleye","black_crappie",
                                 "yellow_perch","smallmouth_bass", "largemouth_bass","bluegill"),
                    labels=c("Cisco","Northern pike","Walleye","Black crappie",
                             "Yellow perch","Smallmouth bass", "Largemouth bass","Bluegill")))

st_drop_geometry(probex_sf) %>%
  filter(site_id != "nhdhr_155417216",
         pred21 > 0.01|smpl==1,
         !is.na(MC)) %>%
  select(spp,MC,LC) %>%
  pivot_longer(c(MC,LC),names_to = "era",values_to = "probex") %>%
  group_by(spp,era) %>%
  summarise(avg=mean(probex),
            med=median(probex),
            iqr25=quantile(probex,0.25),
            iqr75=quantile(probex,0.75),
            hiext=mean(probex>0.9),
            hiext_n=sum(probex>0.9)) %>%
  arrange(desc(era),spp) %>%
  mutate(tabtxt = paste0(sprintf("%.1f",med),
                         " (",
                         sprintf("%.1f",iqr25),
                         ", ",
                         sprintf("%.1f",iqr75),
                         ")"),
         hptab = paste0(sprintf("%.1f",hiext*100),"% (n=",hiext_n,")")) %>%
  select(era,spp,avg,tabtxt,hptab) %>%
  knitr::kable(format="latex",booktabs=TRUE,linesep="",digits=1)






usa <-st_as_sf(maps::map("state",fill=TRUE,plot=FALSE))

mypal <- c("#fee090","#f46d43","#d73027","#a50026","#313695", "#4575b4","#74add1","#e6f598")

p2 <- {ggplot() +
    geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
    geom_sf(data=filter(probex_sf,pred21 > 0.01|smpl==1),
            aes(color=LC)) +
    coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
    facet_wrap(~spp,nrow=2) +
    scale_color_viridis_c(na.value = "#5E5D61",begin = 0.1) +
    scale_shape_manual(values=c("Y"=4,"N"=16),guide=NULL) +
    scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
    theme_bw() +
    labs(color="Probability of \n extinction") +
    theme(text=element_text(size=18),
          strip.text = element_text(face = "bold"))}

g2 <- ggplot_gtable(ggplot_build(p2))
strip_both <- which(grepl('strip-', g2$layout$name))
k <- 1
for (i in strip_both) {
  j2 <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$fill <- mypal[k]
  k <- k+1
}
plot(g2)
gg2 <- gridExtra::arrangeGrob(g2)
ggsave("Figures/probex_LC.jpeg", plot=gg2,width = 13, height = 7, units = 'in', dpi = 600)


cisco_ext <- lakelocs %>%
  filter(spp=="cisco") %>%
  distinct() %>%
  left_join(filter(pctextinct,spp=="cisco"),
             by=c("site_id","spp")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(LC)

# minimum predicted relative abundance amongst sampled lakes where cisco observed ~ 0.02

st_bbox(filter(cisco_ext,pred21>0.02))

# cisco_gg <- bind_rows(cisco_ext %>%
#                         filter(smpl==1) %>%
#                         mutate(gg="Sampled"),
#                       cisco_ext %>%
#                         filter(pred21 > 0.029) %>%
#                         mutate(gg="Predicted")) %>%
#   mutate(gg=factor(gg,levels=c("Sampled","Predicted")))

cisco_gg <- cisco_ext %>%
  filter(pred21>0.01|smpl==1) %>%
  select(site_id,MC,LC) %>%
  pivot_longer(c(MC,LC),names_to = "era",values_to = "probex") %>%
  mutate(era=factor(era,levels=c("MC","LC"),labels=c("Mid-century","Late-century")))


cisco_LC <- ggplot() +
  geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
  geom_sf(data=cisco_gg,
          aes(fill=probex),
          size=4,shape=21) +
  coord_sf(xlim=c(-96.055,-83),y=c(39,49.1)) +
  facet_wrap(~era,nrow=1) +
  scale_fill_viridis_c(na.value = "#5E5D61",breaks=c(0,1)) +
  scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
  theme_bw() +
  labs(fill="Probability of extinction") +
  theme(text=element_text(size=18),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

cisco_LC
ggsave("Figures/cisco_ext.jpeg", plot=cisco_LC, width = 13, height = 7, units = 'in', dpi = 600)


#################################################### quick summaries for Gretchen ####################################################

lammean %>%
  filter(Current > 0.01) %>%
  filter(spp=="walleye") %>%
  mutate(pLC=100*(LC-Current)/Current,
         pMC=100*(MC-Current)/Current) %>%
  summarise(LC=mean(pLC,na.rm=TRUE),
            MC=mean(pMC,na.rm=TRUE))

lampc_LC %>%
  filter(site_id != "nhdhr_155417216") %>%
  filter(Current > 0.01) %>%
  group_by(spp) %>%
  summarise(min=min(pc),
            max=max(pc),
            med=median(pc),
            avg=mean(pc))

st_drop_geometry(probex_LC_sf) %>%
  filter(!is.na(LC)) %>%
  filter(pred21 > 0.21) %>%
  group_by(spp) %>%
  summarise(n=n(),
            avg=mean(LC)*100,
            n90=sum(LC >= 0.9),
            p90=mean(LC >= 0.9)*100)







