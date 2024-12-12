library(rstan)
library(tidyverse)
library(kableExtra)
library(scoringRules)
library(gridExtra)

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
dat <- readRDS("data/MWfish_standat.rds")
# 
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
chains <- qs::qread("output/MWchains.qs")

lapply(chains,dim)

sppv <- c("cisco","northern_pike",
          "walleye","black_crappie",
          "yellow_perch","smallmouth_bass",
          "largemouth_bass","bluegill")

# quick function to summarize parameters from posterior chains
chainfun <- function(x){

  my975 <- function(a) quantile(a,0.975)
  my025 <- function(a) quantile(a,0.025)
  tmpf <- function(y,func){
    switch(as.character(length(dim(y))),
           "1" = eval(call(func,y)),
           "2" = eval(call("apply",y,2,func)),
           "3" = eval(call("apply",y,c(2,3),func)))
  }
  return(list(
    avg = tmpf(x,func="mean"),
    median = tmpf(x,func="median"),
    sd = tmpf(x,func="sd"),
    u95 = tmpf(x,func="my975"),
    l95 = tmpf(x,func="my025")
  )
  )
}

chaintab <- lapply(chains,chainfun)


#################################################################################### beta ##############################################################################################################################
betagg <- bind_rows(as_tibble(chaintab$beta$median,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="median"),
                    as_tibble(chaintab$beta$l95,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="l95"),
                    as_tibble(chaintab$beta$u95,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="u95")
) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(Prmtr=factor(Prmtr,levels=c("beta0","secchi.z","elevation.z","depth.z","lakearea.z","total.dev.z","total.ag.z"),
                      labels=c("Intercept","Secchi","Elevation","Max Depth","Lake area","Developed","Agriculture")),
         zero=ifelse(l95<0 & u95>0,"0","1"))

# summary for manual control of x-axis
betagg %>%
  filter(Prmtr != "Intercept") %>%
  group_by(Species) %>%
  summarise(min=min(l95),
            max=max(u95))

p <- betagg %>%
  filter(Prmtr != "Intercept") %>%
  mutate(Species=factor(Species,levels=c("cisco","northern_pike",
                                         "walleye","black_crappie",
                                         "yellow_perch","smallmouth_bass",
                                         "largemouth_bass","bluegill"),
                        labels=c("Cisco","Northern pike",
                                 "Walleye","Black crappie",
                                 "Yellow perch","Smallmouth bass",
                                 "Largemouth bass","Bluegill"))) %>%
  ggplot() +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(aes(y=Prmtr,x=median,color=zero,shape=zero),size=3) +
  geom_errorbar(aes(y=Prmtr,xmin=l95,xmax=u95,color=zero)) +
  facet_wrap(~Species,scales = "free_x",
             nrow=2) +
  scale_color_manual(values = c("0"="#8E8E8E","1"='black')) +
  scale_shape_manual(values = c("0"=17,"1"=19)) +
  labs(x=expression(beta)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=14,angle=45,hjust=1),
        axis.text.y = element_text(size=20),
        #panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18,face="bold",color="black"),
        legend.position = "none") +
  ggh4x::facetted_pos_scales(x = list(
    Species == "Cisco" ~ scale_x_continuous(limits = c(-2.1,3),
                                            breaks=c(-1,0,2),
                                            labels=scales::number_format(accuracy = 0.1,
                                                                         decimal.mark = '.')),
    Species == "Northern pike" ~ scale_x_continuous(limits = c(-0.7,0.6),
                                                    breaks=c(-0.5,0,0.5)),
    Species == "Walleye" ~ scale_x_continuous(limits = c(-0.35,1.1),
                                              breaks=c(-0.2,0,0.5,1)),
    Species == "Black crappie" ~ scale_x_continuous(limits = c(-0.6,0.3),
                                                    breaks=c(-0.4,0,0.2)),
    Species == "Yellow perch" ~ scale_x_continuous(limits = c(-0.55,0.9),
                                                   breaks=c(-0.4,0,0.6)),
    Species == "Smallmouth bass" ~ scale_x_continuous(limits = c(-0.75,0.85),
                                                      breaks=c(-0.5,0,0.5)),
    Species == "Largemouth bass" ~ scale_x_continuous(limits = c(-0.2,0.6),
                                                      breaks=c(-0.1,0,0.5)),
    Species == "Bluegill" ~ scale_x_continuous(limits = c(-0.3,0.6),
                                               breaks=c(-0.2,0,0.5))))

mypal <- c("#fee090","#f46d43","#d73027","#a50026","#313695", "#4575b4","#74add1","#e6f598")
#mypal <- c("#fddbc7","#f4a582","#d6604d","#C41C0E","#2166ac", "#4393c3", "#92c5de", "#d1e5f0")
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- mypal[k]
  k <- k+1
}
grid::grid.draw(g)
betaplot <- arrangeGrob(g)
ggsave(filename = "figures/betaplot.jpeg",plot = betaplot,
       units = "in",width = 14,height = 10,dpi = 600)

betagg %>%
  filter(Prmtr != "Intercept") %>%
  summarise(x=sum(as.numeric(zero)),y=mean(as.numeric(zero)))

betagg %>%
  filter(Prmtr != "Intercept") %>%
  group_by(Prmtr) %>%
  summarise(x=sum(as.numeric(zero)),
            y=mean(as.numeric(zero)))





#################################################################################### theta #################################################################################
gearcodes <- c("IA_BS_2n" = "IA_DC Boat Shocker Day 2 Netters",
               "IA_BS" = "IA_DC Boat Shocker Day",
               "IA_FN" = "IA_Fyke Netting Unspecified",
               "IA_GN" = "IA_Gill Net Unspecified",
               "IL_EF_AC" = "IL_Boat electrofishing AC (Day)",
               "IL_EF_DC" = "IL_Boat electrofishing DC (Day)",
               "IN_GN_exp" = "IN_experimental gill net",
               "IN_EF" = "IN_dc nighttime electrofishing",
               "IN_GN_std" = "IN_gill net std exp",
               "IN_TN" = "IN_standard trap net",
               "MI_GN" = "MI_inland_gill_net",
               "MI_FN" = "MI_large_mesh_fyke_net",
               "MI_TN" = "MI_trap_net",
               "MI_BS" = "MI_boomshocking",
               "MN_GN" = "MN_Standard gill net sets",
               "MN_TN2f" = "MN_Standard 3/4-in mesh, double frame trap net sets",
               "MN_GNsh" = "MN_Standard gill nets, set shallow in stratified assessment",
               "SD_BS_night" = "SD_boat shocker (night)",
               "SD_FN34" = "SD_frame net (std 3/4 in)",
               "SD_GN_exp" = "SD_std exp gill net",
               "SD_GN_std" = "SD_AFS std gill net",
               "WI_FN" = "WI_fyke_net",
               "WI_BS" = "WI_boom_shocker",
               "WI_VGN" = "WI_vertical_gill_net")

gc_df <- tibble(theta=1:length(gearcodes),
                code=names(gearcodes),
                fullgear=as.character(gearcodes))

gc_df %>%
  mutate(theta=paste0("$\\theta_{",theta,"}$")) %>%
  knitr::kable(format = "latex", digits=2, 
               linesep="",escape = FALSE,booktabs=TRUE,
               caption = "Sampling gear types and their respective codes 
                          used within manuscript.")
  


thetachain <- tibble()
for(i in 1:8){
  x <- as_tibble(chains$theta[,i,]) %>%
    rownames_to_column("iter") %>%
    mutate(iter=as.numeric(iter),
           spp = dimnames(chains$theta)[[2]][i])
  thetachain <- bind_rows(thetachain,x)
}

thetagg <- list()
for(i in 1:8){
  thetagg[[sppv[i]]] <- thetachain %>%
          pivot_longer(-c(spp,iter),names_to = "gears",values_to = "theta") %>%
          filter(spp == sppv[i]) %>%
          ggplot(aes(x=iter,y=theta)) +
          geom_line() +
          ggtitle(sppv[i]) +
          facet_wrap(~gears,scales="free",nrow=6) +
          theme_bw()
          
}

thetagg[[1]]


theta_df <- bind_rows(lapply(chaintab$theta,function(x) as_tibble(x,rownames="spp")),
                      .id="metric") %>%
  pivot_longer(-c(metric,spp),names_to = "gears",values_to = "x") %>%
  mutate(x=x*dat$J) %>%
  pivot_wider(names_from="metric",values_from = "x")

theta_df %>%
  group_by(spp) %>%
  summarise(x=max(median),
            y=which(median==x),
            z=gears[y])




#################################################################################### Pi and Sigma ########################################################################################


Sigmagg <- bind_rows(as_tibble(chaintab$Sigma$median,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="median"),
                     as_tibble(chaintab$Sigma$l95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="l95"),
                     as_tibble(chaintab$Sigma$u95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="u95")) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(median=ifelse(l95<0 & u95>0,sprintf("%.1f",median),paste0(sprintf("%.1f",median),"*"))) %>%
  mutate(txt=ifelse(Sp1==Sp2,"",
                    paste0(median,
                           "\n (",
                           sprintf("%.1f",l95),
                           ", ",
                           sprintf("%.1f",u95),
                           ")")))

taugg <- bind_rows(as_tibble(chaintab$tau$median,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="median"),
                   as_tibble(chaintab$tau$l95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="l95"),
                   as_tibble(chaintab$tau$u95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="u95")) %>%
  select(-Sp2) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(txt=paste0(sprintf("%.1f",median),"\n (",sprintf("%.1f",l95),", ",sprintf("%.1f",u95),")"))

Pigg <- bind_rows(as_tibble(chaintab$Pi$median,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="median"),
                   as_tibble(chaintab$Pi$l95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="l95"),
                   as_tibble(chaintab$Pi$u95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="u95")) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(median=ifelse(l95<0 & u95>0,sprintf("%.1f",median),paste0(sprintf("%.1f",median),"*"))) %>%
  mutate(txt=ifelse(Sp1==Sp2,"",
                    paste0(median,
                           "\n (",
                           sprintf("%.1f",l95),
                           ", ",
                           sprintf("%.1f",u95),
                           ")")))

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat,diag=TRUE)]<- NA
  return(cormat)
}

uppermat <- function(x){
  y = get_upper_tri(as.matrix(x))
  return(y)
}

sigcormed <- as_tibble(chaintab$Sigma$median,rownames="Sp1") %>%
  pivot_longer(-Sp1,names_to = "Sp2",values_to = "sigma") %>%
  arrange(Sp1,Sp2) %>%
  pivot_wider(names_from = Sp2,values_from = sigma) %>%
  column_to_rownames("Sp1") %>%
  uppermat()
  

taumed <- diag(chaintab$tau$median)
taumed[lower.tri(taumed,diag=FALSE)] = NA
taumed[upper.tri(taumed,diag=FALSE)] <- NA
dimnames(taumed) = dimnames(sigcormed)

Pimed <- t(as_tibble(chaintab$Pi$median,rownames="Sp1") %>%
             pivot_longer(-Sp1,names_to = "Sp2",values_to = "sigma") %>%
             arrange(Sp1,Sp2) %>%
             pivot_wider(names_from = Sp2,values_from = sigma) %>%
             column_to_rownames("Sp1") %>%
             uppermat())

melted_sigcormed <- reshape2::melt(sigcormed, na.rm = FALSE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  left_join(select(Sigmagg,Sp1,Sp2,txt),
            by=c("Sp1","Sp2")) %>%
  arrange(Sp1) %>%
  mutate(Sp1 = factor(Sp1,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")),
         Sp2 = factor(Sp2,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")))

melted_taumed <- reshape2::melt(taumed, na.rm = TRUE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  left_join(select(taugg,Sp1,txt),
            by=c("Sp1")) %>%
  mutate(Sp1 = factor(Sp1,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")),
         Sp2 = factor(Sp2,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")))


melted_Pimed <- reshape2::melt(Pimed, na.rm = TRUE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  left_join(select(Pigg,Sp1,Sp2,txt),
            by=c("Sp1","Sp2")) %>%
  mutate(Sp1 = factor(Sp1,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")),
         Sp2 = factor(Sp2,
                      levels=c("black_crappie","bluegill","cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch"),
                      labels=c("Black crappie","Bluegill","Cisco","Largemouth bass","Northern pike","Smallmouth bass","Walleye","Yellow perch")))


jpeg("Figures/speccov.jpeg", width = 15, height = 13, units = 'in', res = 600)
ggplot() +
  geom_tile(data = melted_sigcormed, 
            aes(Sp2, Sp1, fill = value),color="black") +
  geom_tile(data=melted_Pimed,aes(Sp2,Sp1),fill="white",color='black') +
  geom_tile(data=melted_taumed,aes(Sp2,Sp1),fill="gray",color="black",linewidth=1) +
  geom_text(data = filter(melted_sigcormed,!is.na(value)), aes(Sp2, Sp1,label=txt),size=5) +
  geom_text(data = melted_taumed, aes(Sp2, Sp1,label=txt),size=5) +
  geom_text(data = melted_Pimed, aes(Sp2, Sp1,label=txt),size=5) +
  scale_fill_gradient2(low="red",
                       mid="white",
                       high="blue",
                       limits = c(-1,1),) +
  theme_void() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 18, hjust = 1),
        axis.text.y = element_text(size=18),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  coord_fixed() +
  labs(fill="Species \ndependency")
dev.off()


as_tibble(chaintab$A$median,rownames="A") %>%
  pivot_longer(-A,names_to = "spp",values_to = "median") %>%
  group_by(spp) %>%
  summarise(mx=max(median),
            mn=min(median))


#Overlapping zero
colSums(((chaintab$A$l95) < 0 & (chaintab$A$u95) > 0))
colMeans(((chaintab$A$l95) < 0 & (chaintab$A$u95) > 0))


############################################## Full parameter kable ##############################################


gearcodes <- c("IA_BS_2n" = "IA_DC Boat Shocker Day 2 Netters",
               "IA_BS" = "IA_DC Boat Shocker Day",
               "IA_FN" = "IA_Fyke Netting Unspecified",
               "IA_GN" = "IA_Gill Net Unspecified",
               "IL_EF_AC" = "IL_Boat electrofishing AC (Day)",
               "IL_EF_DC" = "IL_Boat electrofishing DC (Day)",
               "IN_GN_exp" = "IN_experimental gill net",
               "IN_EF" = "IN_dc nighttime electrofishing",
               "IN_GN_std" = "IN_gill net std exp",
               "IN_TN" = "IN_standard trap net",
               "MI_GN" = "MI_inland_gill_net",
               "MI_FN" = "MI_large_mesh_fyke_net",
               "MI_TN" = "MI_trap_net",
               "MI_BS" = "MI_boomshocking",
               "MN_GN" = "MN_Standard gill net sets",
               "MN_TN2f" = "MN_Standard 3/4-in mesh, double frame trap net sets",
               "MN_GNsh" = "MN_Standard gill nets, set shallow in stratified assessment",
               "SD_BS_night" = "SD_boat shocker (night)",
               "SD_FN34" = "SD_frame net (std 3/4 in)",
               "SD_GN_exp" = "SD_std exp gill net",
               "SD_GN_std" = "SD_AFS std gill net",
               "WI_FN" = "WI_fyke_net",
               "WI_BS" = "WI_boom_shocker",
               "WI_VGN" = "WI_vertical_gill_net")

gc_df <- tibble(theta=1:length(gearcodes),
                code=names(gearcodes),
                fullgear=as.character(gearcodes))




theta_med <-t(chaintab$theta$median)
theta_l95 <- t(chaintab$theta$l95)
theta_u95 <- t(chaintab$theta$u95)
thetamat <- as_tibble(matrix(paste0(sprintf("%.2f",theta_med*nrow(gc_df)),
                          " (",
                          sprintf("%.2f",theta_l95*nrow(gc_df)),
                          ", ",
                          sprintf("%.2f",theta_u95*nrow(gc_df)),
                          ")"),
                   nrow=nrow(theta_med),
                   dimnames=dimnames(theta_med)),
                   rownames="Parameter") %>%
  mutate(Parameter = factor(Parameter,
                            levels=gc_df$code,
                            labels=paste0("\\theta_{",gc_df$theta,"}"))) %>%
  arrange(Parameter)

beta_med <- chaintab$beta$median
beta_l95 <- chaintab$beta$l95
beta_u95 <- chaintab$beta$u95
betamat <- as_tibble(matrix(paste0(sprintf("%.2f",beta_med),
                         " (",
                         sprintf("%.2f",beta_l95),
                         ", ",
                         sprintf("%.2f",beta_u95),
                         ")"),
                  nrow=nrow(beta_med),
                  dimnames=dimnames(beta_med)),
                  rownames="Parameter") %>%
  mutate(Parameter = factor(Parameter,
                            levels=c("beta0","secchi.z","elevation.z","depth.z","lakearea.z","total.dev.z","total.ag.z"),
                            labels=c("\\beta_0","\\beta_{secchi}","\\beta_{elevation}","\\beta_{depth}","\\beta_{lake area}","\\beta_{developed}","\\beta_{agriculture}"))) %>%
  arrange(Parameter)

sum(!(beta_l95<0 & beta_u95>0))
mean(!(beta_l95<0 & beta_u95>0))

A_med <- chaintab$A$median
A_l95 <- chaintab$A$l95
A_u95 <- chaintab$A$u95
Amat <- as_tibble(matrix(paste0(
  sprintf("%.2f",A_med),
  " (",
  sprintf("%.2f",A_l95),
  ", ",
  sprintf("%.2f",A_u95),
  ")"),
  nrow=nrow(A_med),
  dimnames=dimnames(A_med))) %>%
  mutate(Parameter=paste0("A_{",1:32,"}"))

#Species summaries
summary((A_med))

#Overlapping zero
colSums(((A_l95) < 0 & (A_u95) > 0))
colMeans(((A_l95) < 0 & (A_u95) > 0))

spp_order <- c("bluegill","black_crappie", "cisco","largemouth_bass","northern_pike","smallmouth_bass","walleye","yellow_perch")

sigma_mat <- matrix(paste0(sprintf("%.2f",chaintab$Sigma$median),
                           " (",
                           sprintf("%.2f",chaintab$Sigma$l95),
                           ", ",
                           sprintf("%.2f",chaintab$Sigma$u95),
                           ")"),
                    nrow=nrow(chaintab$Sigma$median),
                    dimnames = dimnames(chaintab$Sigma$avg))[spp_order,spp_order]


sigma_mat[upper.tri(sigma_mat,diag = FALSE)] <- "-"
diag(sigma_mat) <- 1

sigma_mat <- as_tibble(sigma_mat,rownames = "Parameter") %>%
  mutate(Parameter = factor(Parameter,
                            levels=spp_order,
                            labels=c("\\Sigma_{bluegill}","\\Sigma_{black crappie}", "\\Sigma_{cisco}","\\Sigma_{largemouth bass}",
                                     "\\Sigma_{northern pike}","\\Sigma_{smallmouth bass}","\\Sigma_{walleye}","\\Sigma_{yellow perch}"))) %>%
  arrange(Parameter)


tau_med <- chaintab$tau$median
tau_l95 <- chaintab$tau$l95
tau_u95 <- chaintab$tau$u95
taumat <- matrix(paste0(sprintf("%.2f",tau_med),
                        " (",
                        sprintf("%.2f",tau_l95),
                        ", ",
                        sprintf("%.2f",tau_u95),
                        ")"),
                 nrow=1)

dimnames(taumat)[[1]] <- "T"
dimnames(taumat)[[2]] <- names(tau_med)

taumat <- as_tibble(taumat,rownames="Parameter")


fullkab <- bind_rows(betamat,
                     thetamat,
                     Amat,
                     taumat,
                     sigma_mat) %>%
  select(Parameter,bluegill,black_crappie,cisco,largemouth_bass,northern_pike,smallmouth_bass,walleye,yellow_perch)

fullkab %>%
  mutate(Parameter=paste0("$",Parameter,"$")) %>%
  knitr::kable(format = "latex", digits=2, linesep="",escape = FALSE,
               booktabs=TRUE,longtable=TRUE,
               caption = "Median (95\\% credible intervals) parameter estimates for each species of Minnesota lake fish modeled under the jsPGA framework.") %>%
  kableExtra::kable_styling(latex_options=c("repeat_header")) %>%
  kableExtra::landscape()

################################################################################ Map ################################################################################
library(sf)
fishdat <- readRDS("data/MWfishz.rds") %>%
  select(-starts_with("comp"))

usa <-st_as_sf(maps::map("state",fill=TRUE,plot=FALSE))


fishsf <- fishdat %>%
  select(COMMON_NAME,TOTAL_CATCH,lon,lat) %>%
  filter(TOTAL_CATCH > 0) %>%
  group_by(COMMON_NAME,lon,lat) %>%
  summarise(avg=mean(TOTAL_CATCH)) %>%
  ungroup() %>%
  st_as_sf(coords=c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(avg) %>%
  mutate(COMMON_NAME = factor(COMMON_NAME,levels=c("cisco","northern_pike",
                                                   "walleye","black_crappie",
                                                   "yellow_perch","smallmouth_bass",
                                                   "largemouth_bass","bluegill"),
                              labels=c("Cisco","Northern pike",
                                       "Walleye","Black crappie",
                                       "Yellow perch","Smallmouth bass",
                                       "Largemouth bass","Bluegill")))

mypal <- c("#a50026","#d73027","#f46d43","#fee090","#e6f598","#74add1", "#4575b4","#313695")
mypal <- mypal[8:1]
#mypal <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0","#fddbc7","#f4a582","#d6604d","#C41C0E")

sppv <- c("Cisco","Northern pike",
          "Walleye","Black crappie",
          "Yellow perch","Smallmouth bass",
          "Largemouth bass","Bluegill")


fishgg <- list()
for(k in 1:length(sppv)){

  mybreaks = as.integer(quantile(filter(fishsf,COMMON_NAME==sppv[k])$avg,c(0,0.5,1)))
  
  p <- {ggplot() +
      geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
      geom_sf(data=filter(fishsf,COMMON_NAME==sppv[k]),aes(color=avg)) +
      coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
      facet_wrap(~COMMON_NAME,nrow=2) +
      scale_color_viridis_c(na.value = "#5E5D61",trans="log",breaks=mybreaks) +
      scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
      theme_bw() +
      labs(color="Mean catch") +
      theme(text=element_text(size=16),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            strip.text = element_text(face = "bold"),
            plot.margin = unit(c(0,0,0,0),"cm"),
            legend.position = "inside",
            legend.position.inside = c(0.23,0.15),
            legend.background = element_blank(),
            legend.direction = "horizontal",
            legend.title.position = "top")}

fishgg[[k]] <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', fishgg[[k]]$layout$name))
for (i in strip_both) {
  j <- which(grepl('rect', fishgg[[k]]$grobs[[i]]$grobs[[1]]$childrenOrder))
  fishgg[[k]]$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- mypal[k]
}
}
gridExtra::grid.arrange(grobs=fishgg,nrow=2)
gg1 <- gridExtra::arrangeGrob(grobs=fishgg,nrow=2)
ggsave("Figures/fishcatchmap.jpeg", plot=gg1,width = 15, height = 9, units = 'in', dpi = 600)

fishdat %>%
  filter(TOTAL_CATCH > 0) %>%
  group_by(COMMON_NAME) %>%
  summarise(n=n_distinct(site_id),
            max=max(TOTAL_CATCH),
            avg=mean(TOTAL_CATCH),
            med=median(TOTAL_CATCH))



lakedat <- readRDS("Data/MWlakesall.rds")

covv <- c("max_depth","elevation","lakearea","total_ag","total_dev","secchi")
covlabs <- c("Max depth (m)","Elevation (m)","Lake area (ha)","Agriculture (%)","Developed (%)","Secchi depth (m)")

lakesf <- lakedat %>%
  mutate(total_ag=total_ag/100,
         total_dev=total_dev/100) %>%
  select(site_id,all_of(covv),lon,lat) %>%
  distinct() %>%
  pivot_longer(-c(site_id,lon,lat),names_to = "covs",values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(covs=factor(covs,levels=covv,labels=covlabs)) %>%
  st_as_sf(coords=c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(value)


mylogit <- scales::trans_new("mylogit",
                             transform = function(x) ifelse(x>=1,log(0.999/(1-0.999)),
                                                            ifelse(x<=0,log(0.0001/(1-0.0001)),
                                                                   log(x/(1-x)))),
                             inverse = function(x) exp(x)/(exp(x) + 1)
                             )

lakegg <- list()
for(k in 1:length(covlabs)){
  
  mybreaks = quantile(filter(lakesf,covs==covlabs[k])$value,c(0,0.5,1))
  
  if(covlabs[k] %in% c("Lake area (ha)","Max depth (m)")){
    myscale <- scale_color_viridis_c(na.value = "#5E5D61",trans="log",labels=round,breaks=mybreaks)
  }else if(covlabs[k] %in% c("Agriculture (%)","Developed (%)")){
    myscale <- scale_color_viridis_c(na.value = "#5E5D61",trans=mylogit,breaks=mybreaks,labels=function(x) sprintf("%.0f",x*100))
  }else{
    myscale <- scale_color_viridis_c(na.value = "#5E5D61",labels=round,breaks=mybreaks)
  }
  
  lakegg[[k]] <- {ggplot() +
      geom_sf(data=filter(usa,ID %in% c("minnesota","michigan","wisconsin","south dakota","illinois","indiana","iowa"))) + 
      geom_sf(data=filter(lakesf,covs==covlabs[k]),aes(color=value)) +
      coord_sf(xlim=c(-97.47815,-83.0606),y=c(37.15983,49.07641)) +
      facet_wrap(~covs,nrow=2) +
      myscale +
      scale_x_continuous(breaks=c(-96,-92,-88,-84)) +
      theme_bw() +
      labs(color="Value") +
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            strip.text = element_text(face = "bold",size=14),
            plot.margin = unit(c(0,0,0,0),"cm"),
            legend.position = "inside",
            legend.title = element_text(size=14),
            legend.position.inside = c(0.23,0.15),
            legend.background = element_blank(),
            legend.direction = "horizontal",
            legend.title.position = "top",
            legend.text = element_text(size=10))}
}
gridExtra::grid.arrange(grobs=lakegg,nrow=2)
gg1 <- gridExtra::arrangeGrob(grobs=lakegg,nrow=2)
ggsave("Figures/covsmap.jpeg", plot=gg1,width = 10, height = 8, units = 'in', dpi = 600)




st_drop_geometry(lakesf) %>%
  group_by(covs) %>%
  summarise(n=n_distinct(site_id),
            min=min(value),
            avg=mean(value),
            med=median(value),
            max=max(value))














