#climate change, land use change, and traits -- model selection

# set up -----------------------------------------------------------------------
# load some libraries
library(ggplot2); library(cowplot); library(bbsBayes); library(tidybayes)
library(ggpattern); library(coda); library(spdep); library(brinla); library(scico)
library(inlatools); library(sf); library(INLA); library(tidyr); library(dplyr)
library(stringr); library(ggpubr); library(glmmTMB);library(egg); library(inlabru);
library(ape); library(phytools); library(MCMCglmm)


# set some options
inla.setOption(inla.mode="experimental")
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

#simple scaling function
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


# climate by voltinism interactions (prepare data) ---------------------------------------
#link urbanization, agriculture, and climate change to trends
#import data
df.lulc_means <- read.csv("data/climate and land use means and trends_wide.csv")
df.strata<- read.csv("data/trends/trends by strata_1996-2018.csv")
#load covariate data
env.vars2 <- c("tmp","pre","agriculture","natural","urban")
map1 <- load_map(stratify_by="bbs_usgs") %>% rename_all(tolower) %>% 
  dplyr::select(2, 3) %>% rename(strat_name=st_12, strat_A=area_1)
df.lulc <- read.csv("data/land use/lulc_by year and strata.csv") %>%
  mutate(land_type = car::recode(lulc_state, "0='empty'; 22='cropland';33='pasture'; 40 = 'forest'; 41 = 'forest'; 55 = 'grass_shrubland';66 = 'other';77 = 'water';11 = 'urban';43 = 'forest';45 = 'forest';44 = 'forest'"),land_type2 = car::recode(land_type, "'empty' = 'empty';'cropland' = 'agriculture';'pasture' = 'agriculture';'forest' = 'natural';'grass_shrubland' = 'natural';'urban' = 'urban';'other' = 'other';'water' = 'water'"), freq = as.numeric(freq),year = as.numeric(year)) %>% 
  left_join(map1[,1:2],by="strat_name") %>% 
  st_drop_geometry() %>% 
  dplyr::select(-c(lulc_state,X)) %>% 
  group_by(strat_name,year,land_type2,strat_A) %>%
  summarize(freq=sum(freq))%>%
  ungroup() %>%
  mutate(freq= freq/strat_A)%>%
  spread(land_type2,freq)
df.clim<- read.csv("data/climate/climate_by strata and year.csv")
df.env <- df.lulc %>% left_join(df.clim,by=c("strat_name","year")) %>% 
  select(strat_name,year,env.vars2) %>%
  group_by(strat_name) %>% 
  mutate_at(c(env.vars2),scale_this) %>% 
  ungroup()%>%
  left_join(df.lulc_means[,c("strat_name","tmp_mean","pre_mean","urb_mean","ag_mean")],by="strat_name")

#create single dataframe with all species counts*strat*year
sp.vect <- str_split(list.files("data/sampled from model_strata"),"_",simplify = T)[,1]
sampled.abun.list <- list()
for(i in 1:length(sp.vect)){
  d0 <- read.csv(paste0("data/sampled from model_strata/",sp.vect[i],"_samps3.csv"))
  if("strat_nAz" %in% colnames(d0)){
    sampled.abun.list[[i]] <- d0 %>% 
      mutate(species = sp.vect[i],count = strat_nAz)%>%
      group_by(strat_name,year,species)%>%
      summarize(y = quantile(count,0.5)[1],ucl = quantile(count,0.975)[1],lcl =  quantile(count,0.025)[1], log_y=log(y), log_prec_y=1/(((log(ucl)-log(lcl))/(1.96*2))^2)) %>% 
      select(species,strat_name,year,log_y,log_prec_y)
  } else( 
    sampled.abun.list[[i]] <- d0%>%
      mutate(species = sp.vect[i],count = aggregate_N)%>%
      group_by(strat_name,year,species)%>%
      summarize(y = quantile(count,0.5)[1],ucl = quantile(count,0.975)[1],lcl =  quantile(count,0.025)[1], log_y=log(y), log_prec_y=1/(((log(ucl)-log(lcl))/(1.96*2))^2)) %>% 
      select(species,strat_name,year,log_y,log_prec_y)
  )
}

# create map for inla model -----------------------------------------------
sites.list <- list()
for(i in 1:length(sp.vect)){
  sites.list[[i]] <- read.csv(paste0("data/site strat key/",sp.vect[i],"_site_strat_key.csv")) %>% 
    filter(n_sites >= 3 & strat_z > 0 & !(is.na(site))) %>% 
    mutate(lon = str_split(site,"_",simplify = T)[,2],lat = str_split(site,"_",simplify = T)[,3]) %>% 
    select(site,lon,lat)
}
sites1 <- do.call("rbind",sites.list) %>% 
  as.data.frame() %>% 
  distinct(site,lat,lon)%>%
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
  st_transform(crs=st_crs(map1))

# pair each site with its respective bbs stratum
within1 <- sapply(st_within(sites1, map1), function(z) if (length(z)==0) 
  NA_integer_ else z[1])
near1 <- st_nearest_feature(sites1, map1) 
poly_per_site <- ifelse(is.na(within1), near1, within1)
site_strat_key <- sites1 %>% 
  mutate(strat_name=unique(map1$strat_name)[poly_per_site]) %>% 
  st_drop_geometry() %>% 
  right_join(map1 %>% st_drop_geometry()) %>% 
  select(site, strat_name, strat_A)

area1 <- st_convex_hull(st_combine(sites1)) %>% st_as_sf()
strat_in_area <- sapply(st_intersects(map1, area1), 
                        function(z) if (length(z)==0) 
                          NA_integer_ else z[1])
map2 <- map1[which(strat_in_area==1), ] %>% 
  mutate(strat_idx1=as.integer(factor(strat_name))) %>% 
  select(strat_name, strat_idx1, strat_A)

# make the neighborhood matrix and graph for INLA using map2
# ggplot() + geom_sf(data=map2) + 
#   scale_color_scico(palette = "vik", direction=-1,
#                     begin=0.4, end=0.1) +
#   geom_sf(data=area1, col="black", lty=2, fill=NA) + theme_bw()
nb1 <- poly2nb(as(map2, "Spatial"))
summary(nb1)
nb2INLA("map_allspecies.adj", nb1)
nb2INLA("data/neighborhood graphs/map_allspecies.adj", nb1)
g1 <- inla.read.graph(filename = "data/neighborhood graphs/map_allspecies.adj") #load graph

# append data for abundance, covariate  (z-scaled), and strat_id  ---------
abund.df <- do.call("rbind",sampled.abun.list) %>% 
  as.data.frame() %>% 
  left_join(df.env[,c("strat_name","year",env.vars2,"tmp_mean","pre_mean","urb_mean","ag_mean")],by=c("strat_name","year")) %>% 
  left_join(map2 %>% st_drop_geometry() %>% select(-strat_A),by=c("strat_name")) %>% 
  left_join(df.global[,c("Species","Family")],by=c("species"="Species"))
write.csv(abund.df,"data/all species counts_by strata and year_with covariates.csv")


# run trait by climate models  --------------------------------------------
#load all models
load("model selection list.RData")
#define priors
bym_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                  phi = list(prior = "pc", param = c(0.5, 0.5)))
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

#load inla graph and abundance data (with covariates)
g1 <- inla.read.graph(filename = "neighborhood graphs/map_allspecies.adj")
d1 <- read.csv("data/all species counts_by strata and year_with covariates.csv") %>% 
  left_join(sp.order,by = c("species" = "Species"))
d1$strat_idx2 = d1$strat_idx3 = d1$strat_idx4 = d1$strat_idx5 = d1$strat_idx1
d1$sp_idx2 = d1$sp_idx3 = d1$sp_idx4 = d1$sp_idx5 = d1$sp_idx1 
# d1$family_idx1 = as.integer(factor(d1$Family))
d1<-d1 %>% left_join(df.strata[,c("species","strat_name","brood")],by = c("species","strat_name"))

waic.tab <- data.frame(model = model_names,model_type,waic = 1:57)
for(i in 1:57){
  res1 <- inla(mod.list[[i]], data=d1,scale = log_prec_y,
               family = "gaussian",
               control.predictor=list(compute=T, link=1),
               control.compute=list(config=T, waic = T, cpo=T,
                                    control.gcpo=list(enable=T)),
               verbose=T, inla.mode="experimental",
               control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
               control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))
  waic.tab$waic[i]<-ifelse(is.null(res1$waic),NA,res1$waic$waic[[1]])
}
waic.tab$dwaic <- waic.tab$waic-min(waic.tab$waic)

write.csv(waic.tab %>% mutate_at(c("waic","dwaic"), round,1) %>% arrange(dwaic),"tables/trends vs covariates and traits_aic table_long.csv")


