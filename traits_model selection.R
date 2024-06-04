#model selection for butterfly traits

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

load("traits and phylo.RData")

# create map for inla model (traits complete) -----------------------------------------------
sp.vect <- str_split(list.files("data/sampled from model_strata"),"_",simplify = T)[,1]
map1 <- load_map(stratify_by="bbs_usgs") %>% rename_all(tolower) %>% 
  dplyr::select(2, 3) %>% rename(strat_name=st_12, strat_A=area_1)
sites.list <- list()
for(i in 1:length(sp.vect)){
  sites.list[[i]] <- read.csv(paste0("data/site strat key/",sp.vect[i],"_site_strat_key.csv")) %>% 
    filter(n_sites >= 3 & strat_z > 0 & !(is.na(site))) %>% 
    mutate(lon = str_split(site,"_",simplify = T)[,2],lat = str_split(site,"_",simplify = T)[,3]) %>% 
    dplyr::select(site,lon,lat)
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
  dplyr::select(site, strat_name, strat_A)

area1 <- st_convex_hull(st_combine(sites1)) %>% st_as_sf()
strat_in_area <- sapply(st_intersects(map1, area1), 
                        function(z) if (length(z)==0) 
                          NA_integer_ else z[1])
map2 <- map1[which(strat_in_area==1), ] %>% 
  mutate(strat_idx1=as.integer(factor(strat_name))) %>% 
  dplyr::select(strat_name, strat_idx1, strat_A)

# make the neighborhood matrix and graph for INLA using map2
# ggplot() + geom_sf(data=map2) + 
#   scale_color_scico(palette = "vik", direction=-1,
#                     begin=0.4, end=0.1) +
#   geom_sf(data=area1, col="black", lty=2, fill=NA) + theme_bw()
nb1 <- poly2nb(as(map2, "Spatial"))
summary(nb1)
nb2INLA("data/neighborhood graphs/traits_allspecies.adj", nb1)
g1 <- inla.read.graph(filename = "data/neighborhood graphs/traits_allspecies.adj") #load graph

# link trends to traits (strata level -- complete trait list) ------------------------------------
#link trends to traits while accounting for strata level variation
#load strata trends, add all of the traits above
df.strata <- read.csv("data/summary trends/trends by strata_1996-2018.csv") %>% 
  mutate(trend_prec = (1/((ci_range/3.92)^2))) %>% 
  left_join(df.global[,-c(1:6,8)],by=c("species" = "Species"))%>%
  left_join(df.traits4[,c("Species","strat_name","brood")],by=c("species" = "Species","strat_name"))  %>% 
  dplyr::select(mean_trend,disturbance,edge,canopy,moisture,diet_breadth,brood,hairs,clutch,diapause,adult_aposematism,larvae_aposematism,range_size,flight_duration,size,species,strat_name,trend_prec) %>%
  left_join(sp.order.sub,by = c("species" = "Species")) %>% 
  left_join(map2 %>% st_drop_geometry() %>% dplyr::select(strat_name,strat_idx1),by = "strat_name") %>%
  na.omit()

#set priors
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

#define models
traits.res <- list()
traits.models <- c("Disturbance","Edge","Moisture","Canopy","Diet Breadth","Diapause","Oviposition","Voltinism","Aposematism (adult)","Aposematism (larvae)","Hairs","Size","Flight Duration","Range Size","Intercept")
traits.waic <- data.frame(Models = traits.models,waic=1:length(traits.models))
g1 <- inla.read.graph(filename = "map_allspecies.adj")

set.seed(032393)
#disturbance
sum(!is.na(df.global$disturbance)) #225
form1 <-  mean_trend ~ disturbance +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[1]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Disturbance", trait.value = str_split(rownames(res1$summary.fixed),"disturbance",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.res[[1]]$trait.value[1] = "Associated (strong)"
traits.waic$waic[1] <- res1$waic$waic

#edge
sum(!is.na(df.global$edge)) #245
form1 <-  mean_trend ~ edge +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[2]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Edge", trait.value = str_split(rownames(res1$summary.fixed),"edge",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.res[[2]]$trait.value[1] = "Associated (strong)"
traits.waic$waic[2] <- res1$waic$waic

#moisture
sum(!is.na(df.global$moisture)) #281
form1 <-  mean_trend ~ moisture +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[3]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Moisture", trait.value = str_split(rownames(res1$summary.fixed),"moisture",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.res[[3]]$trait.value[1] = "Mesic (strong)"
traits.waic$waic[3] <- res1$waic$waic

#canopy
sum(!is.na(df.global$canopy)) #305
form1 <-  mean_trend ~ canopy +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[4]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Canopy", trait.value = str_split(rownames(res1$summary.fixed),"canopy",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.res[[4]]$trait.value[1] = "Closed (weak)"
traits.waic$waic[4] <- res1$waic$waic

#diet_breadth
sum(!is.na(df.global$diet_breadth)) #286
form1 <-  mean_trend ~ diet_breadth +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[5]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Diet Breadth", trait.value = rownames(res1$summary.fixed)[2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[5] <- res1$waic$waic

#diapause
sum(!is.na(df.global$diapause)) #162
form1 <-  mean_trend ~ diapause +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[6]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Diapause", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'diapauseA'"),trait.value = str_split(trait.value,"diapause",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[6] <- res1$waic$waic

#clutch
sum(!is.na(df.global$clutch)) #156
form1 <-  mean_trend ~ clutch +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[7]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Oviposition", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Cluster'; 'clutchS' = 'Single'; 'clutchSC' = 'Both'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[7] <- res1$waic$waic

#voltinism
sum(!is.na(df.global$voltinism)) #307
df.strata$brood <- as.factor(df.strata$brood)
form1 <-  mean_trend ~ brood +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[8]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Voltinism", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = '1'; 'brood2' = '2'; 'brood3' = '3'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[8] <- res1$waic$waic

#adult_aposematism
sum(!is.na(df.global$adult_aposematism)) #313
form1 <-  mean_trend ~ adult_aposematism +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[9]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Aposematism (adult)", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Aposematic (adult)'; 'adult_aposematismcryptic' = 'Cryptic (adult)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[9] <- res1$waic$waic

#larvae_aposematism
sum(!is.na(df.global$larvae_aposematism)) #248
form1 <-  mean_trend ~ larvae_aposematism +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[10]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Aposematism (larva)", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Aposematic (larva)'; 'larvae_aposematismcryptic' = 'Cryptic (larva)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[10] <- res1$waic$waic

#hairs
sum(!is.na(df.global$hairs)) #246
form1 <-  mean_trend ~ hairs +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[11]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Hairs", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'No'; 'hairsyes' = 'Yes'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[11] <- res1$waic$waic

#size
sum(!is.na(df.global$size)) #313
form1 <-  mean_trend ~ size +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[12]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Size", trait.value = rownames(res1$summary.fixed)[2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[12] <- res1$waic$waic

#flight duration
sum(!is.na(df.global$flight_duration)) #307
form1 <-  mean_trend ~ flight_duration +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[13]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Flight Duration", trait.value = rownames(res1$summary.fixed)[2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[13] <- res1$waic$waic

#range size
sum(!is.na(df.global$range_size)) #282
form1 <-  mean_trend ~ range_size +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res[[14]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Range Size", trait.value = rownames(res1$summary.fixed)[2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic$waic[14] <- res1$waic$waic

#intercept
form1 <-  mean_trend ~ 1 +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.sub)

res1 <- inla(form1, data=df.strata, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))
traits.waic$waic[15] <- res1$waic$waic
traits.waic$dwaic <- traits.waic$waic-min(traits.waic$waic)
write.csv(traits.waic %>% mutate_at(c("waic","dwaic"), round,1) %>% arrange(waic),"tables/butterfly trends vs traits_strata_aic table_with phylogeny_all traits.csv")

# sample size / distribution for each trait -------------------------------
gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "disturbance",stat = "count",fill = "black") +
  labs(y = "Count",x = "Habitat Preference (Disturbance)")
ggsave("plots/trait histograms/disturbance.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "edge",stat = "count",fill = "black") +
  labs(y = "Count",x = "Habitat Preference (Edge)")
ggsave("plots/trait histograms/edge.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "canopy",stat = "count",fill = "black") +
  labs(y = "Count",x = "Habitat Preference (Canopy)")
ggsave("plots/trait histograms/canopy.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "moisture",stat = "count",fill = "black") +
  labs(y = "Count",x = "Habitat Preference (Moisture)")
ggsave("plots/trait histograms/moisture.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "diet_breadth",fill = "black") +
  labs(y = "Count",x = "Diet Breadth")
ggsave("plots/trait histograms/dietbreadth.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "clutch",stat = "count",fill = "black") +
  labs(y = "Count",x = "Oviposition Behavior")
ggsave("plots/trait histograms/oviposition_behavior.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "diapause",stat = "count",fill = "black") +
  labs(y = "Count",x = "Diapause Stage")
ggsave("plots/trait histograms/diapause.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "hairs",stat = "count",fill = "black") +
  labs(y = "Count",x = "Larval Hairs")
ggsave("plots/trait histograms/hairs.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "adult_aposematism",stat = "count",fill = "black") +
  labs(y = "Count",x = "Adult Aposematism")
ggsave("plots/trait histograms/adult_aposematism.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "larvae_aposematism",stat = "count",fill = "black") +
  labs(y = "Count",x = "Larval Aposematism")
ggsave("plots/trait histograms/larval_aposematism.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "range_size",fill = "black",bins = 10) +
  labs(y = "Count",x = "Range Size (z-transformed)")
ggsave("plots/trait histograms/rangesize.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "flight_duration",fill = "black",bins = 10) +
  labs(y = "Count",x = "Flight Duration (z-transformed)")
ggsave("plots/trait histograms/flight_duration.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.global %>% filter(Species %in% unique(df.strata$species)),x = "size",fill = "black",bins = 10) +
  labs(y = "Count",x = "Size (z-transformed)")
ggsave("plots/trait histograms/size.jpg",bg = "white",height = 6,width = 10)

gghistogram(df.strata,x = "brood",fill = "black",stat = "count") +
  labs(y = "Count",x = "Voltinism")
ggsave("plots/trait histograms/voltnism.jpg",bg = "white",height = 6,width = 10)

# get plot of each trait ranked by WAIC (nas omitted) ----------------------
traits.res.df <- do.call(rbind,traits.res) %>% 
  left_join(traits.waic,by=c("trait.cat"="Models")) %>% 
  mutate(trait.cat = car::recode(trait.cat,"'Aposematism (larva)' = 'Aposematism';'Aposematism (adult)' = 'Aposematism'"))

trait.levels<-traits.waic %>% mutate(Models = car::recode(Models,"'Aposematism (larva)' = 'Aposematism';'Aposematism (adult)' = 'Aposematism'")) %>% group_by(Models) %>% summarize(dwaic = mean(dwaic)) %>% arrange(dwaic)

traits.res.df$trait.value2 <- factor(traits.res.df$trait.value, levels = c("1", "2", "3","Single","Both","Cluster","Xeric (strong)","Xeric (weak)","Associated (strong)","Associated (weak)","Variable","Mesic (weak)","Mesic (strong)","Yes","No","Avoidant (weak)","Avoidant (strong)","E","EL","ELP","EP","ELPA","L","LP","LA","LPA","P","A","Open","Open (weak)","Mixed","Closed (weak)","Closed","Cryptic (larvae)","Cryptic (adult)","Aposematic (larva)","Aposematic (adult)",""))

#plot traits results
ggplot(traits.res.df %>% filter(!(trait.cat %in% "Family")),aes(y=mean,x=trait.value2))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,position=position_dodge(.9))+
  theme_classic()+
  geom_hline(yintercept = 0,color="red")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~factor(trait.cat,levels=trait.levels$Models),scales = "free")+
  labs(x="",y = "Mean Trend (%/yr)")+
  coord_flip()
ggsave("plots/trends vs traits_strata_nas removed.jpg",bg="white",height = 7,width = 10)


# link trend to traits--show trait effects without subsetting -----------------------
#link trends to traits while accounting for strata level variation
#load strata trends, add all of the traits above
df.strata2 <- read.csv("data/summary trends/trends by strata_1996-2018.csv") %>% 
  mutate(trend_prec = (1/((ci_range/3.92)^2))) %>% 
  left_join(df.global[,-c(1:6,8)],by=c("species" = "Species"))%>%
  left_join(df.traits4[,c("Species","strat_name","brood")],by=c("species" = "Species","strat_name"))  %>% 
  dplyr::select(mean_trend,disturbance,edge,canopy,moisture,diet_breadth,brood,hairs,clutch,diapause,adult_aposematism,larvae_aposematism,range_size,flight_duration,size,species,strat_name,trend_prec) %>%
  left_join(sp.order,by = c("species" = "Species")) %>% 
  left_join(map2 %>% st_drop_geometry() %>% dplyr::select(strat_name,strat_idx1),by = "strat_name")

#set priors
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

#define models
traits.res2 <- list()
traits.models <- c("Disturbance","Edge","Moisture","Canopy","Diet Breadth","Diapause","Oviposition","Voltinism","Aposematism (adult)","Aposematism (larvae)","Hairs","Size","Flight Duration","Range Size","Intercept")
traits.waic2 <- data.frame(Models = traits.models,waic=1:length(traits.models))
g1 <- inla.read.graph(filename = "map_allspecies.adj")
set.seed(032393)

#disturbance
sum(!is.na(df.global$disturbance)) #225
form1 <-  mean_trend ~ disturbance +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[1]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Disturbance", trait.value = str_split(rownames(res1$summary.fixed),"disturbance",simplify = T)[,2],trait.value = car::recode(trait.value,"'' = 'Associated (strong)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[1] <- res1$waic$waic

#edge
sum(!is.na(df.global$edge)) #245
form1 <-  mean_trend ~ edge +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[2]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Edge", trait.value = str_split(rownames(res1$summary.fixed),"edge",simplify = T)[,2],trait.value = car::recode(trait.value,"'' = 'Associated (strong)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[2] <- res1$waic$waic

#moisture
sum(!is.na(df.global$moisture)) #281
form1 <-  mean_trend ~ moisture +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[3]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Moisture", trait.value = str_split(rownames(res1$summary.fixed),"moisture",simplify = T)[,2],trait.value = car::recode(trait.value,"'' = 'Mesic (strong)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[3] <- res1$waic$waic

#canopy
sum(!is.na(df.global$canopy)) #305
form1 <-  mean_trend ~ canopy +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[4]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Canopy", trait.value = str_split(rownames(res1$summary.fixed),"canopy",simplify = T)[,2],trait.value = car::recode(trait.value,"'' = 'Closed'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[4] <- res1$waic$waic

#diet_breadth
sum(!is.na(df.global$diet_breadth)) #286
form1 <-  mean_trend ~ diet_breadth +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[5]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Diet Breadth", trait.value = rownames(res1$summary.fixed)[2],trait.value = car::recode(trait.value,"'diet_breadth' = ''")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[5] <- res1$waic$waic

#diapause
sum(!is.na(df.global$diapause)) #162
form1 <-  mean_trend ~ diapause +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[6]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Diapause", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'diapauseA'"),trait.value = str_split(trait.value,"diapause",simplify = T)[,2]) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[6] <- res1$waic$waic

#clutch
sum(!is.na(df.global$clutch)) #156
form1 <-  mean_trend ~ clutch +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[7]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Oviposition", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Cluster'; 'clutchS' = 'Single'; 'clutchSC' = 'Both'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[7] <- res1$waic$waic

#voltinism
sum(!is.na(df.global$voltinism)) #307
df.strata2$brood <- as.factor(df.strata2$brood)
df.strata2$brood <- as.factor(df.strata2$brood)

form1 <-  mean_trend ~ brood +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[8]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Voltinism", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = '1'; 'brood2' = '2'; 'brood3' = '3'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[8] <- res1$waic$waic

#adult_aposematism
sum(!is.na(df.global$adult_aposematism)) #313
form1 <-  mean_trend ~ adult_aposematism +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[9]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Aposematism (adult)", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Aposematic (adult)'; 'adult_aposematismcryptic' = 'Cryptic (adult)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[9] <- res1$waic$waic

#larvae_aposematism
sum(!is.na(df.global$larvae_aposematism)) #248
form1 <-  mean_trend ~ larvae_aposematism +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[10]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Aposematism (larva)", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'Aposematic (larva)'; 'larvae_aposematismcryptic' = 'Cryptic (larva)'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[10] <- res1$waic$waic

#hairs
sum(!is.na(df.global$hairs)) #246
form1 <-  mean_trend ~ hairs +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[11]] <- res1$summary.fixed[,c(1,3,5)] %>% 
  mutate(trait.cat = "Hairs", trait.value = rownames(res1$summary.fixed),trait.value = car::recode(trait.value,"'(Intercept)' = 'No'; 'hairsyes' = 'Yes'")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[11] <- res1$waic$waic

#size
sum(!is.na(df.global$size)) #313
form1 <-  mean_trend ~ size +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[12]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Size", trait.value = rownames(res1$summary.fixed)[2],trait.value = car::recode(trait.value,"'size' = ''")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[12] <- res1$waic$waic

#flight duration
sum(!is.na(df.global$flight_duration)) #307
form1 <-  mean_trend ~ flight_duration +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[13]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Flight Duration", trait.value = rownames(res1$summary.fixed)[2],trait.value = car::recode(trait.value,"'flight_duration' = ''")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[13] <- res1$waic$waic

#range size
sum(!is.na(df.global$range_size)) #282
form1 <-  mean_trend ~ range_size +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.res2[[14]] <- res1$summary.fixed[2,c(1,3,5)] %>% 
  mutate(trait.cat = "Range Size", trait.value = rownames(res1$summary.fixed)[2],trait.value = car::recode(trait.value,"'range_size' = ''")) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`)
traits.waic2$waic[14] <- res1$waic$waic

#intercept
form1 <-  mean_trend ~ 1 +   
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

res1 <- inla(form1, data=df.strata2, scale = trend_prec,
             family = "gaussian", 
             control.predictor=list(compute=T, link=1),
             control.compute=list(config=T, waic=T, cpo=T,
                                  control.gcpo=list(enable=T)),
             verbose=T, inla.mode="experimental",
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
             control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))

traits.waic2$waic[15] <- res1$waic$waic
traits.waic2$dwaic <- traits.waic2$waic-min(traits.waic2$waic)
write.csv(traits.waic2 %>% mutate_at(c("waic","dwaic"), round,1) %>% arrange(waic),"tables/butterfly trends vs traits_strata_aic table_with phylogeny_all traits_no subsetting.csv")

# traits--tables and plots ------------------------------------------------
traits.res2.df <- do.call(rbind,traits.res2) %>% 
  left_join(traits.waic,by=c("trait.cat"="Models")) %>% 
  mutate(trait.cat = car::recode(trait.cat,"'Aposematism (larva)' = 'Aposematism';'Aposematism (adult)' = 'Aposematism'"))

trait.levels<-traits.waic %>% mutate(Models = car::recode(Models,"'Aposematism (larva)' = 'Aposematism';'Aposematism (adult)' = 'Aposematism'")) %>% group_by(Models) %>% summarize(dwaic = mean(dwaic)) %>% arrange(dwaic)

traits.res2.df$trait.value2 <- factor(traits.res2.df$trait.value, levels = c("1", "2", "3","Single","Both","Cluster","Xeric (strong)","Xeric (weak)","Associated (strong)","Associated (weak)","Variable","Mesic (weak)","Mesic (strong)","Yes","No","Avoidant (weak)","Avoidant (strong)","E","EL","ELP","EP","ELPA","L","LP","LA","LPA","P","A","Open","Open (weak)","Mixed","Closed (weak)","Closed","Cryptic (larvae)","Cryptic (adult)","Aposematic (larva)","Aposematic (adult)",""))

#plot traits results
ggplot(traits.res2.df %>% filter(!(trait.cat %in% "Family")),aes(y=mean,x=trait.value2))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,position=position_dodge(.9))+
  theme_classic()+
  geom_hline(yintercept = 0,color="red")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~factor(trait.cat,levels=trait.levels$Models),scales = "free")+
  labs(x="",y = "Mean Trend (%/yr)")+
  coord_flip()
ggsave("plots/trends vs traits_strata_no subsetting.jpg",bg="white",height = 7,width = 10)

#voltinism only
p.volt <- ggplot(traits.res2.df %>% filter(trait.cat %in% "Voltinism"),aes(y=mean,x=trait.value2))+
  geom_hline(yintercept = 0,color="red",linewidth=1.25)+
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,position=position_dodge(.9))+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12,color="black"))+
  # facet_wrap(~factor(trait.cat,levels=trait.levels$Models),scales = "free")+
  labs(x="Number of Broods",y = "Mean Trend (%/yr)")+
  coord_flip()
p.volt
ggsave("plots/voltinism effects on strata level trends.jpg",bg="white",height = 4,width = 5)

#compare parameter estimates between the full data and the subet (NAs removed / complete trait data available)
traits.comb = traits.res.df %>% left_join(traits.res2.df,by = c("trait.cat","trait.value"))
ggplot(traits.comb,aes(y = mean.y, x = mean.x))+
  geom_point(alpha = .5)+
  geom_errorbar(aes(ymin = lcl.y,ymax = ucl.y),alpha = .5) + 
  geom_errorbar(aes(xmin = lcl.x,xmax = ucl.x),alpha = .5)+
  # geom_smooth(method = "lm")+
  geom_abline(slope = 1, intercept = 0,color = "red")+
  theme_classic()+
  labs(y = "Trait Estimate (Subsetted Data)", x = "Trait Estimate (Full Data)")
cor.test(traits.comb$mean.y,traits.comb$mean.x)
ggsave("plots/effect of subsetting on trait effect estimates.jpg",bg="white",height = 4,width = 5)



