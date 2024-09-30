#objective: extract land use data for bbs strata
#load packages
library(dplyr); library(raster); library(sf); library(bbsBayes); library(ncdf4); library(tidyr); library(stringr); library(RCurl); library(XML); library(ncdf4)

#get strata that contain NABA sites
#load bbs stratum shapes
map1 <- load_map(stratify_by="bbs_usgs") %>% rename_all(tolower) %>% 
  dplyr::select(2, 3) %>% rename(strat_name=st_12, strat_A=area_1)
#load butterfly site data
sites1 <- 
  read.table("data/raw/butterfly_data_NoMerge.txt",as.is=T,check.names=F,header=T)%>%
  rename_all(tolower)%>%
  filter(between(julian.date, 152, 243))%>%
  distinct(site,longitude,latitude) %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>%
  st_transform(crs=st_crs(map1))

# pair each site with its respective bbs stratum
within1 <- sapply(st_within(sites1, map1), function(z) if (length(z)==0) 
  NA_integer_ else z[1])
near1 <- st_nearest_feature(sites1, map1) 
poly_per_site <- ifelse(is.na(within1), near1, within1)
strata <- sites1 %>% 
  mutate(strat_name=unique(map1$strat_name)[poly_per_site]) %>% 
  distinct(strat_name)

map2<-map1 %>% filter(strat_name %in% strata$strat_name) %>% st_transform(map1, crs = "+proj=longlat +datum=WGS84 +no_defs ")

#load land use data
lulc.brick <- brick("data/land use/hildap_vGLOB-1.0-f_netcdf/hildaplus_GLOB-1-0-f_states.nc",varname = "LULC_states") 
lulc.list.yr <- list()
lulc.list.strat <- list()
for(i in 98:120){ #I think there was an error and so I resumed from the error onward --change to start and end index for strata extractions
  for(j in 1:length(map2$strat_name)){
    lulc.list.strat[[j]] <-table(raster::extract(x=lulc.brick[[i]],y=map2[j,],df=T,na.rm=T)) %>% as.data.frame() %>% mutate(strat_name = map2[j,1]$strat_name,lulc_state = .[,2],year = lulc.brick[[i]]@z[[1]]) %>% rename_all(tolower) %>% dplyr::select(strat_name,year,lulc_state,freq)
  }
  lulc.list.yr[[i-97]] <- do.call("rbind",lulc.list.strat)
}
lulc.df <- do.call("rbind",lulc.list.yr)
write.csv(lulc.df,"data/land use/lulc_by year and strata.csv")
