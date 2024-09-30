#objective: extract climate data for bbs strata
#load packages
library(dplyr); library(raster); library(sf); library(bbsBayes); library(ncdf4); library(sf); library(tidyr); library(stringr)

#get strata that contain NABA sites
#load bbs stratum shapes
map1 <- load_map(stratify_by="bbs_usgs") %>% rename_all(tolower) %>% 
  dplyr::select(2, 3) %>% rename(strat_name=st_12, strat_A=area_1)
#load butterfly site data
sites1 <- 
  read.table("./data/butterfly_data_NoMerge.txt",as.is=T,check.names=F,header=T)%>%
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
  
# #download climate data
# #download data
# vars <- c("cld","dtr","frs","pet","pre","tmp","tmn","tmx","vap","wet")
# years<-c("1991","2001","2011")
# 
# for(i in 1:length(vars)){
#   clim.url <- paste0("http://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.06/cruts.2205201912.v4.06/",vars[i],"/")
#   clim.curl <- getURL(clim.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
#   file.vector <- getHTMLLinks(clim.curl)
#   for(k in 1:length(years)){
#     file.ncs<-file.vector[grep(".nc", file.vector)]
#     download.file(url = paste0(clim.url,file.ncs[grep(years[k],file.ncs)]),destfile = paste0("data/climate/raw/",file.ncs[grep(years[k],file.ncs)]))
#   }
# }

#load climate data 
vars <- c("cld","dtr","frs","pet","pre","tmp","tmn","tmx","vap","wet")
for(i in 1:length(vars)){
  var.files <- list.files("climate/raw")[grep(vars[i],list.files("climate/raw"))]
  assign(paste0(vars[i],".stack"),stack(brick(paste0("data/climate/raw/",var.files[1])),brick(paste0("data/climate/raw/",var.files[2])),brick(paste0("data/climate/raw/",var.files[3]))))
}

#subset map to include only strata with butterfly census circles
map2<-map1 %>% filter(strat_name %in% strata$strat_name) %>% st_transform(map1, crs = st_crs(tmp.stack))

#extract mean for each strata
for(i in 1:length(vars)){
  assign(paste0(vars[i],".strat.df"),raster::extract(x=get(paste0(vars[i],".stack"), envir = .GlobalEnv),y=map2,fun=mean,df=T) %>% mutate(strat_name = map2$strat_name) %>% relocate(strat_name))
}

#select months/years of interest only
#convert wide to long
#get average climate for each year
cld.long <- gather(cld.strat.df,date,cld,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(cld = mean(cld)) %>%
  filter(year2 %in% c(1996:2018))

dtr.long <- gather(dtr.strat.df,date,dtr,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(dtr = mean(dtr)) %>%
  filter(year2 %in% c(1996:2018))

frs.long <- gather(frs.strat.df,date,frs,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(frs = mean(frs)) %>%
  filter(year2 %in% c(1996:2018))

pet.long <- gather(pet.strat.df,date,pet,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(pet = sum(pet)) %>%
  filter(year2 %in% c(1996:2018))

pre.long <- gather(pre.strat.df,date,pre,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(pre = sum(pre)) %>%
  filter(year2 %in% c(1996:2018))

tmp.long <- gather(tmp.strat.df,date,tmp,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(tmp = mean(tmp)) %>%
  filter(year2 %in% c(1996:2018))

tmn.long <- gather(tmn.strat.df,date,tmn,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(tmn = mean(tmn)) %>%
  filter(year2 %in% c(1996:2018))

tmx.long <- gather(tmx.strat.df,date,tmx,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(tmx = mean(tmx)) %>%
  filter(year2 %in% c(1996:2018))

vap.long <- gather(vap.strat.df,date,vap,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(vap = mean(vap)) %>%
  filter(year2 %in% c(1996:2018))

wet.long <- gather(wet.strat.df,date,wet,X1991.01.16:X2020.12.16,factor_key=TRUE) %>% 
  mutate(year=as.numeric(substr(str_split(date,"X", simplify=T)[,2],1,4)),
         month=as.numeric(substr(str_split(date,"X", simplify=T)[,2],6,7)),
         year2=ifelse(month %in% c(7:12),year+1,year))%>%
  dplyr::select(-date)%>%
  group_by(strat_name,year2)%>%
  summarise(wet = mean(wet)) %>%
  filter(year2 %in% c(1996:2018))


#save files
clim.df <- cld.long %>% left_join(dtr.long) %>% left_join(frs.long) %>% left_join(pet.long) %>% left_join(pre.long) %>% left_join(tmp.long) %>% left_join(tmn.long) %>% left_join(tmx.long) %>% left_join(vap.long) %>% left_join(wet.long) %>% rename(year = year2) 
clim.df %>% write.csv(.,"data/climate/climate by strata and year.csv")

#--------------------------#plot climate trends#------------------------------#
library(ggplot2)
clim.df <- read.csv("data/climate/climate by strata and year.csv")
# map climate data 
clim.df.long <- clim.df %>% gather(.,clim_var,value,cld:wet,factor_key = TRUE) 
map3 <- map2 %>% rename_all(tolower) %>% 
  right_join(clim.df.long) 
ggplot() + 
  geom_sf(data=map3%>%filter(year == 2018), aes(fill=value)) + 
  # geom_sf(data=sites2, alpha=0.4, aes(size=n_years), pch=16) +
  # scale_size_continuous("number\nyears", range=c(0.2, 2), guide="none") +
  scale_fill_scico(palette = "vik") + 
  theme_bw()+
  facet_wrap(~clim_var,dir = "h")
