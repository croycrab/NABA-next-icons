# objective: get 10 year change rates for each species

# set up -----------------------------------------------------------------------
# load some libraries
library(tidyverse); library(ggridges);library(ggpubr)


# load trait data ---------------------------------------------------------
load("traits and phylo.RData")


# estimate range size from strata modeled ---------------------------------
sp.list = df.global$Species
df.widespread = data.frame(Species = sp.list)
for(i in 1:length(sp.list)){
    df.widespread$strat_area[i] = read_csv(paste0("data/site strat key/",sp.list[i],"_site_strat_key.csv")) %>% 
    filter(strat_z > 0) %>% 
    group_by(strat_idx1) %>% 
    reframe(strat_A = mean(strat_A)) %>% 
    pull(strat_A) %>% 
    sum()
}

# view trends of widespread species ---------------------------------------
df.global %>% 
  arrange(-range_size) %>% 
  slice(1:75) %>% 
  pull(Species) -> widespread.sp
  
df.widespread %>% 
  arrange(-strat_area) %>% 
  slice(1:75) %>% 
  pull(Species) -> widespread.sp2

intersect(widespread.sp,widespread.sp2) %>% length() #some minor discrepancies between range sizes estimated from brood maps vs from the areas we modeled -- probably worth using definition of "widespread" using the area that we were able to model trends for. I think this serves the additional benefit of highlighting trends for species that are both relatively widespread and for which we had decent spatial coverage for. 

# ggplot(df.global %>% filter(Species %in% widespread.sp), aes(x = reorder(Species,mean_trend,decreasing=T), y = mean_trend, col = as.factor(sign(mean_trend)))) +
#   geom_point(size=2)+
#   geom_errorbar(aes(ymin = lcl_trend, ymax = ucl_trend), position = "dodge", width = 0.25)+
#   scale_color_manual(values=c("-1" = "red","1" = "black"),name = "",labels=c("-1" = "-","1" = "+")) +
#   labs(title = '',x='',y="Trend (%/yr)")+
#   geom_hline(yintercept = 0,color="black",linetype="dashed")+
#   theme_classic()+
#   theme(axis.text.y = element_text(face = "italic",color = "black"),axis.title.x = element_text(hjust = .5),
#         panel.grid.major.y = element_line(color = "gray",linewidth = 0.15))+
#   coord_flip()


#define IUCN limits
cr_lim = log(2)/(abs((log(.2)/10)))
vu_lim = log(2)/(abs((log(.7)/10)))
trend.df3 <- read.csv("data/trends_all species_all draws.csv") %>% 
  group_by(species) %>% 
  mutate(q2.5 = quantile(slope,prob=.025),
         q97.5= quantile(slope,prob=.975),
         q50= quantile(slope,prob=.50),
         q25=quantile(slope,prob=.25),
         q75=quantile(slope,prob=.75),
         q90=quantile(slope,prob=.9),
         q10=quantile(slope,prob=.1)) %>% 
  ungroup() %>% 
  filter(slope < q75 & slope > q25) %>% 
  mutate(half_life = (log(2)/(abs(slope)/100)),q25_hl = (log(2)/(abs(q25)/100)),iucn_lims = ifelse(q25_hl<cr_lim,"CR",ifelse(cr_lim < q25_hl & q25_hl < 10,"EN",ifelse(10 < q25_hl & q25_hl < vu_lim,"VU","NT"))),
         species = case_when(species == "Colias eurydice" ~ "Zerene eurydice",
                             species == "Everes amyntula" ~ "Cupidio amyntula",
                             species == "Megisto rubricata" ~ "Cissia rubricata",
                             species == "Enodia anthedon" ~ "Lethe anthedon",
                             .default = species))

ggplot(trend.df3 %>% 
         filter(q75 < 0 & species %in% widespread.sp2) %>% 
         left_join(df.widespread[,c("Species","strat_area")],by = c("species" = "Species")), aes(x = half_life, y = reorder(species,strat_area,decreasing = F), fill = factor(iucn_lims,levels = c("CR","EN","VU","NT")))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_d("") +
  labs(title = '',y='',x="Half Life (years)")+
  # geom_vline(xintercept = c(cr_lim,10,vu_lim),color=c("red4","firebrick","darkorange"),linetype="solid",size = 1.25)+
  theme_classic()+
  theme_ridges(grid = T,font_size = 11)+
  theme(axis.text.y = element_text(face = "italic"),axis.title.x = element_text(hjust = .5),legend.position = c(.75,.85),legend.background = element_rect(fill = "white"))+
  lims(x=c(0,150))
ggsave("plots/IUCN criteria and widespread species.png", bg = "white")

#get top 10 widely distributed species that also appear to be declining
trend.df3 %>% 
  filter(q75 < 0 & species %in% widespread.sp2) %>% 
  group_by(species) %>% 
  reframe(trend = mean(slope)) %>% 
  ungroup() %>% 
  left_join(df.widespread[,c("Species","strat_area")],by = c("species" = "Species")) %>% 
  select(species,strat_area,trend) %>% 
  arrange(-strat_area) %>% 
  select(species) %>%
  slice(1:30) %>% 
  print(n = 30)

trend.df3 %>% 
  filter(q75 < 0 & species %in% widespread.sp2) %>% 
  group_by(species) %>% 
  reframe(trend = mean(slope)) %>% 
  ungroup() %>% 
  left_join(df.widespread[,c("Species","strat_area")],by = c("species" = "Species")) %>% 
  select(species,strat_area,trend) %>% 
  arrange(-strat_area) %>% 
  slice(1:10) %>% 
  arrange(trend) %>% 
  print(n=30)

trend.df3 %>% 
  filter(q75 < 0 & species %in% widespread.sp2) %>% 
  group_by(species) %>% 
  reframe(trend = mean(slope)) %>% 
  ungroup() %>% 
  left_join(df.widespread[,c("Species","strat_area")],by = c("species" = "Species")) %>% 
  select(species,strat_area,trend) %>% 
  arrange(-strat_area) %>% 
  pull(species) -> wide.dec

# view distributions of traits for widespread declining species -----------
gghistogram(df.global,x = "canopy2",stat = "count",fill = "black") +
  labs(y = "Count",x = "Canopy Preference")
gghistogram(df.global %>% filter(Species %in% wide.dec),x = "canopy2",stat = "count",fill = "black") +
  labs(y = "Count",x = "Canopy Preference")

(sum(df.global$canopy %in% c("Open (strong)","Open (weak)"),na.rm = T)/sum(!is.na(df.global$canopy),na.rm = T))*100
(sum(df.global$canopy[which(df.global$Species %in% wide.dec)] %in% c("Open (strong)","Open (weak)"),na.rm = T)/sum(!is.na(df.global$canopy[which(df.global$Species %in% wide.dec)]),na.rm = T))*100

gghistogram(df.global,x = "voltinism",stat = "count",fill = "black") +
  labs(y = "Count",x = "Canopy Preference")
gghistogram(df.global %>% filter(Species %in% wide.dec),x = "voltinism",stat = "count",fill = "black") +
  labs(y = "Count",x = "Canopy Preference") 
# REVIEW THIS CODE BELOW--not sure if it is necessary, most likely needs to be deleted
  
# get all 5000 slope estimates --------------------------------------------
sp.vect <- str_split(list.files("data/sampled from model_strata"),"_",simplify = T)[,1]
trend.list <- list()
for(i in 1:length(sp.vect)){
  d0 <- read.csv(paste0("data/sampled from model_strata/",sp.vect[i],"_samps3.csv"))
  if("strat_nAz" %in% colnames(d0)){
    trend.list[[i]] <- d0 %>% 
      mutate(species = sp.vect[i]) %>% 
      group_by(.draw,year) %>% 
      reframe(count = sum(strat_nAz))%>%
      ungroup() %>% 
      group_by(.draw) %>%
      do(mod = glm(round(count, 0) ~ year, data = ., family="poisson")) %>%
      mutate(slope = summary(mod)$coeff[2],slope = ((exp(slope)-1)*100)) %>%
      select(-mod) %>% ungroup()

  } else( 
    trend.list[[i]] <- d0 %>%
      mutate(species = sp.vect[i],count = aggregate_N)%>%
      group_by(.draw) %>%
      do(mod = glm(round(count, 0) ~ year, data = ., family="poisson")) %>%
      mutate(slope = summary(mod)$coeff[2],slope = ((exp(slope)-1)*100)) %>%
      select(-mod) %>% ungroup()
  )
}
for(i in 1:length(sp.vect)){
  trend.list[[i]]$species <- sp.vect[i]
  
}
trend.df<-do.call("rbind",trend.list)
write.csv(trend.df,"data/trends_all species_all draws.csv")


# prepare trend data ------------------------------------------------------
trend.df2 <- read.csv("data/trends_all species_all draws.csv") %>%
  group_by(species) %>%
  mutate(q2.5 = quantile(slope,prob=.025),q97.5= quantile(slope,prob=.975),q50= quantile(slope,prob=.50),q25=quantile(slope,prob=.25),q75=quantile(slope,prob=.75)) %>%
  ungroup() %>%
  filter(slope < q75 & slope > q25) %>%
  mutate(half_life = (log(2)/(abs(slope)/100)),
         species = case_when(species == "Colias eurydice" ~ "Zerene eurydice",
                             species == "Everes amyntula" ~ "Cupidio amyntula",
                             species == "Megisto rubricata" ~ "Cissia rubricata",
                             species == "Enodia anthedon" ~ "Lethe anthedon",
                             .default = species))

# ggplot(trend.df2 %>% filter(q75 < 0), aes(x = log(half_life), y = reorder(species,q50,decreasing = T), fill = q50)) +
#   geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
#   scale_fill_viridis_c("") +
#   labs(title = '',y='',x="Half Life")+
#   guides(fill = "none")+
#   # scale_x_continuous(labels = scales::percent_format(scale = 1))+
#   geom_vline(xintercept = log(10),color="red",linetype="solid")+
#   theme_classic()+
#   theme_ridges(grid = T,font_size = 11)+
#   theme(axis.text.y = element_text(face = "italic"),axis.title.x = element_text(hjust = .5))
# # lims(x=c(-0.35,0.35))
# ggsave('./plots/half life_ridges.jpg',bg="white",height = 9,width = 5)

# #repeat but for 2008-2018 trends
# trend.list.2008 <- list()
# for(i in 1:length(sp.vect)){
#   d0 <- read.csv(paste0("data/sampled from model_strata/",sp.vect[i],"_samps3.csv"))
#   if("strat_nAz" %in% colnames(d0)){
#     trend.list.2008[[i]] <- d0 %>%
#       filter(year>=2008) %>%
#       mutate(species = sp.vect[i],count = strat_nAz)%>%
#       group_by(.draw) %>%
#       do(mod = glm(round(count, 0) ~ year, data = ., family="poisson")) %>%
#       mutate(slope = summary(mod)$coeff[2],slope = ((exp(slope)-1)*100)) %>%
#       select(-mod) %>%
#       ungroup()
# 
#   } else(
#     trend.list.2008[[i]] <- d0 %>%
#       filter(year>=2008) %>%
#       mutate(species = sp.vect[i],count = aggregate_N)%>%
#       group_by(.draw) %>%
#       do(mod = glm(round(count, 0) ~ year, data = ., family="poisson")) %>%
#       mutate(slope = summary(mod)$coeff[2],slope = ((exp(slope)-1)*100)) %>%
#       select(-mod) %>%
#       ungroup()
#     )
#   trend.list.2008[[i]]$species <- sp.vect[i]
# }
# # for(i in 1:length(sp.vect)){
# #   trend.list.2008[[i]]$species <- sp.vect[i]
# # }
# trend.2008.df<-do.call("rbind",trend.list.2008)
# write.csv(trend.2008.df,"data/trends_2008_all species_all draws.csv")


#define IUCN limits
cr_lim = log(2)/(abs((log(.2)/10)))
vu_lim = log(2)/(abs((log(.7)/10)))
trend.df3 <- read.csv("data/trends_all species_all draws.csv") %>% 
  group_by(species) %>% 
  mutate(q2.5 = quantile(slope,prob=.025),q97.5= quantile(slope,prob=.975),q50= quantile(slope,prob=.50),q25=quantile(slope,prob=.25),q75=quantile(slope,prob=.75)) %>% 
  ungroup() %>% 
  filter(slope < q75 & slope > q25) %>% 
  mutate(half_life = (log(2)/(abs(slope)/100)),q25_hl = (log(2)/(abs(q25)/100)),iucn_lims = ifelse(q25_hl<cr_lim,"CR",ifelse(cr_lim < q25_hl & q25_hl < 10,"EN",ifelse(10 < q25_hl & q25_hl < vu_lim,"VU","NT"))),
         species = case_when(species == "Colias eurydice" ~ "Zerene eurydice",
                             species == "Everes amyntula" ~ "Cupidio amyntula",
                             species == "Megisto rubricata" ~ "Cissia rubricata",
                             species == "Enodia anthedon" ~ "Lethe anthedon",
                             .default = species))

ggplot(trend.df3 %>% filter(species %in% trend.df2$species[which(trend.df2$q97.5<0)] & slope > q25 & slope < q75), aes(x = half_life, y = reorder(species,q50,decreasing = T), fill = factor(iucn_lims,levels = c("CR","EN","VU","NT")))) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_d("") +
  labs(title = '',y='',x="Half Life (years)")+
  # geom_vline(xintercept = c(cr_lim,10,vu_lim),color=c("red4","firebrick","darkorange"),linetype="solid",size = 1.25)+
  theme_classic()+
  theme_ridges(grid = T,font_size = 11)+
  theme(axis.text.y = element_text(face = "italic"),axis.title.x = element_text(hjust = .5),legend.position = c(.75,.85),legend.background = element_rect(fill = "white"))
  # lims(x=c(0,50))
ggsave('plots/halflife from 10 year change_ridges_declining species.jpg',bg="white")

# ggplot(tyc.df  %>% filter(species %in% trend.df2$species[which(trend.df2$q97.5<0)]), aes(x = TYC, y = reorder(species,q50,decreasing = T), fill = q50)) +
#   geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
#   scale_fill_viridis_c("") +
#   labs(title = '',y='',x="Ten Year Change")+
#   guides(fill = "none")+
#   scale_x_continuous(labels = scales::percent_format(scale = 1))+
#   geom_vline(xintercept = c(-80,-50,-30),color="red",linetype="solid")+
#   theme_classic()+
#   theme_ridges(grid = T,font_size = 11)+
#   theme(axis.text.y = element_text(face = "italic"),axis.title.x = element_text(hjust = .5))
# # lims(x=c(-0.35,0.35))
# ggsave('./plots/ten year change_ridges_declining species.jpg',bg="white",height = 15,width = 5)
# 

#create table of iucn/esa status
#create table of relevant information
sp.trends <- read.csv("data/trends_1996-2018.csv")
iucn.list <- read.csv("data/IUCN and ESA/assessments.csv")
iucn.list2 <- iucn.list %>% 
  rename(Species = scientificName,status = redlistCategory) %>% 
  mutate(source = "IUCN") %>% 
  select(Species,status,source,populationTrend,yearPublished) %>% 
  filter(Species %in% sp.trends$species) %>% 
  group_by(Species) %>% 
  filter(yearPublished == max(yearPublished))

esa.list <- read.csv("data/IUCN and ESA/esa animals.csv") %>% 
  filter(Species.Group %in% "Insects")
rows.keep<-unique(c(grep("butterfly",esa.list$Common.Name),
                    grep("Butterfly",esa.list$Common.Name),
                    grep("skipper",esa.list$Common.Name),
                    grep("skipperling",esa.list$Common.Name),
                    grep("Skipper",esa.list$Common.Name)))
esa.list2 <- esa.list[rows.keep,] %>% 
  mutate(genus = str_split(Scientific.Name," ",simplify = T)[,1],
         species = ifelse(str_detect(str_split(Scientific.Name," ",simplify = T)[,2],"^[a-zA-Z0-9]+")
                          ,str_split(Scientific.Name," ",simplify = T)[,2],str_split(Scientific.Name," ",simplify = T)[,3]),
         subspecies = ifelse(str_detect(str_split(Scientific.Name," ",simplify = T)[,2],"^[a-zA-Z0-9]+")
                             ,str_split(Scientific.Name," ",simplify = T)[,3],str_split(Scientific.Name," ",simplify = T)[,4]),
         Species = str_trim(paste(genus,species," ")), source = "ESA") %>% 
  select(Species,subspecies,Federal.Listing.Status,Where.Listed,source) %>% 
  rename(status = Federal.Listing.Status,location_listed = Where.Listed) %>% 
  filter(Species %in% sp.trends$species & status != "Experimental Population, Non-Essential")

tyc.sum <- tyc.df %>% 
  group_by(species) %>% 
  summarize(med_tyc = median(TYC),q975_TYC = quantile(TYC,prob=0.975),q025_TYC = quantile(TYC,prob=0.025))

trends_10yr.sum <- trend.df3 %>% 
  group_by(species) %>% 
  summarize(med_trend = mean(q50),q975_trend = mean(q97.5),q025_trend = mean(q2.5))

status.list <- bind_rows(esa.list2,iucn.list2) %>% 
  left_join(sp.trends[,c("species","mean_trend","lcl_trend","ucl_trend","diff_zero")],by=c("Species"="species")) %>% 
  left_join(tyc.sum[,c("species","med_tyc","q975_TYC","q025_TYC")],by=c("Species" = "species")) %>% 
  left_join(trends_10yr.sum[,c("species","med_trend","q975_trend","q025_trend")],by=c("Species" = "species")) %>% 
  mutate(agreement = ifelse((source == "ESA" & mean_trend < 0) | (source == "IUCN" & mean_trend <0 & populationTrend == "Decreasing"), "yes","no"))

sum(status.list[which(!(status.list$populationTrend %in% c("Stable","Unknown"))),"agreement"] == "yes")/length(status.list[which(!(status.list$populationTrend %in% c("Stable","Unknown"))),"agreement"])


status_decline <- sp.trends %>% 
  filter(ucl_trend < 0) %>% 
  left_join(status.list[,c("Species","status","source","populationTrend")],by=c("species"="Species")) %>% 
  left_join(trends_10yr.sum[,c("species","med_trend","q975_trend","q025_trend")],by="species") %>% 
  arrange(med_trend) %>% 
  select(species,mean_trend,lcl_trend,ucl_trend,med_trend,q975_trend,q025_trend,status,source,populationTrend) %>% 
  mutate_if(is.numeric,round,2)
write.csv(status_decline,"species of concern/conservation status_declining species.csv")

#categories
df.classify <- trend.df3 %>% 
  filter(species %in% unique(trend.df2$species[which(trend.df2$q97.5 < 0)])) %>% 
  select(species,q25,q75) %>% 
  distinct() %>% 
  mutate(CR = ifelse((log(2)/(abs(q25)/100)) < log(2)/(abs((log(.2)/10))),1,0),
        EN = ifelse((log(2)/(abs(q25)/100)) < 10,1,0),
        VU = ifelse((log(2)/(abs(q25)/100)) < log(2)/(abs((log(.7)/10))),1,0))

sum(df.classify$CR) #14
sum(df.classify$EN)-sum(df.classify$CR) #14
sum(df.classify$VU)-sum(df.classify$EN) #11

sum(trend.df3$half_life[which(trend.df3$species == "Cercyonis pegala")] <= 50)/2500
sum(trend.df3$half_life[which(trend.df3$species == "Asterocampa celtis")] <= 50)/2500
sum(trend.df3$half_life[which(trend.df3$species == "Junonia genoveva")] <= 50)/2500
sum(trend.df3$half_life[which(trend.df3$species == "Megisto cymela")] <= 50)/2500


