# set up -----------------------------------------------------------------------
# load some libraries
library(ggplot2);library(dplyr); library(tidyr)

# expand, append, insert zeros, and save ---------------------------------------
butterflies <- 
  read.table("data/raw/butterfly_data_NoMerge.txt",as.is=T,check.names=F,header=T)%>%
  filter(between(Julian.date, 152, 243))%>%
  arrange(Site,Year)

butterflies_exp <- 
  read.table("data/raw/butterfly_data_NoMerge.txt",as.is=T,check.names=F,header=T) %>%
  filter(between(Julian.date, 152, 243))%>%
  expand(nesting(Site,Year),Species)%>%
  arrange(Site,Year,Species)

but.sum = butterflies[,c(1,2,6:9)] %>%
  group_by(Site,Year)%>%
  summarize_all(.funs = "unique")

sp.list <- unique(butterflies$Species)
for(i in 1:length(sp.list)){
  butterflies_exp %>% filter(Species == sp.list[i]) %>%
    left_join(but.sum,by=c("Site","Year"))%>%
    left_join(butterflies[which(butterflies$Species == sp.list[i]),names(butterflies)[1:4]],by=c("Site","Year","Species"))%>%
    replace_na(list(N.butterflies = 0))%>%
    rename(hours=Party_Hours, j_date=Julian.date, count = N.butterflies, lat = Latitude, lon = Longitude) %>%
    rename_all(tolower) %>% 
    write.csv(paste0("data/zeros added/",sp.list[i],"_zeros added.csv"),row.names = F)
}