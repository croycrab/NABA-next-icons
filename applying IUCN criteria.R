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

