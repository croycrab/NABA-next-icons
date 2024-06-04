#phylogeny plot


# setup -------------------------------------------------------------------
# load some libraries
library(phytools)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)

load("traits and phylo.RData")

# plot tree with trends ---------------------------------------------------

tips.rm3 = tree.full$tip.label[which(!(tree.full$tip.label %in% df.global$Species))]
tree.sub3 = tree.full %>% drop.tip(tips.rm3) %>% as.phylo()

trends = df.global[match(tree.sub3$tip.label, df.global$Species),] %>% 
  dplyr::select(Species,mean_trend,ucl_trend,lcl_trend,Family,diff_zero) %>% 
  left_join(df.traits2[,c("Species","Subfamily")], by = "Species")
# filter(is.na(Family) == F & mean_trend < 15 & mean_trend > -15)

ggtree(tree.sub3,layout = "circular") %<+% trends +
  geom_tippoint(size = .5,aes(color = Family))+
  scale_color_viridis_d("Family",begin = 0,end = .9,)+
  geom_tiplab(size=1.5,align = T,aes(color = Family))+
  geom_rootedge()+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))+
  theme_tree2(legend.position = "top",legend.direction = "horizontal")+
  xlim(c(0,200))
ggsave("plots/phylogeny/phylogeny of species modeled.jpg",bg = "white",height = 12, width = 12)