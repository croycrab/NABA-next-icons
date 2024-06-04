#examine and plot results from top model
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


# examine top model -------------------------------------------------------
#define priors
bym_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                  phi = list(prior = "pc", param = c(0.5, 0.5)))
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

#load inla graph and abundance data (with covariates)
#load all models, graph, and data
load("model selection list.RData")
g1 <- inla.read.graph(filename = "neighborhood graphs/map_allspecies.adj")
d1 <- read.csv("data/all species counts_by strata and year_with covariates.csv") %>% 
  left_join(sp.order,by = c("species" = "Species"))
d1$strat_idx2 = d1$strat_idx3 = d1$strat_idx4 = d1$strat_idx5 = d1$strat_idx1
d1$sp_idx2 = d1$sp_idx3 = d1$sp_idx4 = d1$sp_idx5 = d1$sp_idx1 
# d1$family_idx1 = as.integer(factor(d1$Family))
d1<-d1 %>% left_join(df.strata2[,c("species","strat_name","brood")],by = c("species","strat_name"))

top.model <- inla(mod.list[[53]], data=d1,scale = log_prec_y,
                  family = "gaussian",
                  control.predictor=list(compute=T, link=1),
                  control.compute=list(config=T, waic = T, cpo=T,
                                       control.gcpo=list(enable=T)),
                  verbose=T, inla.mode="experimental",
                  control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
                  control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))


# species level covariate effects -----------------------------------------
length(unique(d1$strat_idx1))
mod.list[[53]]
#most/least susceptible to warming
top.model$summary.random$sp_idx2[which(top.model$summary.random$sp_idx2$mean == min(top.model$summary.random$sp_idx2$mean)),"ID"]
d1$species[which(d1$sp_idx2 == "810")] %>% unique() #Thymelicus lineola

top.model$summary.random$sp_idx2[which(top.model$summary.random$sp_idx2$mean == max(top.model$summary.random$sp_idx2$mean)),"ID"]
d1$species[which(d1$sp_idx2 == "275")] %>% unique() #Hemiargus ceraunus

#most/least susceptible to drying
top.model$summary.random$sp_idx3[which(top.model$summary.random$sp_idx3$mean == min(top.model$summary.random$sp_idx3$mean)),"ID"]
d1$species[which(d1$sp_idx3 == "514")] %>% unique() #Speyeria cybele

top.model$summary.random$sp_idx3[which(top.model$summary.random$sp_idx3$mean == max(top.model$summary.random$sp_idx3$mean)),"ID"]
d1$species[which(d1$sp_idx3 == "115")] %>% unique() #Ascia monuste

#most/least susceptible to urbanization
top.model$summary.random$sp_idx4[which(top.model$summary.random$sp_idx4$mean == min(top.model$summary.random$sp_idx4$mean)),"ID"]
d1$species[which(d1$sp_idx4 == "482")] %>% unique() #Chlosyne harrisii

top.model$summary.random$sp_idx4[which(top.model$summary.random$sp_idx4$mean == max(top.model$summary.random$sp_idx4$mean)),"ID"]
d1$species[which(d1$sp_idx4 == "798")] %>% unique() #Cymaenes tripunctus

#most/least susceptible to ag expansion
top.model$summary.random$sp_idx5[which(top.model$summary.random$sp_idx5$mean == min(top.model$summary.random$sp_idx5$mean)),"ID"]
d1$species[which(d1$sp_idx5 == "300")] %>% unique() #Calephelis arizonensis

top.model$summary.random$sp_idx5[which(top.model$summary.random$sp_idx5$mean == max(top.model$summary.random$sp_idx5$mean)),"ID"]
d1$species[which(d1$sp_idx5 == "781")] %>% unique() #Amblyscirtes texanae

#pieris rapae covariate effects
sp.order$sp_idx1[which(sp.order$Species == "Pieris rapae")]
top.model$summary.fixed
top.model$summary.random$sp_idx2[which(top.model$summary.random$sp_idx2$ID == "107"),] #warming
top.model$summary.random$sp_idx3[which(top.model$summary.random$sp_idx3$ID == "107"),] #drying
top.model$summary.random$sp_idx4[which(top.model$summary.random$sp_idx4$ID == "107"),] #urbanization
top.model$summary.random$sp_idx5[which(top.model$summary.random$sp_idx5$ID == "107"),] #agricultural expansion
hist(top.model$summary.random$sp_idx2$mean)


p.rapae.tmp <- gghistogram(top.model$summary.random$species_idx2$mean+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "tmp")],fill = "gray")+
  annotate("rect",xmin = (top.model$summary.random$species_idx2$`0.025quant`[which(top.model$summary.random$species_idx2$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "tmp")]),xmax = (top.model$summary.random$species_idx2$`0.975quant`[which(top.model$summary.random$species_idx2$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "tmp")]),ymin = 0,ymax = Inf,fill="blue", alpha = .5)+
  # geom_vline(xintercept = 0,color="black",size = 1.25,lty = "solid")+
  labs(y= "Frequency", x="Temperature Effect (slope)")

p.rapae.pre <- gghistogram(top.model$summary.random$species_idx3$mean+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "pre")],fill = "gray")+
  annotate("rect",xmin = (top.model$summary.random$species_idx3$`0.025quant`[which(top.model$summary.random$species_idx3$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "pre")]),xmax = (top.model$summary.random$species_idx3$`0.975quant`[which(top.model$summary.random$species_idx3$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "pre")]),ymin = 0,ymax = Inf,fill="blue", alpha = .5)+
  # geom_vline(xintercept = 0,color="black",size = 1.25,lty = "solid")+
  labs(y= "Frequency", x="Precipitation Effect (slope)")

p.rapae.urb <- gghistogram(top.model$summary.random$species_idx4$mean+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "urban")],fill = "gray")+
  annotate("rect",xmin = (top.model$summary.random$species_idx4$`0.025quant`[which(top.model$summary.random$species_idx4$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "urban")]),xmax = (top.model$summary.random$species_idx4$`0.975quant`[which(top.model$summary.random$species_idx4$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "urban")]),ymin = 0,ymax = Inf,fill="blue", alpha = .5)+
  # geom_vline(xintercept = 0,color="black",size = 1.25,lty = "solid")+
  labs(y= "Frequency", x="Urbanization Effect (slope)")

p.rapae.ag <- gghistogram(top.model$summary.random$species_idx5$mean+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "agriculture")],fill = "gray")+
  annotate("rect",xmin = (top.model$summary.random$species_idx5$`0.025quant`[which(top.model$summary.random$species_idx5$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "agriculture")]),xmax = (top.model$summary.random$species_idx5$`0.975quant`[which(top.model$summary.random$species_idx5$ID == "225")]+top.model$summary.fixed$mean[which(rownames(top.model$summary.fixed) == "agriculture")]),ymin = 0,ymax = Inf,fill="blue", alpha = .5)+
  # geom_vline(xintercept = 0,color="black",size = 1.25,lty = "solid")+
  labs(y= "Frequency", x="Agriculture Effect (slope)")

ggpubr::ggarrange(p.rapae.tmp,p.rapae.pre,p.rapae.urb,p.rapae.ag)
ggsave("plots/pieris rapae vs other species.jpg",bg="white",width = 7.5,height = 5)
ggsave("plots/pieris rapae vs other species.pdf",bg="white",width = 7.5,height = 5)

# get model predictions ---------------------------------------------------
#set up climate predictions
all.pred <- data.frame(tmp = as.data.frame(quantile(unique(d1$tmp_mean),probs = c(0.05,0.5,0.95)))[,1],pre = as.data.frame(quantile(unique(d1$pre),probs = c(0.05,0.5,0.95)))[,1],urban = as.data.frame(quantile(unique(d1$urban),probs = c(0.05,0.5,0.95)))[,1],agriculture = as.data.frame(quantile(unique(d1$agriculture),probs = c(0.05,0.5,0.95)))[,1],tmp_mean = as.data.frame(quantile(unique(d1$tmp_mean),probs = c(0.05,0.5,0.95)))[,1],pre_mean = as.data.frame(quantile(unique(d1$pre_mean),probs = c(0.05,0.5,0.95)))[,1],  urb_mean = as.data.frame(quantile(unique(d1$urb_mean),probs = c(0.05,0.5,0.95)))[,1],ag_mean = as.data.frame(quantile(unique(d1$ag_mean),probs = c(0.05,0.5,0.95)))[,1],brood = c("1","2","3")) %>% 
  mutate(tmp_idx = as.integer(as.factor(tmp)),pre_idx = as.integer(as.factor(pre)),urb_idx = as.integer(as.factor(urban)),ag_idx = as.integer(as.factor(agriculture)),tmp_mean_idx = as.integer(as.factor(tmp_mean)),pre_mean_idx = as.integer(as.factor(pre_mean)),urb_mean_idx = as.integer(as.factor(urb_mean)),ag_mean_idx = as.integer(as.factor(ag_mean)),volt_idx = as.integer(as.factor(brood)))

#create all combinations that we need for prediction
clim_combs <- all.pred %>% 
  mutate(urb_idx = 2,urb_mean_idx = 2,ag_idx = 2,ag_mean_idx = 2) %>% 
  expand(tmp_idx,pre_idx,tmp_mean_idx,pre_mean_idx,urb_idx,urb_mean_idx,ag_idx,ag_mean_idx,volt_idx)

# clim_combs2 <- all.pred %>% 
#   mutate(tmp_idx = 2, tmp_mean_idx = 2,urb_idx = 2,urb_mean_idx = 2,ag_idx = 2,ag_idx = 2) %>% 
#   expand(tmp_idx,pre_idx,tmp_mean_idx,pre_mean_idx,urb_idx,urb_mean_idx,ag_idx,ag_mean_idx)

urb_combs <- all.pred %>% 
  mutate(tmp_idx = 2,tmp_mean_idx = 2,pre_idx = 2,pre_mean_idx = 2,ag_idx = 2,ag_mean_idx = 2) %>% 
  expand(tmp_idx,pre_idx,tmp_mean_idx,pre_mean_idx,urb_idx,urb_mean_idx,ag_idx,ag_mean_idx)

ag_combs <- all.pred %>% 
  mutate(tmp_idx = 2,tmp_mean_idx = 2,pre_idx = 2,pre_mean_idx = 2,urb_idx = 2,urb_mean_idx = 2) %>% 
  expand(tmp_idx,pre_idx,tmp_mean_idx,pre_mean_idx,urb_idx,urb_mean_idx,ag_idx,ag_mean_idx)

cov.combs <- bind_rows(clim_combs,urb_combs,ag_combs)

d2 <- cov.combs %>% 
  left_join(all.pred[,c("tmp","tmp_idx")],by=c("tmp_idx")) %>% 
  left_join(all.pred[,c("pre","pre_idx")],by=c("pre_idx")) %>% 
  left_join(all.pred[,c("urban","urb_idx")],by=c("urb_idx")) %>% 
  left_join(all.pred[,c("agriculture","ag_idx")],by=c("ag_idx")) %>% 
  left_join(all.pred[,c("tmp_mean","tmp_mean_idx")],by=c("tmp_mean_idx")) %>% 
  left_join(all.pred[,c("pre_mean","pre_mean_idx")],by=c("pre_mean_idx")) %>% 
  left_join(all.pred[,c("urb_mean","urb_mean_idx")],by=c("urb_mean_idx")) %>% 
  left_join(all.pred[,c("ag_mean","ag_mean_idx")],by=c("ag_mean_idx")) %>% 
  left_join(all.pred[,c("brood","volt_idx")],by=c("volt_idx")) %>% 
  mutate(log_y = NA,species=NA,strat_idx1 = NA, sp_idx1 = NA,sp_idx2 = NA,sp_idx3 = NA,sp_idx4 = NA,sp_idx5 = NA)
d3 <- d1 %>% 
  bind_rows(d2 %>% mutate(brood = as.integer(brood))  %>% dplyr::select(-c(tmp_idx,pre_idx,urb_idx,ag_idx,tmp_mean_idx,pre_mean_idx,urb_mean_idx,ag_mean_idx,volt_idx)))

pred_idx1 <- (nrow(d1)+1):nrow(d3)

#run top model
top.model <- inla(mod.list[[53]], data=d3,scale = log_prec_y,
                  family = "gaussian",
                  control.predictor=list(compute=T, link=1),
                  control.compute=list(config=T, waic = T, cpo=T,
                                       control.gcpo=list(enable=T)),
                  verbose=T, inla.mode="experimental",
                  control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
                  control.family=list(hyper=list(prec=list(initial=log(1), fixed=TRUE))))
top.model$summary.fixed
top.model$summary.fixed %>% 
  select(1,3,5) %>% 
  rename(lcl = `0.025quant`,ucl = `0.975quant`) %>% 
  mutate(diff_zero = ifelse(lcl > 0 | ucl < 0,1,0)) %>% 
  mutate_if(is.numeric,round,3) %>% 
  write.csv("final results/top model covariate coefficients.csv")
# extract predicted abundance indices from the prediction rows
# added above to d3. then multiply by stratum scalars A and z
d3$abun_fit <- top.model$summary.linear.predictor$mean
d3$abun_lcl <- top.model$summary.linear.predictor$`0.025quant`
d3$abun_ucl <- top.model$summary.linear.predictor$`0.975quant`

all.pred$brood = as.factor(all.pred$brood)
d3$brood = as.factor(d3$brood)

pred1 <- d3 %>% 
  slice(pred_idx1) %>% 
  left_join(all.pred[,c("tmp","tmp_idx")],by=c("tmp")) %>% 
  left_join(all.pred[,c("pre","pre_idx")],by=c("pre")) %>% 
  left_join(all.pred[,c("urban","urb_idx")],by=c("urban")) %>% 
  left_join(all.pred[,c("agriculture","ag_idx")],by=c("agriculture")) %>% 
  left_join(all.pred[,c("tmp_mean","tmp_mean_idx")],by=c("tmp_mean")) %>% 
  left_join(all.pred[,c("pre_mean","pre_mean_idx")],by=c("pre_mean")) %>% 
  left_join(all.pred[,c("urb_mean","urb_mean_idx")],by=c("urb_mean")) %>% 
  left_join(all.pred[,c("ag_mean","ag_mean_idx")],by=c("ag_mean")) %>% 
  left_join(all.pred[,c("brood","volt_idx")],by=c("brood"))

# write.csv(pred1,"data/predictions_top_model.csv")

#build covariate plot
pred2 <- pred1 %>% 
  mutate_if(is.numeric,round,2) %>% 
  dplyr::select(tmp,pre,urban,agriculture,tmp_mean,pre_mean,urb_mean,ag_mean,tmp_idx,pre_idx,urb_idx,ag_idx,tmp_mean_idx,pre_mean_idx,urb_mean_idx,ag_mean_idx,abun_fit,abun_lcl,abun_ucl,brood,volt_idx)
# mutate(brood2 = car::recode(brood,"'U' = 'Univoltine'; 'B' = 'Bivoltine'; 'M' = 'Multivoltine'"))

pred2$brood = factor(pred2$brood,levels = c("1","2","3"))
pred2 = pred2 %>% mutate(Brood = brood)

p.clim <- ggplot(pred2 %>% filter(pre_idx == 2 & urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & pre_mean_idx == 2 & tmp_mean_idx %in% c(1,3) & brood %in% c("1","2","3")),aes(x=tmp, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl,group=tmp_mean_idx,fill=tmp_mean)) +
  geom_point(inherit.aes = F,data=d3 %>% filter(!(is.na(brood)) & brood %in% c("1","2","3")),aes(y=abun_fit,x=tmp,col=tmp_mean),alpha=.01)+
  geom_ribbon(alpha=0.75,col="black") +
  # geom_line() +
  scale_color_gradient("Regional\nMAT (sd)",high="red",low="blue")+ 
  scale_fill_gradient(high="red",low="blue")+
  guides(fill="none") +
  labs(x="Annual Temperature (sd)", y="ln-Abundance (95% CI)") +
  facet_grid(~Brood,scales = "free",labeller = label_both)+
  theme_classic()+
  theme(axis.title = element_text(size = 14))
p.clim

p.clim2 <- ggplot(pred2 %>% filter(tmp_idx == 2 & urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & tmp_mean_idx == 2 & pre_mean_idx %in% c(1,3) & brood %in% c("1","2","3")),aes(x=pre, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl,group=pre_mean_idx,fill=pre_mean)) +
  geom_point(inherit.aes = F,data=d3 %>% filter(!(is.na(brood))& brood %in% c("1","2","3")),aes(y=abun_fit,x=pre,col=pre_mean),alpha=.01)+
  geom_ribbon(alpha=0.75,col="black") +
  # geom_line() +
  scale_color_gradient("Regional\nMAP (sd)",high="blue",low="red")+
  scale_fill_gradient(high="blue",low="red")+
  guides(fill="none") +
  labs(x="Annual Precipitation (sd)", y="ln-Abundance (95% CI)") +
  facet_grid(~Brood,scales = "free",labeller = label_both)+
  theme_classic()+
  theme(axis.title = element_text(size = 14))
p.clim2

#view temp by precip interaction
# ggplot(pred2 %>% filter(urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & tmp_mean_idx %in% c(1,3) & pre_mean_idx %in% c(1,3) & pre_idx %in% c(2) & Brood == "1"),aes(x=tmp, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl)) +
#   # geom_point(inherit.aes = F,data=d3%>% filter(!(is.na(brood)) & brood %in% c("1","2","3")),aes(y=abun_fit,x=pre,col=tmp),alpha=.01)+
#   # geom_ribbon(alpha=0.75,col="black") +
#   geom_line() +
#   # scale_color_gradient("Annual\nMAP (sd)",high="blue",low="red")+
#   # scale_fill_gradient(high="blue",low="red")+
#   guides(fill="none") +
#   labs(x="Annual Temperature (sd)", y="ln-Abundance (95% CI)") +
#   facet_grid(pre_mean_idx~tmp_mean_idx,scales = "free",labeller = label_both)+
#   theme_classic()+
#   theme(axis.title = element_text(size = 14))
# 
# ggplot(pred2 %>% filter(urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & tmp_mean_idx %in% c(1,3) & pre_mean_idx %in% c(1,3) & pre_idx %in% c(2) & Brood == "3"),aes(x=tmp, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl)) +
#   # geom_point(inherit.aes = F,data=d3%>% filter(!(is.na(brood)) & brood %in% c("1","2","3")),aes(y=abun_fit,x=pre,col=tmp),alpha=.01)+
#   # geom_ribbon(alpha=0.75,col="black") +
#   geom_line() +
#   # scale_color_gradient("Annual\nMAP (sd)",high="blue",low="red")+
#   # scale_fill_gradient(high="blue",low="red")+
#   guides(fill="none") +
#   labs(x="Annual Temperature (sd)", y="ln-Abundance (95% CI)") +
#   facet_grid(pre_mean_idx~tmp_mean_idx,scales = "free",labeller = label_both)+
#   theme_classic()+
#   theme(axis.title = element_text(size = 14))
# 
# ggplot(pred2 %>% filter(urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & tmp_mean_idx %in% c(1,3) & pre_mean_idx %in% c(1,3) & tmp_idx %in% c(2) & Brood == "1"),aes(x=pre, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl)) +
#   # geom_point(inherit.aes = F,data=d3%>% filter(!(is.na(brood)) & brood %in% c("1","2","3")),aes(y=abun_fit,x=pre,col=tmp),alpha=.01)+
#   # geom_ribbon(alpha=0.75,col="black") +
#   geom_line() +
#   # scale_color_gradient("Annual\nMAP (sd)",high="blue",low="red")+
#   # scale_fill_gradient(high="blue",low="red")+
#   guides(fill="none") +
#   labs(x="Annual Precipitation (sd)", y="ln-Abundance (95% CI)") +
#   facet_grid(pre_mean_idx~tmp_mean_idx,scales = "free",labeller = label_both)+
#   theme_classic()+
#   theme(axis.title = element_text(size = 14))
# 
# ggplot(pred2 %>% filter(urb_idx == 2 & ag_idx == 2 & urb_mean_idx == 2 & ag_mean_idx == 2 & tmp_mean_idx %in% c(1,3) & pre_mean_idx %in% c(1,3) & tmp_idx %in% c(2) & Brood == "3"),aes(x=pre, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl)) +
#   # geom_point(inherit.aes = F,data=d3%>% filter(!(is.na(brood)) & brood %in% c("1","2","3")),aes(y=abun_fit,x=pre,col=tmp),alpha=.01)+
#   # geom_ribbon(alpha=0.75,col="black") +
#   geom_line() +
#   # scale_color_gradient("Annual\nMAP (sd)",high="blue",low="red")+
#   # scale_fill_gradient(high="blue",low="red")+
#   guides(fill="none") +
#   labs(x="Annual Precipitation (sd)", y="ln-Abundance (95% CI)") +
#   facet_grid(pre_mean_idx~tmp_mean_idx,scales = "free",labeller = label_both)+
#   theme_classic()+
#   theme(axis.title = element_text(size = 14))

p.urb <- ggplot(pred2 %>% filter(pre_idx == 2 & tmp_idx == 2 & ag_idx == 2 & urb_mean_idx %in% c(1,3) & ag_mean_idx == 2 & pre_mean_idx == 2 & tmp_mean_idx == 2),aes(x=urban, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl,group=urb_mean_idx,fill=urb_mean)) +
  geom_point(inherit.aes = F,data=d3,aes(y=abun_fit,x=urban,col=urb_mean),alpha=.01)+
  geom_ribbon(alpha=0.75,col="black") +
  # geom_line() +
  scale_color_scico("Regional\nUrban (sd)",palette = "vik")+ 
  scale_fill_scico(palette = "vik")+
  guides(fill="none") +
  labs(x="Annual Urban (sd)", y="ln-Abundance (95% CI)") +
  # facet_grid(~tmp_mean,scales = "free",labeller = label_both)+
  theme_classic()+
  theme(axis.title = element_text(size = 14))
p.urb

p.ag <- ggplot(pred2 %>% filter(pre_idx == 2 & tmp_idx == 2 & urb_idx == 2 & urb_mean_idx == 2 & tmp_mean_idx == 2 & pre_mean_idx == 2 & ag_mean_idx %in% c(1,3)),aes(x=agriculture, y=abun_fit, ymin=abun_lcl,ymax=abun_ucl,group=ag_mean_idx,fill=ag_mean)) +
  geom_point(inherit.aes = F,data=d3,aes(y=abun_fit,x=agriculture,col=ag_mean),alpha=.01)+
  geom_ribbon(alpha=0.75,col="black") +
  # geom_line() +
  scale_color_viridis_c("Regional\nAgriculture (sd)")+ 
  scale_fill_viridis_c()+
  guides(fill="none") +
  labs(x="Annual Agriculture (sd)", y="ln-Abundance (95% CI)") +
  # facet_grid(~tmp_mean,scales = "free",labeller = label_both)+
  theme_classic()+
  theme(axis.title = element_text(size = 14))
p.ag

cowplot::plot_grid(p.clim,p.urb,labels = c("A","B"),nrow = 2)
ggsave("plots/temperature and urbanization.jpg",bg="white",height = 8,width = 6)

cowplot::plot_grid(p.clim2,p.ag,labels = c("A","B"),nrow = 2)
ggsave("plots/precip and ag effects.jpg",bg="white",height = 8,width = 6)

