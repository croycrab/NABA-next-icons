#climate and land use change model selection -- generate list of models to compare

# library(package)

# create list of relevant models ------------------------
mod.list <- list()
model_names <- c("intercept","climate1","climate2","urbanization","agriculture","climate1 + urbanization","climate2 + urbanization","climate2 + agriculture","climate2 + agriculture","urbanization + agriculture","climate1 + urbanization + agrictulture","climate2 + urbanization + agrictulture","climate1*urbanization","climate2*urbanization","climate1*agriculture","climate2*agriculture","urbanization*agriculture","climate1","climate2","urbanization","agriculture","climate1 + urbanization","climate2 + urbanization","climate2 + agriculture","climate2 + agriculture","urbanization + agriculture","climate1 + urbanization + agrictulture","climate2 + urbanization + agrictulture","climate1*urbanization","climate2*urbanization","climate1*agriculture","climate2*agriculture","urbanization*agriculture","brood*climate1","brood*climate2","brood*climate1 + urbanization","brood*climate2 + urbanization","brood*climate2 + agriculture","brood*climate2 + agriculture","brood*climate1 + urbanization + agrictulture","brood*climate2 + urbanization + agrictulture","brood*climate1*urbanization","brood*climate2*urbanization","brood*climate1*agriculture","brood*climate2*agriculture","brood*climate1","brood*climate2","brood*climate1 + urbanization","brood*climate2 + urbanization","brood*climate2 + agriculture","brood*climate2 + agriculture","brood*climate1 + urbanization + agrictulture","brood*climate2 + urbanization + agrictulture","brood*climate1*urbanization","brood*climate2*urbanization","brood*climate1*agriculture","brood*climate2*agriculture")



model_type <- c("intercept",rep("no_means/no_traits",16),rep("with_means/no_traits",16),rep("no_means/with_traits",12),rep("with_means/with_traits",12))

#intercept
mod.list[[1]]<- log_y ~ 1 +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate only
mod.list[[2]]<- log_y ~ tmp + pre +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate only 2 --temp*precip
mod.list[[3]]<- log_y ~ tmp*pre +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urban only
mod.list[[4]]<- log_y ~ urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#agriculture only
mod.list[[5]]<- log_y ~ agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx5,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and urbanization
mod.list[[6]]<- log_y ~ tmp + pre + urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization
mod.list[[7]]<- log_y ~ tmp*pre + urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and agriculture
mod.list[[8]]<- log_y ~ tmp + pre + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and agriculture
mod.list[[9]]<- log_y ~ tmp*pre + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urbanization and agriculture
mod.list[[10]]<- log_y ~ urban + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)

#climate and urbanization and agriculture
mod.list[[11]]<- log_y ~ tmp + pre + urban + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization and agriculture
mod.list[[12]]<- log_y ~ tmp*pre + urban + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*urbanization
mod.list[[13]]<- log_y ~ tmp*urban + pre*urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*urbanization
mod.list[[14]]<- log_y ~ tmp*pre*urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*agriculture
mod.list[[15]]<- log_y ~ tmp*agriculture + pre*agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*agriculture
mod.list[[16]]<- log_y ~ tmp*pre*agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urbanization*agriculture
mod.list[[17]]<- log_y ~ urban*agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate only (w/means)
mod.list[[18]]<- log_y ~ tmp*tmp_mean + pre*pre_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 only (w/means)
mod.list[[19]]<- log_y ~ tmp*tmp_mean*pre*pre_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urban only (w/means)
mod.list[[20]]<- log_y ~ urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#agriculture only (w/means)
mod.list[[21]]<- log_y ~ agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  f(sp_idx5,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and urbanization (w/means)
mod.list[[22]]<- log_y ~ tmp*tmp_mean + pre*pre_mean + urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization (w/means)
mod.list[[23]]<- log_y ~ tmp*pre*tmp_mean*pre_mean + urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and agriculture (w/means)
mod.list[[24]]<- log_y ~ tmp*tmp_mean + pre*pre_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and agriculture (w/means)
mod.list[[25]]<- log_y ~ tmp*pre*tmp_mean*pre_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urbanization and agriculture (w/means)
mod.list[[26]]<- log_y ~ urban*urb_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)

#climate and urbanization and agriculture (w/means)
mod.list[[27]]<- log_y ~ tmp*tmp_mean + pre*pre_mean + urban*urb_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization and agriculture (w/means)
mod.list[[28]]<- log_y ~ tmp*pre*tmp_mean*pre_mean + urban*urb_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)


#climate*urbanization (w/means)
mod.list[[29]]<- log_y ~ tmp*tmp_mean*urban*urb_mean + pre*pre_mean*urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*urbanization (w/means)
mod.list[[30]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*agriculture (w/means)
mod.list[[31]]<- log_y ~ tmp*tmp_mean*agriculture*ag_mean + pre*pre_mean*agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*agriculture (w/means)
mod.list[[32]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#urbanization*agriculture (w/means)
mod.list[[33]]<- log_y ~ urban*urb_mean*agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#with brood
#climate only
mod.list[[34]]<- log_y ~ tmp*brood + pre*brood + 
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate only 2 --temp*precip
mod.list[[35]]<- log_y ~ tmp*pre*brood +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)


#climate and urbanization
mod.list[[36]]<- log_y ~ tmp*brood + pre*brood + urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization
mod.list[[37]]<- log_y ~ tmp*pre*brood + urban +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and agriculture
mod.list[[38]]<- log_y ~ tmp*brood + pre*brood + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and agriculture
mod.list[[39]]<- log_y ~ tmp*pre*brood + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and urbanization and agriculture
mod.list[[40]]<- log_y ~ tmp*brood + pre*brood + urban + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization and agriculture
mod.list[[41]]<- log_y ~ tmp*pre*brood + urban + agriculture +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*urbanization
mod.list[[42]]<- log_y ~ tmp*urban + pre*urban + tmp*brood + pre*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*urbanization
mod.list[[43]]<- log_y ~ tmp*pre*urban + tmp*pre*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*agriculture
mod.list[[44]]<- log_y ~ tmp*agriculture + pre*agriculture + tmp*brood + pre*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*agriculture
mod.list[[45]]<- log_y ~ tmp*pre*agriculture + tmp*pre*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate only (w/means)
mod.list[[46]]<- log_y ~ tmp*tmp_mean*brood+ + pre*pre_mean*brood+ +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 only (w/means)
mod.list[[47]]<- log_y ~ tmp*tmp_mean*pre*pre_mean*brood+ +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+ #random intercept for strata
  #random species intercept
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+ #random slope of tmp within species
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+ #random slope of pre within species
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and urbanization (w/means)
mod.list[[48]]<- log_y ~ tmp*tmp_mean*brood+ + pre*pre_mean*brood+ + urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization (w/means)
mod.list[[49]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*brood+ + urban*urb_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and agriculture (w/means)
mod.list[[50]]<- log_y ~ tmp*tmp_mean*brood+ + pre*pre_mean*brood+ + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and agriculture (w/means)
mod.list[[51]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*brood+ + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate and urbanization and agriculture (w/means)
mod.list[[52]]<- log_y ~ tmp*tmp_mean*brood+ + pre*pre_mean*brood+ + urban*urb_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2 and urbanization and agriculture (w/means)
mod.list[[53]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*brood+ + urban*urb_mean + agriculture*ag_mean +
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)


#climate*urbanization (w/means)
mod.list[[54]]<- log_y ~ tmp*tmp_mean*urban*urb_mean + pre*pre_mean*urban*urb_mean + tmp*tmp_mean*brood + pre*pre_mean*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*urbanization (w/means)
mod.list[[55]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*urban*urb_mean + tmp*pre*tmp_mean*pre_mean*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx4,urban, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate*agriculture (w/means)
mod.list[[56]]<- log_y ~ tmp*tmp_mean*agriculture*ag_mean + pre*pre_mean*agriculture*ag_mean + tmp*tmp_mean*brood + pre*pre_mean*brood+ 
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

#climate2*agriculture (w/means)
mod.list[[57]]<- log_y ~ tmp*pre*tmp_mean*pre_mean*agriculture*ag_mean +tmp*pre*tmp_mean*pre_mean*brood+
  f(strat_idx1, model="besag", graph=g1, constr=T, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  
  f(sp_idx2,tmp, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx3,pre, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx5,agriculture, model="iid", constr=T, hyper=prec_prior)+
  f(sp_idx1, model = "generic0", Cmatrix = phylo_prec_mat.full)

save(mod.list,file = "model selection list.RData")
