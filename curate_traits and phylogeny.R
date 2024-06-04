#curate trait data
#combines various sources of trait data, and cleans up trait data by fixing species names, renaming traits, etc. 
#in hindsight, this is super messy, but the  main jist of it is that species level traits were cleaned up and appended to the dataframe containing species-level abundance trends (vs strata-level species abundance trends). I did this originally to plot the species-level trends, so this dataframe is used for that. There was no need to append traits for thsi however, but that's the way the cookie crumbles. I use the this global trend + trait dataframe later to append species level traits to strata-level trends, and then I append the strata-level voltinism data. I know, messy, but I don't have time to clean it all up. Hopefully it is clear however.


# global trends with traits -----------------------------------------------
#import trends 
df.but.global <- read.csv("data/trends/trends_1996-2018.csv") %>%
  rename('Species' = "species") %>%
  mutate(trend_prec = (1/((ci_range/3.92)^2)))dddd

#import LepTraits1.0 data and clean up
df.traits <-read.csv("data/traits/consensus.csv",header = T,na.strings=c(""," ","NA")) 
df.traits$Species[which(df.traits$Species == "Apodemia palmerii")] = "Apodemia palmeri"
df.traits$Species[which(df.traits$Species == "Glutophrissa drusilla")] = "Appias drusilla"
df.traits$Species[which(df.traits$Species == "Cecropterus casica")] = "Achalarus casica"
df.traits$Species[which(df.traits$Species == "Cecropterus lyciades")] = "Achalarus lyciades"
df.traits$Species[which(df.traits$Species == "Dione vanillae")] = "Agraulis vanillae"
df.traits$Species[which(df.traits$Species == "Apodemia palmerii")] = "Apodemia palmeri"
df.traits$Species[which(df.traits$Species == "Glutophrissa drusilla")] = "Appias drusilla"
df.traits$Species[which(df.traits$Species == "Brephidium exilis")] = "Brephidium exile"
df.traits$Species[which(df.traits$Species == "Zerene cesonia")] = "Colias cesonia"
df.traits$Species[which(df.traits$Species == "Zerene eurydice")] = "Colias eurydice"
df.traits$Species[which(df.traits$Species == "Apodemia zela")] = "Emesis zela"
df.traits$Species[which(df.traits$Species == "Lethe anthedon")] = "Enodia anthedon"
df.traits$Species[which(df.traits$Species == "Lethe creola")] = "Enodia creola"
df.traits$Species[which(df.traits$Species == "Lethe portlandia")] = "Enodia portlandia"
df.traits$Species[which(df.traits$Species == "Ephyriades brunnea")] = "Ephyriades brunneus"
df.traits$Species[which(df.traits$Species == "Gesta baptisiae")] = "Erynnis baptisiae"
df.traits$Species[which(df.traits$Species == "Gesta funeralis")] = "Erynnis funeralis"
df.traits$Species[which(df.traits$Species == "Gesta horatius")] = "Erynnis horatius"
df.traits$Species[which(df.traits$Species == "Gesta juvenalis")] = "Erynnis juvenalis"
df.traits$Species[which(df.traits$Species == "Gesta lucilius")] = "Erynnis lucilius"
df.traits$Species[which(df.traits$Species == "Gesta pacuvius")] = "Erynnis pacuvius"
df.traits$Species[which(df.traits$Species == "Gesta persius")] = "Erynnis persius"
df.traits$Species[which(df.traits$Species == "Gesta propertius")] = "Erynnis propertius"
df.traits$Species[which(df.traits$Species == "Gesta tristis")] = "Erynnis tristis"
df.traits$Species[which(df.traits$Species == "Gesta zarucco")] = "Erynnis zarucco"
df.traits$Species[which(df.traits$Species == "Abaeis mexicana")] = "Eurema mexicana"
df.traits$Species[which(df.traits$Species == "Pyrisitia lisa")] = "Eurema lisa"
df.traits$Species[which(df.traits$Species == "Abaeis nicippe")] = "Eurema nicippe"
df.traits$Species[which(df.traits$Species == "Pyrisitia nise")] = "Eurema nise"
df.traits$Species[which(df.traits$Species == "Pyrisitia proterpia")] = "Eurema proterpia"
df.traits$Species[which(df.traits$Species == "Cupido amyntula")] = "Everes amyntula"
df.traits$Species[which(df.traits$Species == "Cupido comyntas")] = "Everes comyntas"
df.traits$Species[which(df.traits$Species == "Echinargus isola")] = "Hemiargus isola"
df.traits$Species[which(df.traits$Species == "Junonia nigrosuffusa")] = "Junonia evarete"
df.traits$Species[which(df.traits$Species == "Plebejus idas")] = "Lycaeides idas"
df.traits$Species[which(df.traits$Species == "Plebejus melissa")] = "Lycaeides melissa"
df.traits$Species[which(df.traits$Species == "Cissia rubricata")] = "Megisto rubricata"
df.traits$Species[which(df.traits$Species == "Neonympha areolatus")] = "Neonympha areolata"
df.traits$Species[which(df.traits$Species == "Aglais milberti")] = "Nymphalis milberti"
df.traits$Species[which(df.traits$Species == "Nymphalis l-album")] = "Nymphalis vaualbum"
df.traits$Species[which(df.traits$Species == "Pterourus canadensis")] = "Papilio canadensis"
df.traits$Species[which(df.traits$Species == "Heraclides cresphontes")] = "Papilio cresphontes"
df.traits$Species[which(df.traits$Species == "Pterourus eurymedon")] = "Papilio eurymedon"
df.traits$Species[which(df.traits$Species == "Pterourus glaucus")] = "Papilio glaucus"
df.traits$Species[which(df.traits$Species == "Pterourus multicaudata")] = "Papilio multicaudata"
df.traits$Species[which(df.traits$Species == "Pterourus palamedes")] = "Papilio palamedes"
df.traits$Species[which(df.traits$Species == "Pterourus rutulus")] = "Papilio rutulus"
df.traits$Species[which(df.traits$Species == "Pterourus troilus")] = "Papilio troilus"
df.traits$Species[which(df.traits$Species == "Phyciodes pulchella")] = "Phyciodes campestris"
df.traits$Species[which(df.traits$Species == "Phyciodes cocyta")] = "Phyciodes selenis"
df.traits$Species[which(df.traits$Species == "Anthanassa texana")] = "Phyciodes texana"
df.traits$Species[which(df.traits$Species == "Phyciodes graphica")] = "Phyciodes vesta"
df.traits$Species[which(df.traits$Species == "Icaricia acmon")] = "Plebejus acmon"
df.traits$Species[which(df.traits$Species == "Icaricia icarioides")] = "Plebejus icarioides"
df.traits$Species[which(df.traits$Species == "Icaricia lupini")] = "Plebejus lupini"
df.traits$Species[which(df.traits$Species == "Icaricia saepiolus")] = "Plebejus saepiolus"
df.traits$Species[which(df.traits$Species == "Lon hobomok")] = "Poanes hobomok"
df.traits$Species[which(df.traits$Species == "Lon melane")] = "Poanes melane"
df.traits$Species[which(df.traits$Species == "Lon taxiles")] = "Poanes taxiles"
df.traits$Species[which(df.traits$Species == "Lon zabulon")] = "Poanes zabulon"
df.traits$Species[which(df.traits$Species == "Limochores mystic")] = "Polites mystic"
df.traits$Species[which(df.traits$Species == "Limochores origenes")] = "Polites origenes"
df.traits$Species[which(df.traits$Species == "Limochores sonora")] = "Polites sonora"
df.traits$Species[which(df.traits$Species == "Hedone vibex")] = "Polites vibex"
df.traits$Species[which(df.traits$Species == "Vernia verna")] = "Pompeius verna"
df.traits$Species[which(df.traits$Species == "Atrytone bulenta")] = "Problema bulenta"
df.traits$Species[which(df.traits$Species == "Atrytone byssus")] = "Problema byssus"
df.traits$Species[which(df.traits$Species == "Burnsius albescens")] = "Pyrgus albescens"
df.traits$Species[which(df.traits$Species == "Burnsius communis")] = "Pyrgus communis"
df.traits$Species[which(df.traits$Species == "Burnsius oileus")] = "Pyrgus oileus"
df.traits$Species[which(df.traits$Species == "Burnsius philetas")] = "Pyrgus philetas"
df.traits$Species[which(df.traits$Species == "Apyrrothrix araxes")] = "Pyrrhopyge araxes"
df.traits$Species[which(df.traits$Species == "Satyrium caryaevorus")] = "Satyrium caryaevorum"
df.traits$Species[which(df.traits$Species == "Lethe appalachia")] = "Satyrodes appalachia"
df.traits$Species[which(df.traits$Species == "Lethe appalachia")] = "Satyrodes appalachia"
df.traits$Species[which(df.traits$Species == "Lethe eurydice")] = "Satyrodes eurydice"
df.traits$Species[which(df.traits$Species == "Cecropterus bathyllus")] = "Thorybes bathyllus"
df.traits$Species[which(df.traits$Species == "Cecropterus confusis")] = "Thorybes confusis"
df.traits$Species[which(df.traits$Species == "Cecropterus drusius")] = "Thorybes drusius"
df.traits$Species[which(df.traits$Species == "Cecropterus mexicana")] = "Thorybes mexicanus"
df.traits$Species[which(df.traits$Species == "Cecropterus pylades")] = "Thorybes pylades"
df.traits$Species[which(df.traits$Species == "Cecropterus dorantes")] = "Urbanus dorantes"
df.traits$Species[which(df.traits$Species == "Brephidium pseudofea")] = "Brephidium isophthalma"
df.traits$Species[which(df.traits$Species == "Copaeodes aurantiaca")] = "Copaeodes aurantiacus"
df.traits$Species[which(df.traits$Species == "Copaeodes minima")] = "Copaeodes minimus"
df.traits$Species[which(df.traits$Species == "Hemiargus hanno")] = "Hemiargus ceraunus"
df.traits$Species[which(df.traits$Species == "Prolibythea carinenta")] = "Libytheana carinenta"
df.traits$Species[which(df.traits$Species == "Mestra dorcas")] = "Mestra amymone"
df.traits$Species[which(df.traits$Species == "Pheobis statira")] = "Phoebis statira"
df.traits$Species[which(df.traits$Species == "Argynnis aphrodite")] = "Speyeria aphrodite"
df.traits$Species[which(df.traits$Species == "Argynnis atlantis")] = "Speyeria atlantis"
df.traits$Species[which(df.traits$Species == "Argynnis callippe")] = "Speyeria callippe"
df.traits$Species[which(df.traits$Species == "Argynnis coronis")] = "Speyeria coronis"
df.traits$Species[which(df.traits$Species == "Argynnis cybele")] = "Speyeria cybele"
df.traits$Species[which(df.traits$Species == "Argynnis diana")] = "Speyeria diana"
df.traits$Species[which(df.traits$Species == "Argynnis edwardsii")] = "Speyeria edwardsii"
df.traits$Species[which(df.traits$Species == "Argynnis egleis")] = "Speyeria egleis"
df.traits$Species[which(df.traits$Species == "Argynnis hydaspe")] = "Speyeria hydaspe"
df.traits$Species[which(df.traits$Species == "Argynnis idalia")] = "Speyeria idalia"
df.traits$Species[which(df.traits$Species == "Argynnis mormonia")] = "Speyeria mormonia"
df.traits$Species[which(df.traits$Species == "Argynnis zerene")] = "Speyeria zerene"
df.traits$Species[which(df.traits$Species == "Thorybes mexicana")] = "Thorybes mexicanus"

df.traits = df.traits %>% 
  mutate(id = seq(1,length(Family),1)) %>% 
  filter(Species %in% df.but.global$Species) %>% 
  mutate(disturbance = car::recode(DisturbanceAffinity, "'Disturbance-associated (strong)'='Associated (strong)'; 'Disturbance-associated (weak)'='Associated (weak)';'Disturbance-avoidant (strong)'='Avoidant (strong)'; 'Disturbance-avoidant (weak)' = 'Avoidant (weak)'; 'Disturbance association varies' = 'Variable'; 'Seen near and away from disturbed habitat' = 'Variable'"),
         moisture = car::recode(MoistureAffinity, "'Xeric-associated (strong)'='Xeric (strong)'; 'Xeric-associated (weak)'='Xeric (weak)';'Mesic-associated (strong)'='Mesic (strong)'; 'Mesic-associated (weak)' = 'Mesic (weak)'; 'Moisture association varies' = 'Variable'; 'No evidence for moisture association' = NA; 'Both'='Variable'"),
         edge = car::recode(EdgeAffinity, "'Edge-associated (strong)'='Associated (strong)'; 'Edge-associated (weak)'='Associated (weak)';'Edge-avoidant (strong)'='Avoidant (strong)'; 'Edge-avoidant (weak)' = 'Avoidant (weak)'; 'Edge association varies' = 'Variable'; 'Seen near and away from edges' = 'Variable'; 'No evidence for edge associations' = NA"),
         canopy = car::recode(CanopyAffinity,"'Mixed canopy (open affinity)' = 'Open (weak)'; 'Canopy generalist' = 'Mixed'; 'Open canopy' = 'Open'; 'Mixed canopy (closed affinity)' = 'Closed (weak)'; 'Edge associated' = 'Mixed'; 'Mixed canopy' = 'Mixed'; 'Closed canopy' = 'Closed'"),
         disturbance2 = car::recode(DisturbanceAffinity, "'Disturbance-associated (strong)'='Associated'; 'Disturbance-associated (weak)'='Associated';'Disturbance-avoidant (strong)'='Avoidant'; 'Disturbance-avoidant (weak)' = 'Avoidant'; 'Disturbance association varies' = 'Variable'; 'Seen near and away from disturbed habitat' = 'Variable'"),
         moisture2 = car::recode(MoistureAffinity, "'Xeric-associated (strong)'='Xeric'; 'Xeric-associated (weak)'='Xeric';'Mesic-associated (strong)'='Mesic'; 'Mesic-associated (weak)' = 'Mesic'; 'Moisture association varies' = 'Variable'; 'No evidence for moisture association' = NA; 'Both'='Variable'"),
         diapause = car::recode(DiapauseStage,"'PL' = 'LP' ; 'PLL' = 'LP'; 'PPL' = 'LP'; 'LLEE' = 'EL'; 'PE' = 'EP'; 'PLE' = 'ELP'; 'PLEAA' = 'ELPA'; 'PPLAAA' = 'LPA'; 'LE' = 'EL'"),
         clutch = car::recode(OvipositionStyle,"'C, S' = 'SC';'S, C' = 'SC'; 'S,C' = 'SC'"),
         edge2 = car::recode(EdgeAffinity, "'Edge-associated (strong)'='Associated'; 'Edge-associated (weak)'='Associated';'Edge-avoidant (strong)'='Avoidant'; 'Edge-avoidant (weak)' = 'Avoidant'; 'Edge association varies' = 'Variable'; 'Seen near and away from edges' = 'Variable'; 'No evidence for edge associations' = NA"),
         canopy2 = car::recode(CanopyAffinity,"'Mixed canopy (open affinity)' = 'Mixed'; 'Canopy generalist' = 'Mixed'; 'Open canopy' = 'Open'; 'Mixed canopy (closed affinity)' = 'Mixed'; 'Edge associated' = 'Mixed'; 'Mixed canopy' = 'Mixed'; 'Closed canopy' = 'Closed'"))%>% 
  rename(diet_breadth = NumberOfHostplantFamilies,voltinism = Voltinism, flight_duration = FlightDuration)

#import Crossley trait data
df.traits2 <-openxlsx::read.xlsx("data/traits/NABA_traits_Crossley.xlsx",sheet = "no_header")
df.traits2$Species[which(df.traits2$Species == "Cupido amyntula")] = "Everes amyntula"
df.traits2 = df.traits2 %>%
  mutate(Species = paste(Genus,Species,sep=" ")) %>% 
  filter(Species %in% df.but.global$Species) %>% 
  rename(size = `Size.(mid-range)`,adult_aposematism = `cryptic,.aposematic`,larvae_aposematism = `aposematic/cryptic?`,hairs = `hairs/yes?`) %>%  
  mutate(across(where(is.character), str_trim),
         adult_aposematism = car::recode(adult_aposematism,"'crytpic' = 'cryptic';'cyptic' = 'cryptic'"),
         larvae_aposematism = car::recode(larvae_aposematism,"'aposemantic' = 'aposematic'"),
         hairs = car::recode(hairs,"'yes and yes' = 'yes'; 'both' = 'yes'"))

#import range size data
df.traits3 <- read.csv("data/traits/rangesize.csv")%>% 
  mutate(Species = str_to_sentence(species)) %>% 
  dplyr::select(-species,-X)

#import voltinism data
df.traits4 <- read.csv("data/traits/voltinism.csv") %>% 
  mutate(Species = str_to_sentence(species)) %>% 
  dplyr::select(-species,-X)

#append and clean up trait data durther
# hostplant = ifelse(is.na(SoleHostplantFamily),PrimaryHostplantFamily,SoleHostplantFamily) -- need to account for primary and secondary hosts being equal
id.rm <- c(11456, 11552 ,10808, 963, 1185, 8584, 9527, 5769)
trait1.vars<-c("disturbance","edge","moisture","canopy","disturbance2","edge2","moisture2","canopy2","diet_breadth","diapause","clutch","voltinism","flight_duration")
trait2.vars<-c("size","adult_aposematism","larvae_aposematism","hairs")

df.global<- df.but.global %>%
  left_join(df.traits[,c("Species",trait1.vars,"id")],by="Species") %>%
  filter(!(id %in% id.rm)) %>%
  left_join(df.traits2[,c("Species",trait2.vars,"Family")],by="Species")%>%
  left_join(df.traits3,by="Species")%>%
  # left_join(df.traits4,by="Species")%>%
  mutate_at(c("size","diet_breadth","flight_duration","range_size"),scale_this) %>% 
  # mutate(trend_prec = (1/((ci_range/3.92)^2)))%>%
  mutate_at("Family",str_to_title)

#NOTE: for LepTraits1.0. duplicate species = Colias alexandra, Erynnis funeralis, Hesperia comma, Apodemia mormo, Calephelis arizonensis (kept record with data vs empty record), Junonia evarete (kept record with same verbatim and species names), Polygonia oreas (kept record with multiple plant families)

# create tree with species where we have full trait data ------------------------
#this section is essentially getting a list of species for which we have trait data for, loading in the phylogenetic tree and cleaning up species names, subsetting the tree to only include species with data, then computing a distance matrix for inclusion in model

#view how many species we have with data for each trait 
complete.trait.df = df.global %>% 
  dplyr::select(disturbance,edge,canopy,moisture,diet_breadth,clutch,diapause,hairs,adult_aposematism,larvae_aposematism,range_size,flight_duration,size,voltinism)
#get sample size for each trait
trait.count = data.frame(trait = colnames(complete.trait.df)[1:14])
for(i in 1:14){
  trait.count$sample_size[i] = sum(!is.na(complete.trait.df[,trait.count$trait[i]]))
}

#vector of species with complete trait data
sp.all.traits = read.csv("data/trends/trends by strata_1996-2018.csv") %>% 
  mutate(trend_prec = (1/((ci_range/3.92)^2))) %>% 
  left_join(df.global[,-c(1:6,8)],by=c("species" = "Species"))%>%
  left_join(df.traits4[,c("Species","strat_name","brood")],by=c("species" = "Species","strat_name")) %>% 
  dplyr::select(mean_trend,disturbance,edge,canopy,moisture,diet_breadth,clutch,diapause,hairs,brood,larvae_aposematism,adult_aposematism,range_size,flight_duration,size,species,strat_name,trend_prec) %>%
  na.omit() %>% 
  pull(species) %>% 
  unique()

tree = read.tree("data/phylogenetic tree/US_bugs_calibrated.tre")
tree$tip.label = str_split(tree$tip.label,"\\|",simplify = T)[,1]
tree$tip.label = str_split(tree$tip.label,"\\'",simplify = T)[,2]
#fix species names
tree$tip.label[which(tree$tip.label == "Cecropterus casica")] = "Achalarus casica"
tree$tip.label[which(tree$tip.label == "Cecropterus lyciades")] = "Achalarus lyciades"
tree$tip.label[which(tree$tip.label == "Dione vanillae")] = "Agraulis vanillae"
tree$tip.label[which(tree$tip.label == "Apodemia palmerii")] = "Apodemia palmeri"
tree$tip.label[which(tree$tip.label == "Glutophrissa drusilla")] = "Appias drusilla"
tree$tip.label[which(tree$tip.label == "Brephidium exilis")] = "Brephidium exile"
tree$tip.label[which(tree$tip.label == "Zerene cesonia")] = "Colias cesonia"
tree$tip.label[which(tree$tip.label == "Zerene eurydice")] = "Colias eurydice"
tree$tip.label[which(tree$tip.label == "Apodemia zela")] = "Emesis zela"
tree$tip.label[which(tree$tip.label == "Lethe anthedon")] = "Enodia anthedon"
tree$tip.label[which(tree$tip.label == "Lethe creola")] = "Enodia creola"
tree$tip.label[which(tree$tip.label == "Lethe portlandia")] = "Enodia portlandia"
tree$tip.label[which(tree$tip.label == "Ephyriades brunnea")] = "Ephyriades brunneus"
tree$tip.label[which(tree$tip.label == "Gesta baptisiae")] = "Erynnis baptisiae"
tree$tip.label[which(tree$tip.label == "Gesta funeralis")] = "Erynnis funeralis"
tree$tip.label[which(tree$tip.label == "Gesta horatius")] = "Erynnis horatius"
tree$tip.label[which(tree$tip.label == "Gesta juvenalis")] = "Erynnis juvenalis"
tree$tip.label[which(tree$tip.label == "Gesta lucilius")] = "Erynnis lucilius"
tree$tip.label[which(tree$tip.label == "Gesta pacuvius")] = "Erynnis pacuvius"
tree$tip.label[which(tree$tip.label == "Gesta persius")] = "Erynnis persius"
tree$tip.label[which(tree$tip.label == "Gesta propertius")] = "Erynnis propertius"
tree$tip.label[which(tree$tip.label == "Gesta tristis")] = "Erynnis tristis"
tree$tip.label[which(tree$tip.label == "Gesta zarucco")] = "Erynnis zarucco"
tree$tip.label[which(tree$tip.label == "Abaeis mexicana")] = "Eurema mexicana"
tree$tip.label[which(tree$tip.label == "Pyrisitia lisa")] = "Eurema lisa"
tree$tip.label[which(tree$tip.label == "Abaeis nicippe")] = "Eurema nicippe"
tree$tip.label[which(tree$tip.label == "Pyrisitia nise")] = "Eurema nise"
tree$tip.label[which(tree$tip.label == "Pyrisitia proterpia")] = "Eurema proterpia"
tree$tip.label[which(tree$tip.label == "Cupido amyntula")] = "Everes amyntula"
tree$tip.label[which(tree$tip.label == "Cupido comyntas")] = "Everes comyntas"
tree$tip.label[which(tree$tip.label == "Echinargus isola")] = "Hemiargus isola"
tree$tip.label[which(tree$tip.label == "Junonia nigrosuffusa")] = "Junonia evarete"
tree$tip.label[which(tree$tip.label == "Plebejus idas")] = "Lycaeides idas"
tree$tip.label[which(tree$tip.label == "Plebejus melissa")] = "Lycaeides melissa"
tree$tip.label[which(tree$tip.label == "Cissia rubricata")] = "Megisto rubricata"
tree$tip.label[which(tree$tip.label == "Neonympha areolatus")] = "Neonympha areolata"
tree$tip.label[which(tree$tip.label == "Aglais milberti")] = "Nymphalis milberti"
tree$tip.label[which(tree$tip.label == "Nymphalis l-album")] = "Nymphalis vaualbum"
tree$tip.label[which(tree$tip.label == "Pterourus canadensis")] = "Papilio canadensis"
tree$tip.label[which(tree$tip.label == "Heraclides cresphontes")] = "Papilio cresphontes"
tree$tip.label[which(tree$tip.label == "Pterourus eurymedon")] = "Papilio eurymedon"
tree$tip.label[which(tree$tip.label == "Pterourus glaucus")] = "Papilio glaucus"
tree$tip.label[which(tree$tip.label == "Pterourus multicaudata")] = "Papilio multicaudata"
tree$tip.label[which(tree$tip.label == "Pterourus palamedes")] = "Papilio palamedes"
tree$tip.label[which(tree$tip.label == "Pterourus rutulus")] = "Papilio rutulus"
tree$tip.label[which(tree$tip.label == "Pterourus troilus")] = "Papilio troilus"
tree$tip.label[which(tree$tip.label == "Phyciodes pulchella")] = "Phyciodes campestris"
tree$tip.label[which(tree$tip.label == "Phyciodes cocyta")] = "Phyciodes selenis"
tree$tip.label[which(tree$tip.label == "Anthanassa texana")] = "Phyciodes texana"
tree$tip.label[which(tree$tip.label == "Phyciodes graphica")] = "Phyciodes vesta"
tree$tip.label[which(tree$tip.label == "Icaricia acmon")] = "Plebejus acmon"
tree$tip.label[which(tree$tip.label == "Icaricia icarioides")] = "Plebejus icarioides"
tree$tip.label[which(tree$tip.label == "Icaricia lupini")] = "Plebejus lupini"
tree$tip.label[which(tree$tip.label == "Icaricia saepiolus")] = "Plebejus saepiolus"
tree$tip.label[which(tree$tip.label == "Lon hobomok")] = "Poanes hobomok"
tree$tip.label[which(tree$tip.label == "Lon melane")] = "Poanes melane"
tree$tip.label[which(tree$tip.label == "Lon taxiles")] = "Poanes taxiles"
tree$tip.label[which(tree$tip.label == "Lon zabulon")] = "Poanes zabulon"
tree$tip.label[which(tree$tip.label == "Limochores mystic")] = "Polites mystic"
tree$tip.label[which(tree$tip.label == "Limochores origenes")] = "Polites origenes"
tree$tip.label[which(tree$tip.label == "Limochores sonora")] = "Polites sonora"
tree$tip.label[which(tree$tip.label == "Hedone vibex")] = "Polites vibex"
tree$tip.label[which(tree$tip.label == "Vernia verna")] = "Pompeius verna"
tree$tip.label[which(tree$tip.label == "Atrytone bulenta")] = "Problema bulenta"
tree$tip.label[which(tree$tip.label == "Atrytone byssus")] = "Problema byssus"
tree$tip.label[which(tree$tip.label == "Burnsius albescens")] = "Pyrgus albescens"
tree$tip.label[which(tree$tip.label == "Burnsius communis")] = "Pyrgus communis"
tree$tip.label[which(tree$tip.label == "Burnsius oileus")] = "Pyrgus oileus"
tree$tip.label[which(tree$tip.label == "Burnsius philetas")] = "Pyrgus philetas"
tree$tip.label[which(tree$tip.label == "Apyrrothrix araxes")] = "Pyrrhopyge araxes"
tree$tip.label[which(tree$tip.label == "Satyrium caryaevorus")] = "Satyrium caryaevorum"
tree$tip.label[which(tree$tip.label == "Lethe appalachia")] = "Satyrodes appalachia"
tree$tip.label[which(tree$tip.label == "Lethe appalachia")] = "Satyrodes appalachia"
tree$tip.label[which(tree$tip.label == "Lethe eurydice")] = "Satyrodes eurydice"
tree$tip.label[which(tree$tip.label == "Cecropterus bathyllus")] = "Thorybes bathyllus"
tree$tip.label[which(tree$tip.label == "Cecropterus confusis")] = "Thorybes confusis"
tree$tip.label[which(tree$tip.label == "Cecropterus drusius")] = "Thorybes drusius"
tree$tip.label[which(tree$tip.label == "Cecropterus mexicana")] = "Thorybes mexicanus"
tree$tip.label[which(tree$tip.label == "Cecropterus pylades")] = "Thorybes pylades"
tree$tip.label[which(tree$tip.label == "Cecropterus dorantes")] = "Urbanus dorantes"
tree$tip.label[which(tree$tip.label == "Aphrissa statira")] = "Phoebis statira"
tree$tip.label[which(tree$tip.label == "Telegonus cellus")] = "Autochton cellus"
tree$tip.label[which(tree$tip.label == "Brephidium pseudofea")] = "Brephidium isophthalma"

#full tree--all species
tree.full  = force.ultrametric(tree,method = "extend")
phylo_prec_mat.full = inverseA(tree.full,nodes = "TIPS",scale = TRUE)$Ainv
sp.order = data.frame(Species = phylo_prec_mat.full@Dimnames[[1]],sp_idx1 = 1:length(phylo_prec_mat.full@Dimnames[[1]]))

#trim tree -- complete trait data
tips.rm1 = tree$tip.label[which(!(tree$tip.label %in% sp.all.traits))]
tree.sub1 = tree %>% drop.tip(tips.rm1) %>% as.phylo()
# tree.sub1$node.label = 1:length(tree.sub1$node.label)
tree.sub1  = force.ultrametric(tree.sub1,method = "extend")
phylo_prec_mat.sub = inverseA(tree.sub1,nodes = "TIPS",scale = TRUE)$Ainv
sp.order.sub = data.frame(Species = phylo_prec_mat.sub@Dimnames[[1]],sp_idx1 = 1:length(phylo_prec_mat.sub@Dimnames[[1]]))
plot(tree.sub1)

# strata-level data prep --------------------------------------------------
df.strata <- read.csv("data/trends/trends by strata_1996-2018.csv") %>% 
  mutate(trend_prec = (1/((ci_range/3.92)^2))) %>% 
  left_join(df.global[,-c(1:6,8)],by=c("species" = "Species"))%>%
  left_join(df.traits4[,c("Species","strat_name","brood")],by=c("species" = "Species","strat_name"))  %>% 
  dplyr::select(mean_trend,disturbance,edge,canopy,moisture,diet_breadth,brood,hairs,clutch,diapause,adult_aposematism,larvae_aposematism,range_size,flight_duration,size,species,strat_name,trend_prec) %>%
  left_join(sp.order.sub,by = c("species" = "Species")) %>% 
  left_join(map2 %>% st_drop_geometry() %>% dplyr::select(strat_name,strat_idx1),by = "strat_name") %>%
  na.omit()

# save output for modeling ------------------------------------------------
save(df.global,df.strata,phylo_prec_mat.sub,phylo_prec_mat.full,tree.full, file = "traits and phylo.RData")

