########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Calculates the size structure of old-growth and mature secondary forest, and woody residence times between both forest age
## Controls for the size class distribution
## Produces: Figures S5, 3
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggpubr); library(ggplot2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12)) #, face ="bold"

#census-interval data
load("survival_data_assembled.Rdata") 
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                       simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                       !is.na(simple_combined$dbhmm.1.3),]

simple_combined_sp = simple_combined_sp[simple_combined_sp$openyear == 2005 & simple_combined_sp$endyear == 2010 |
                                          simple_combined_sp$openyear == 2003 & simple_combined_sp$endyear == 2008 |
                                          simple_combined_sp$openyear == 2004 & simple_combined_sp$endyear == 2014 |
                                          simple_combined_sp$openyear == 2000 & simple_combined_sp$endyear == 2014,]

#secondary forest data #1: collected between 2011-2014 by Dent et al. (https://datastorre.stir.ac.uk/handle/11667/144)
secondary = read.table("BCNM_secondary_forest_1ha_plots.txt", h = T) #secondary forest data #collected between 2011-2014
baseline = as.data.frame(secondary[secondary$site_age >=90,]) #oldest of the 2dary
baseline$quadrat_plot = baseline$site #Dent data has no finer spatial subdivision than the 1ha site itself, so site is the cluster unit for the bootstrap below

#secondary forest data #2: the Lang-Knight plot (1.5ha) was 90 years old in 1968, then I use the data from 1968
lang.k = read.table("LangKnight_1968to2014_treefile_fixed.txt", h = T)
lang.k = lang.k[!is.na(lang.k$dbh1968),]
lang.k$site = "LAN"; lang.k$site_age = 90

#add secondary forest data from 25ha, AVA (6ha), 10ha; includes all trees because we only have from 2004 - 2014 census
plot25ha = simple_combined_sp[simple_combined_sp$plot == "plot25ha" & simple_combined_sp$mortality == 0,] #5605 trees
plot25ha$site = "P25"; plot25ha$sp_code = plot25ha$sp; plot25ha$dbh = plot25ha$dbhmm.1.3; plot25ha$site_age = 140; #site ages in plot25ha are 144-100 years old
plot25ha = plot25ha[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plotava = simple_combined_sp[simple_combined_sp$plot == "plotAVA" & simple_combined_sp$mortality == 0,] #596 trees
plotava$site = "AVA"; plotava$sp_code = plotava$sp; plotava$dbh = plotava$dbhmm.1.3; plotava$site_age = 140; #site ages in plot25ha are 144-100 years old, mostly 100
plotava = plotava[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

#the 10ha plot has one part that is mature forest (remove it and use it below)
plot10ha.temp = simple_combined_sp[simple_combined_sp$plot == "plot10ha" & simple_combined_sp$mortality == 0,] #5605 trees

#around 420 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 420m there is mature forest, to the right of 420m it is secondary.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 780m there is young forest, to the right of 780m it has both secondary/mature.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and above 50m in the y-axis it is young forest. 
#quadrat coming from Maekem/Cushman
load("all2.rdata")
plot10ha.1 = plot10ha.temp[plot10ha.temp$gx>250 & plot10ha.temp$gx<750,] #young 10ha area 1 #5ha
plot10ha.2 = plot10ha.temp[plot10ha.temp$gx>800 & plot10ha.temp$gy>50,] #young 10ha area 2 #1ha
plot10ha.old = plot10ha.temp[plot10ha.temp$gx<200,] #old 10haplot ~4ha
plot10ha = rbind(plot10ha.1, plot10ha.2)

plot10ha$site = "P10"; plot10ha$sp_code = plot10ha$sp; plot10ha$dbh = plot10ha$dbhmm.1.3; plot10ha$site_age = 140; #site ages in plot25ha are 144-100 years old, mostly 100
plot10ha = plot10ha[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plot10ha.old$site = "P10"; plot10ha.old$sp_code = plot10ha.old$sp; plot10ha.old$dbh = plot10ha.old$dbhmm.1.3; plot10ha.old$site_age = 140; #site ages in plot25ha are 144-100 years old, mostly 100
plot10ha.old = plot10ha.old[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")] #2ha

#bind the 5 datasets: Dent, Lang-Knight- Plot25ha, PlotAVA, Plot10ha
names(baseline); names(lang.k)
lang.k$quadrat_plot = paste(lang.k$q20, "LAN", sep = "_") #q20 is the 20x20m quadrat id in the Lang-Knight treefile
lang.k = lang.k[,c("site", "site_age", "sp6", "dbh1968", "quadrat_plot")]
colnames(lang.k) = c("site", "site_age", "sp_code", "dbh", "quadrat_plot")
baseline = rbind(baseline, lang.k, plot25ha, plotava, plot10ha)
baseline$sp_code = tolower(baseline$sp_code)

#retrieve the species to use based on all of our analysis for the paper
load("ML_mortality_50CI.Rdata")
data = ML_mortality
species.to.use = tolower(unique(data$sp))

## young forest and mixed vicinity forest are total = 2.36 ha (with 1.92ha of young forest and 0.44ha of mixed forest)
# 50-ha plot
##where is the young forest areas in the 50ha plot?
#loads the topographic types for each quadrat from Harms et al., 2001 J.of Ecology
harms = read.table("harms.txt",h=T)
harms$quadrat = paste(harms$GX20/20, (harms$GY20/20), sep = "") #code to name the quadrats on the 50haplot
young.forest.50ha = harms[harms$Habitat == "young",] #48 quadrats = 1.92
mixed.forest.50ha = c("2123","2222","2422","2522","2622","2722","3320","3520","3620","3722","3824") #this is mixed in vicinity of the young forest #11quads 0.44ha

#Add young forest data to the secondary forest data from above (young 50ha is 1.92ha = 48quadrats*400m2)
plot50ha = simple_combined_sp[simple_combined_sp$plot == "plot50ha",]; nrow(plot50ha) #177964
plot50ha.young = plot50ha[which(plot50ha$quadrat %in% young.forest.50ha$quadrat),];nrow(plot50ha.young) #6887
plot50ha.old = plot50ha[-which(plot50ha$quadrat %in% young.forest.50ha$quadrat),]; nrow(plot50ha.old) #171077
plot50ha.old = plot50ha.old[-which(plot50ha.old$quadrat %in% mixed.forest.50ha),]; nrow(plot50ha.old) #169457
plot50ha.young$site = "P50"
plot50ha.young$site_age = 120
plot50ha.young = plot50ha.young[,c("site", "site_age", "sp", "dbhmm.1.3", "quadrat_plot")]
colnames(plot50ha.young) = c("site", "site_age", "sp_code", "dbh", "quadrat_plot")

##mature secondary
baseline = rbind(baseline, plot50ha.young)

## Add taxon to the data
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk_taxon = function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "taxon"]
  return(tx)
}
baseline$taxon = unlist(lapply(baseline$sp, funk_taxon))

funk_family = function(species){ 
  use = match(toString(species), bci.spptable$sp)
  f = bci.spptable[use, "Family"]
  return(f)
}
baseline$family = unlist(lapply(baseline$sp, funk_family))

#adding wood density to the dataframe for later biomass calculations. 
jw_data = read.table("WD_species.txt", h = T)
jw_data = data.frame(jw_data)
jw_data$SP. = tolower(jw_data$sp6)

funk_wd = function(species){
  use = match(toString(species), jw_data$SP.)
  jw_data[use, "wd100.mean"]
}
baseline$wd = unlist(lapply(baseline$sp, funk_wd))
baseline = as.data.frame(baseline)

#extras wood from ctfs data when the species did not have reported wood density in the previous dataset
load("wsg.ctfs.Rdata")
ctfs = wsg.ctfs3
ctfs$taxon = paste(ctfs$genus, ctfs$species, sep=" ")
with.nas = baseline[is.na(baseline$wd),]
unique.with.nas = unique(with.nas$taxon)
ctfs[which(ctfs$taxon == "Acacia melanoceras"), "taxon"] =  "Vachellia melanoceras"

funk_ctfs = function(tax){ 
  use = match(toString(tax), ctfs$taxon)
  wd2 = ctfs[use,"wsg"]
  baseline[which(baseline$taxon == tax), "wd"] <<- wd2
}
q = unlist(lapply(unique.with.nas, funk_ctfs))

#for species with No wood density value: the mean of all other species in the site was used
with.nas = baseline[is.na(baseline$wd),]
unique.with.nas = unique(with.nas$sp)
baseline[which(baseline$sp %in% unique.with.nas), "wd"] = mean(na.omit(baseline$wd))

#calculate their AGB
#when height is unknown - Chave et al., 2015 (Equation 7) for moist tropical forest
#D in "cm"; wsg in "g/cm3"
#Kg/1000 to convert to Mg
#dbhmm.1.3 is already corrected with taper equations when the data was cleaned

#For palm species -> ln(ABG) = -3.3488 + 2.7484*ln(dbh) * CF (in Kilograms)
#CF correction factor by log-transforming the data see Chave et al., 2005 (CF = exp((RSE^2)/2) = 1.188714)
palms = baseline[baseline$family == "Arecaceae",]
palms$AGB = (exp(-3.3488 + 2.7484*log(palms$dbh/10))/1000) * (1.188714)

trees = baseline[which(!baseline$family == "Arecaceae"),]
E = 0.051763982 #for BCI
trees$AGB = (exp(-1.803-0.976*E+(0.976*log(trees$wd))+
                   (2.673*log(trees$dbh/10))-
                   (0.0299*(log(trees$dbh/10))^2))/1000)
baseline = rbind(trees, palms)

#Basal Area
baseline$ba = (pi/4)*(baseline$dbh)^2
#Convert to m2 instead of mm2
baseline$ba = baseline$ba/1000000 

##old-growth

#where is the mixed forest in the vicinity of the above mixed forest?
#remove the quadrats that have young forest
plot50ha.old$site = "P50"
plot50ha.old$site_age = 500
plot50ha.old = plot50ha.old[,c("site", "site_age", "sp", "dbhmm.1.3", "quadrat_plot")]
colnames(plot50ha.old) = c("site", "site_age", "sp_code", "dbh", "quadrat_plot")

#add data from Pearson, Zetek and Drayton as mature forests
plotpearson= simple_combined_sp[simple_combined_sp$plot == "plotPearson" & simple_combined_sp$mortality == 0,] #5605 trees
plotpearson = plotpearson[plotpearson$endyear == 2014,] #select the same census as 25ha, AVA (2004-2014)
plotpearson$site = "PSN"; plotpearson$sp_code = plotpearson$sp; plotpearson$dbh = plotpearson$dbhmm.1.3; plotpearson$site_age = 500; #this is a random number
plotpearson = plotpearson[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plotzetek = simple_combined_sp[simple_combined_sp$plot == "plotZetek" & simple_combined_sp$mortality == 0,] #596 trees
plotzetek$site = "ZET"; plotzetek$sp_code = plotzetek$sp; plotzetek$dbh = plotzetek$dbhmm.1.3; plotzetek$site_age = 500; #site ages in plot25ha are 144-100 years old, mostly 100
plotzetek = plotzetek[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plotdrayton = simple_combined_sp[simple_combined_sp$plot == "plotDrayton" & simple_combined_sp$mortality == 0,] #596 trees
plotdrayton$site = "DRA"; plotdrayton$sp_code = plotdrayton$sp; plotdrayton$dbh = plotdrayton$dbhmm.1.3; plotdrayton$site_age = 500; #site ages in plot25ha are 144-100 years old, mostly 100
plotdrayton = plotdrayton[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plotGigante = simple_combined_sp[simple_combined_sp$plot == "Gigante" & simple_combined_sp$mortality == 0,] #12964 trees
plotGigante$site = "GIG"; plotGigante$sp_code = plotGigante$sp; plotGigante$dbh = plotGigante$dbhmm.1.3; plotGigante$site_age = 600; #site ages in plot25ha are 144-100 years old, mostly 100
plotGigante = plotGigante[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

#bind the with 50ha datase datasets
mature = rbind(plot50ha.old, plotpearson, plotzetek, plotdrayton, plot10ha.old, plotGigante)

## Add taxon to the data
mature$taxon = unlist(lapply(mature$sp_code, funk_taxon))
mature$family = unlist(lapply(mature$sp_code, funk_family))

#adding wood density to the dataframe for later biomass calculations
mature$wd = unlist(lapply(mature$sp_code, funk_wd))
mature = as.data.frame(mature)

#extras wd from ctfs data 
q = unlist(lapply(unique.with.nas, funk_ctfs))

#for species with No wood density value: the mean of all other species in the site was used
with.nas = mature[is.na(mature$wd),]
unique.with.nas = unique(with.nas$sp)

mature[which(mature$sp %in% unique.with.nas), "wd"] =
  mean(na.omit(mature$wd))

#calculate their AGB
#when height is unknown - Chave et al., 2015 (Equation 7) for moist tropical forest
#D in "cm"; wsg in "g/cm3"
#Kg/1000 to convert to Mg
#dbhmm.1.3 is already corrected with taper equations when the data was cleaned

#For palm species -> ln(ABG) = -3.3488 + 2.7484*ln(dbh) * CF (in Kilograms)
#CF correction factor by log-transforming the data see Chave et al., 2005 (CF = exp((RSE^2)/2) = 1.188714)
palms_mature = mature[mature$family == "Arecaceae",]
palms_mature$AGB = (exp(-3.3488 + 2.7484*log(palms_mature$dbh/10))/1000) * (1.188714)

trees_mature = mature[which(!mature$family == "Arecaceae"),]
E = 0.051763982 #for BCI
trees_mature$AGB = (exp(-1.803-0.976*E+(0.976*log(trees_mature$wd))+
                          (2.673*log(trees_mature$dbh/10))-
                          (0.0299*(log(trees_mature$dbh/10))^2))/1000)

mature = rbind(trees_mature, palms_mature)
#Basal Area
mature$ba = (pi/4)*(mature$dbh)^2
#Convert to m2 instead of mm2
mature$ba = mature$ba/1000000 

########################################################################################################################
#calculate crown area using species allometries
#import species allometry dataset
#Martínez- Cano et al. 2019 published the allometry data
species_allometry = read.csv("species_allometry.csv") 
#we do not have allometry for some species so we use the community-level parameters from Martinez-Cano table2.
#No allometries: Hura crepitans 
species_allometry[nrow(species_allometry) + 1,] = c("Hura crepitans", 57,54.5,60, 0.735,0.718,0.752, 21.77,20.70,22.89, 0.56,0.50,0.63, 1.35,1.31,1.38)
#No allometries: "vatear"
species_allometry[nrow(species_allometry)+1,] = c("Vatairea erythrocarpa", 57,54.5,60, 0.735,0.718,0.752, 21.77,20.7,22.89, 0.56,0.50,0.63, 1.35,1.31,1.38, "vataer")

#as.numeric columns 
species_allometry[,2:ncol(species_allometry)] = sapply(species_allometry[,2:ncol(species_allometry)],as.numeric)

##add mnemonic to the allometry data
funk= function(x){ 
  use = match(toString(x), bci.spptable$Latin)
  sp_code = bci.spptable[use,"sp"]
  return(sp_code)
}
species_allometry$sp_code = unlist(lapply(species_allometry$Species, funk))

#Aspidosperma spruceanum is synonym of our Aspidosperma cruentum
species_allometry[species_allometry$Species == "Aspidosperma spruceanum", "sp_code"] = "aspicr"
#Beilschmiedia pendula is synonym of our Beilschmiedia tovarensis
species_allometry[species_allometry$Species == "Beilschmiedia pendula", "sp_code"] = "beilpe"
#Guatteria dumetorum is synonym of our Guatteria lucens
species_allometry[species_allometry$Species == "Guatteria dumetorum", "sp_code"] = "guatdu"
#Quararibea asterolepis does not have the mnemonic
species_allometry[species_allometry$Species == "Quararibea asterolepis", "sp_code"] = "quaras"
#Tabebuia guayacan is now Handroanthus
species_allometry[species_allometry$Species == "Tabebuia guayacan", "sp_code"] = "tab1gu"
#Tachigali versicolor is synonym of our Tachigali panamensis
species_allometry[species_allometry$Species == "Tachigali versicolor", "sp_code"] = "tachve"
#Tetragastris panamensis is synonym of Protium stevensonii
species_allometry[species_allometry$Species == "Tetragastris panamensis", "sp_code"] = "tet2pa"
#Virola surinamensis is synonym of Virola nobilis
species_allometry[species_allometry$Species == "Virola surinamensis", "sp_code"] = "virosu"

species_allometry[, c(2:16)] = lapply(species_allometry[, c(2:16)], as.numeric)

#add allometric information to the secondary forest 
young_allometry<-left_join(baseline,species_allometry, by="sp_code")
young_allometry = as.data.frame(young_allometry)

#add allometric information to the mature forest
old_allometry<-left_join(mature,species_allometry, by="sp_code")
#View(old_allometry)
old_allometry = as.data.frame(old_allometry)

#trim datasets to use only our species
young_allometry = young_allometry[young_allometry$sp_code %in% species.to.use,]
old_allometry = old_allometry[old_allometry$sp_code %in% species.to.use,]

#find the species without allometric data
young_allometry[is.na(young_allometry$Tree.height.a),]$sp_code #all the 39 species have allometric data
old_allometry[is.na(old_allometry$Tree.height.a),]$sp_code #all the 39 species have allometric data

## calculate crown area:
#crown area (dbh in cm)
young_allometry$allometric_crownarea = young_allometry$Crown.area.a*(young_allometry$dbh/10)^young_allometry$Crown.area.b
old_allometry$allometric_crownarea = old_allometry$Crown.area.a*(old_allometry$dbh/10)^old_allometry$Crown.area.b

# #species for which we have young forest data from our 38 species in the study
# species.to.use.young = unique(young_allometry$sp_code) #38 present

# ##graphing for ONLY giant trees (500 mm)
young = young_allometry[young_allometry$dbh >= 500,]
young = young[!young$sp_code== "vataer",] #remove "vataer" as there are not individuals in the old forest, they are mostly in Gigante peninsula
species = (unique(young$sp_code)); length(unique(young$sp_code)) #37 species
old = old_allometry[old_allometry$dbh >= 500,]
old = old[old$sp_code %in% species,]; length(unique(old$sp_code)) #37 species

##graphing for ONLY giant trees species
#create a data.frame that has each species, m_rate, young forest values, old growth forest, so it is easier to graph
data.plots = data.frame("species" = species,
                        "young.biomass" = rep(0, length(species)), "old.biomass" = rep(0, length(species)),
                        "young.abund" = rep(0, length(species)), "old.abund" = rep(0, length(species)),
                        "y_abund.loCI" = rep(0, length(species)), "y_abund.hiCI" = rep(0, length(species)),
                        "o_abund.loCI" = rep(0, length(species)), "o_abund.hiCI" = rep(0, length(species)),
                        "young.crown" = rep(0, length(species)), "old.crown" = rep(0, length(species)),
                        "young.diam" = rep(0, length(species)), "old.diam" = rep(0, length(species)),
                        "y_diam.loCI" = rep(0, length(species)), "y_diam.hiCI" = rep(0, length(species)),
                        "o_diam.loCI" = rep(0, length(species)), "o_diam.hiCI" = rep(0, length(species)),
                        "young.ba" = rep(0, length(species)), "old.ba" = rep(0, length(species)),
                        "y_biomass.loCI" = rep(0, length(species)), "y_biomass.hiCI" = rep(0, length(species)),
                        "o_biomass.loCI" = rep(0, length(species)),"o_biomass.hiCI" = rep(0, length(species)),
                        "o_crown.hiCI" = rep(0, length(species)),
                        "o_crown.loCI" = rep(0, length(species)),"y_crown.hiCI" = rep(0, length(species)), "y_crown.loCI" = rep(0, length(species)),
                        "o_ba.loCI" = rep(0, length(species)), "o_ba.hiCI" = rep(0, length(species)),
                        "y_ba.loCI" = rep(0, length(species)),"y_ba.hiCI" = rep(0, length(species)))

#areas in ha #check this
y.ha = 4+1.5+27.04+4+6+1.92 #(dent+lang+25ha+ava+10ha_young+50ha_young) = 44.46ha
o.ha = 105.88 #4+47.64+6+6+6+36.24 #(10ha_old+50ha_old+pearson+zetek+drayton+gigante_peninsula) = 105.88ha

#pre-split row indices by quadrat_plot once, so each bootstrap draw is a fast lookup
young.quads = unique(young$quadrat_plot)
young.quad.idx = split(seq_len(nrow(young)), young$quadrat_plot)
old.quads = unique(old$quadrat_plot)
old.quad.idx = split(seq_len(nrow(old)), old$quadrat_plot)

funk =  function(x){
  print(x)
  y.use = young[young$sp_code == x,]
  o.use = old[old$sp_code == x,]
  #biomass
  data.plots[data.plots$species == x, "young.biomass"] <<- sum(y.use$AGB)/y.ha
  data.plots[data.plots$species == x, "old.biomass"] <<- sum(o.use$AGB)/o.ha
  #abund
  data.plots[data.plots$species == x, "young.abund"] <<- nrow(y.use)/y.ha
  data.plots[data.plots$species == x, "old.abund"] <<- nrow(o.use)/o.ha
  #crown area
  data.plots[data.plots$species == x, "young.crown"] <<- sum(y.use$allometric_crownarea)/y.ha
  data.plots[data.plots$species == x, "old.crown"] <<- sum(o.use$allometric_crownarea)/o.ha
  #use the 90th percentile of giant tree diameter
  data.plots[data.plots$species == x, "young.diam"] <<- quantile(y.use$dbh, probs = 0.90)
  data.plots[data.plots$species == x, "old.diam"] <<- quantile(o.use$dbh, probs = 0.90)
  #basal area
  data.plots[data.plots$species == x, "young.ba"] <<- sum(y.use$ba)/y.ha
  data.plots[data.plots$species == x, "old.ba"] <<- sum(o.use$ba)/o.ha
}

a = unlist(lapply(species, funk))

##CI of structure
n.boot = 999
mk.boot.mat = function() matrix(NA_real_, nrow = n.boot, ncol = length(species), dimnames = list(NULL, species))

biomass.young.boot = mk.boot.mat(); biomass.old.boot = mk.boot.mat()
abund.young.boot = mk.boot.mat(); abund.old.boot = mk.boot.mat()
crown.young.boot = mk.boot.mat(); crown.old.boot = mk.boot.mat()
diam.young.boot = mk.boot.mat(); diam.old.boot = mk.boot.mat()
ba.young.boot = mk.boot.mat(); ba.old.boot = mk.boot.mat()

for(j in 1:n.boot){
  #resample whole quadrats with replacement (not individual trees)
  boot.quads.young = sample(young.quads, size = length(young.quads), replace = TRUE)
  resampled_data_young = young[unlist(young.quad.idx[boot.quads.young]),]
  boot.quads.old = sample(old.quads, size = length(old.quads), replace = TRUE)
  resampled_data_old = old[unlist(old.quad.idx[boot.quads.old]),]
  
  for (x in species){
    y.use.resampled = resampled_data_young[resampled_data_young$sp_code == x,]
    o.use.resampled = resampled_data_old[resampled_data_old$sp_code == x,]
    
    #rare species can have zero individuals in this draw for one forest type = leave NA 
    if(nrow(y.use.resampled) > 0){
      biomass.young.boot[j, x] = sum(y.use.resampled$AGB)/y.ha
      abund.young.boot[j, x] = nrow(y.use.resampled)/y.ha
      crown.young.boot[j, x] = sum(y.use.resampled$allometric_crownarea)/y.ha
      diam.young.boot[j, x] = quantile(y.use.resampled$dbh, probs = 0.90)
      ba.young.boot[j, x] = sum(y.use.resampled$ba)/y.ha
    }
    
    if(nrow(o.use.resampled) > 0){
      biomass.old.boot[j, x] = sum(o.use.resampled$AGB)/o.ha
      abund.old.boot[j, x] = nrow(o.use.resampled)/o.ha
      crown.old.boot[j, x] = sum(o.use.resampled$allometric_crownarea)/o.ha
      diam.old.boot[j, x] = quantile(o.use.resampled$dbh, probs = 0.90)
      ba.old.boot[j, x] = sum(o.use.resampled$ba)/o.ha
    }
  }
}

q.lo = function(m) apply(m, 2, quantile, probs = 0.025, na.rm = TRUE)
q.hi = function(m) apply(m, 2, quantile, probs = 0.975, na.rm = TRUE)

data.plots$o_biomass.loCI = q.lo(biomass.old.boot)[data.plots$species]
data.plots$o_biomass.hiCI = q.hi(biomass.old.boot)[data.plots$species]
data.plots$y_biomass.loCI = q.lo(biomass.young.boot)[data.plots$species]
data.plots$y_biomass.hiCI = q.hi(biomass.young.boot)[data.plots$species]

data.plots$o_abund.loCI = q.lo(abund.old.boot)[data.plots$species]
data.plots$o_abund.hiCI = q.hi(abund.old.boot)[data.plots$species]

data.plots$y_abund.loCI = q.lo(abund.young.boot)[data.plots$species]
data.plots$y_abund.hiCI = q.hi(abund.young.boot)[data.plots$species]

data.plots$o_crown.loCI = q.lo(crown.old.boot)[data.plots$species]
data.plots$o_crown.hiCI = q.hi(crown.old.boot)[data.plots$species]
data.plots$y_crown.loCI = q.lo(crown.young.boot)[data.plots$species]
data.plots$y_crown.hiCI = q.hi(crown.young.boot)[data.plots$species]

data.plots$o_diam.loCI = q.lo(diam.old.boot)[data.plots$species]
data.plots$o_diam.hiCI = q.hi(diam.old.boot)[data.plots$species]

data.plots$y_diam.loCI = q.lo(diam.young.boot)[data.plots$species]
data.plots$y_diam.hiCI = q.hi(diam.young.boot)[data.plots$species]

data.plots$o_ba.loCI = q.lo(ba.old.boot)[data.plots$species]
data.plots$o_ba.hiCI = q.hi(ba.old.boot)[data.plots$species]
data.plots$y_ba.loCI = q.lo(ba.young.boot)[data.plots$species]
data.plots$y_ba.hiCI = q.hi(ba.young.boot)[data.plots$species]

#taxon
data.plots$taxon = unlist(lapply(data.plots$species, funk_taxon))

#add mortality rates
funk2 = function(x){ 
  use = match(toString(x), ML_mortality$sp)
  sp = ML_mortality[use,"ML_mortality"]
  return(sp)
}
data.plots$ML_mortality = unlist(lapply(species, funk2))

### graphs
# panel graphs
library(tibble)
library(tidyverse)

## what proportion of species have higher vs. lower total giant tree biomass in the older forests
tmp.data.plots = data.plots
tmp.data.plots$change.biomass = rep(NA,nrow(tmp.data.plots))

for(i in 1:nrow(tmp.data.plots)){
  if(tmp.data.plots[i, "young.biomass"]>tmp.data.plots[i, "old.biomass"]){
    tmp.data.plots[i, "change.biomass"] = 1 #decresed
  }else{tmp.data.plots[i, "change.biomass"] = 2} #incresed
}

sum(data.plots$old.biomass > data.plots$young.biomass) 
sum(data.plots$old.abund > data.plots$young.abund) 
sum(data.plots$old.diam > data.plots$young.diam)

#### There are some species that are not in the young forest: eg. Ceiba, Vataer
names(data.plots)
alp = 0.8

### *Figure S5a: abundance in secondary(y) and old-growth(x) with Confidence intervals
model1 = lm(young.abund ~ old.abund,data = data.plots)
summary(model1)

S.g1 = ggplot(data = data.plots, aes(x = old.abund, y = young.abund)) +
  geom_point(data = data.plots, size = 3, aes(x = old.abund, y = young.abund, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  scale_y_log10(limits = c(0.009,5), expression(bold("Abundance secondary forest"~(''~ha^{-1}))))+
  scale_x_log10(limits = c(0.009,5), expression(bold("Abundance old growth forest"~(''~ha^{-1}))))+
  guides(color = "none")+ theme_basis+
  geom_errorbar(aes(x = old.abund, y = young.abund, ymax = y_abund.hiCI, ymin = y_abund.loCI), size = 0.1)+
  geom_errorbar(aes(x = old.abund, y = young.abund, xmax = o_abund.hiCI, xmin = o_abund.loCI), size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )

points_to_label1 <- data.plots %>%
  filter(species %in% c("cordal", "pse1se", "pri2co", "quaras","poular","calolo")) %>%
  mutate(number = c(1:6))

### *Figure 3a: abundance in secondary(y) and old-growth(x)
g1.main = ggplot(data = data.plots, aes(x = old.abund, y = young.abund)) +
  geom_point(data = data.plots, size = 3, aes(x = old.abund, y = young.abund, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  geom_text(data = points_to_label1, aes(label = round(number, 1)), vjust = c(-1,-1,-1,-1,-1,1.7), size = 3)+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  scale_x_log10(limits = c(0.01,2.5), expression(bold("Number of giant trees, old-growth"~(''~ha^{-1}))))+
  scale_y_log10(limits = c(0.01,2.5), expression(bold("Number of giant trees, secondary"~(''~ha^{-1}))))+
  guides(color = "none")+ theme_basis+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )+
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) 

#how many points are to the right of the 1:1 line (i.e not recovered)
sum(data.plots$old.abund > data.plots$young.abund)/
  nrow(data.plots)

### *Figure S5b: max diameter in secondary(y) and old-growth(x) with confidence intervals
model2 = lm(log(young.diam)~log(old.diam),data = data.plots)
summary(model2)
S.g2 = ggplot(data = data.plots, aes(x = old.diam, y = young.diam)) +
  geom_point(data = data.plots, size = 3, aes(x = old.diam, y = young.diam, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  theme_basis + guides(color = "none")+ 
  scale_x_log10(expression(bold("Max diameter old growth forest (mm)")), limits = c(500,3225)) + 
  scale_y_log10(expression(bold("Max diameter secondary forest (mm)")), limits = c(500,3225))+
  geom_errorbar(aes(x = old.diam, y = young.diam, xmax = o_diam.hiCI, xmin = o_diam.loCI), size = 0.1)+
  geom_errorbar(aes(x = old.diam, y = young.diam, ymax = y_diam.hiCI, ymin = y_diam.loCI), size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )

points_to_label2 <- data.plots %>%
  filter(species %in% c("pse1se", "ceibpe", "huracr","anacex", "diptpa")) %>%
  mutate(number = c(7,1,8,9,10))

### *Figure 3b: max diameter in secondary(y) and old-growth(x)
g2.main = ggplot(data = data.plots, aes(x = old.diam, y = young.diam)) +
  geom_point(data = data.plots, size = 3, aes(x = old.diam, y = young.diam, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + 
  geom_text(data = points_to_label2, aes(label = round(number, 1)), vjust = c(-1,-1,-1,-1,-1), size = 3)+
  scale_color_gradient(low= "blue2",high = "red1")+
  theme_basis + guides(color = "none")+
  scale_x_log10(expression(bold("Max diameter, old-growth (mm)")), limits = c(500,3225)) + 
  scale_y_log10(expression(bold("Max diameter, secondary (mm)")), limits = c(500,3225))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )+
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

### *Figure S5d: basal area in m2 in secondary(y) and old-growth(x) with conficende intervals
model3 = lm(log(young.ba)~log(old.ba),data = data.plots)
summary(model3)
S.g3 = ggplot(data = data.plots, aes(x = old.ba, y = young.ba)) +
geom_point(data = data.plots, size = 3, aes(x = old.ba, y = young.ba, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  theme_basis + guides(color = "none")+ 
  scale_x_log10(limits = c(0.002,3), expression(bold("Basal area old growth forest"~(''~m^{2}~ha^{-1})))) + 
  scale_y_log10(limits = c(0.002,3), expression(bold("Basal area secondary forest"~(''~m^{2}~ha^{-1}))))+
  geom_errorbar(aes(x = old.ba, y = young.ba, xmax = o_ba.hiCI, xmin = o_ba.loCI), size = 0.1)+
  geom_errorbar(aes(x = old.ba, y = young.ba, ymax = y_ba.hiCI, ymin = y_ba.loCI), size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )

### *Figure S5c: biomass in secondary(y) and old-growth(x) with confidence intervals
model4 = lm(log(young.biomass)~log(old.biomass),data = data.plots)
summary(model4)
S.g4 = ggplot(data = data.plots, aes(x = old.biomass, y = young.biomass)) +
  geom_point(data = data.plots, size = 3, aes(x = old.biomass, y = young.biomass, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  theme_basis + guides(color = "none")+ 
  scale_x_log10(expression(bold("Biomass old growth forest"~('Mg'~ha^{-1})))) + 
  scale_y_log10(expression(bold("Biomass secondary forest"~('Mg'~ha^{-1}))))+
  geom_errorbar(aes(x = old.biomass, y = young.biomass, xmax = o_biomass.hiCI, xmin = o_biomass.loCI), size = 0.1)+
  geom_errorbar(aes(x = old.biomass, y = young.biomass, ymax = y_biomass.hiCI, ymin = y_biomass.loCI), size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )

points_to_label3 <- data.plots %>%
  filter(species %in% c("cordal", "diptpa", "pri2co", "pse1se","calolo")) %>%
  mutate(number = c(7,1,2,3,4))

### *Figure 3c: biomass in secondary(y) and old-growth(x)
g4.main = ggplot(data = data.plots, aes(x = old.biomass, y = young.biomass)) +
  geom_point(data = data.plots, size = 3, aes(x = old.biomass, y = young.biomass, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  geom_text(data = points_to_label3, aes(label = round(number, 1)), vjust = c(-1,-1,-1,-1,-1), size = 3)+
  scale_color_gradient(low= "blue2",high = "red1",name = expression("Giant tree \nmortality"('%'~yr^{-1}), sep =""))+
  guides(color = guide_colorbar(
    direction = "vertical", title.position = "top", barwidth = 1, barheight = 6))+
  theme(legend.position = "right")+
  theme_basis +
  scale_x_log10(expression(bold("Biomass, old-growth"~('Mg'~ha^{-1})))) + 
  scale_y_log10(expression(bold("Biomass, secondary"~('Mg'~ha^{-1}))))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )+
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

### *Figure S5e: crown area (m2) in secondary(y) and old-growth(x)
model5 = lm(log(young.crown)~log(old.crown),data = data.plots)
summary(model5)
S.g5 = ggplot(data = data.plots, aes(x = old.crown, y = young.crown)) +
  geom_point(data = data.plots, size = 3, aes(x = old.crown, y = young.crown, color = log(ML_mortality)),
             stroke = 0, alpha = alp) + guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"(yr^{-1})))+
  theme_basis + guides(color = "none")+ 
  scale_x_log10(expression(bold("Crown area old growth forest"~(''~m^{2}~ha^{-1}))), limits = c(1,1000)) + 
  scale_y_log10(expression(bold("Crown secondary forest"~(''~m^{2}~ha^{-1}))), limits = c(1,1000))+
  geom_errorbar(aes(x = old.crown, y = young.crown, xmax = o_crown.hiCI, xmin = o_crown.loCI), size = 0.1)+
  geom_errorbar(aes(x = old.crown, y = young.crown, ymax = y_crown.hiCI, ymin = y_crown.loCI), size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )

###1:1 species that are significantly different 
#how many species cross the 1:1 line
value.x = data.plots$old.biomass
value.x.loCI = data.plots$o_biomass.loCI
value.x.hiCI = data.plots$o_biomass.hiCI

value.y = data.plots$young.biomass
value.y.loCI = data.plots$y_biomass.loCI
value.y.hiCI = data.plots$y_biomass.hiCI

crosses_x <- (value.x.loCI <= value.y) & (value.y <= value.x.hiCI)
crosses_y <- (value.y.loCI <= value.x) & (value.x <= value.y.hiCI)
#At least one axis,  each species counted once, no double counting
count <- crosses_x | crosses_y
sum(count)   # total unique species touching the line

#species that do not touch the line in either axis
non_touching <- data.plots[!count, ]
non_touching$below_line <- non_touching$old.biomass > non_touching$young.biomass
all(non_touching$below_line)
non_touching[!non_touching$below_line, ]
sum(non_touching$below_line) # below the line
100*(sum(non_touching$below_line)/length(which(value.x > value.y)))
sum(!non_touching$below_line) # above the line

### woody residence time in old-growth vs. mature secondary 
all.comm.giant.m.rate = 0.0206576

#assigning community mortality rate to the trees lower than 50cm from our species
load("ML_mortality_sp39.less50.Rdata") #from Other.spp.mortality.R
sp39.less50.m.rate = ML_mortality_sp39.less50$ML_mortality

#mortality rate of all individuals that are not part of our 39species; from: 
load("ML_mortality_other.more50.Rdata")
spother.more50.m.rate = ML_mortality_other.more50$ML_mortality

load("ML_mortality_other.less50.Rdata")
spother.less50.m.rate = ML_mortality_other.less50$ML_mortality

#### young forest ####
#divide the data in 39sp and other species and size class, and calculate their biomass, the reference is the quadrats 
#in the 50ha that are young forest
wrts.young = rep(0,1000) #for CIs

#pre-split row indices by quadrat once, so each bootstrap draw is a fast lookup
young.df = as.data.frame(young)
young.quads = unique(young.df$quadrat_plot)
young.quad.idx = split(seq_len(nrow(young.df)), young.df$quadrat_plot)

for(k in 1:1000){
  if(k == 1){
    dataset = young
  }else{
    #resample whole quadrats with replacement (not individual trees),
    #so trees within the same quadrat move together and spatial autocorrelation is preserved
    boot.quads = sample(young.quads, size = length(young.quads), replace = TRUE)
    dataset = young.df[unlist(young.quad.idx[boot.quads]), ]
  }
  biomass.39sp.less50 = data.frame("species" = species.to.use, "biomass.39sp.10_25mm" = rep(0, length(species.to.use)),
                                   "biomass.39sp.25_100mm" = rep(0, length(species.to.use)),
                                   "biomass.39sp.100_200mm" = rep(0, length(species.to.use))) #<less than 20cm
  
  funk.wrt = function(x){ #39 species, this calculates the biomass in each class on the young forest of the 50ha plot, then use this value to extrapolate to other plots (extra biomass)
    
    use = dataset[dataset$sp == x & dataset$plot == "plot50ha",] #species x in the young are of the 50ha plot
    use2 = sum(use[use$dbhmm.1.3>=10 & use$dbhmm.1.3<25, "AGB"]) #between 10-25mm
    use3 = sum(use[use$dbhmm.1.3>=25 & use$dbhmm.1.3<100, "AGB"]) #between 25-100mm
    use4 = sum(use[use$dbhmm.1.3>=100 & use$dbhmm.1.3<200, "AGB"]) #between 100-200mm

    biomass.39sp.less50[biomass.39sp.less50$species == x, "biomass.39sp.10_25mm"] <<- use2
    biomass.39sp.less50[biomass.39sp.less50$species == x, "biomass.39sp.25_100mm"] <<- use3
    biomass.39sp.less50[biomass.39sp.less50$species == x, "biomass.39sp.100_200mm"] <<- use4
  }
  t = lapply(species.to.use, funk.wrt) 
  
  #now, calculate extra biomass for other species for the same class sizes for the 50ha plot young forest
  temp = dataset[dataset$plot == "plot50ha",]
  dataset.other = temp[-(which(temp$sp %in% species.to.use)),]; length(unique(dataset.other$sp)) #5172 rows and 139spp
  biomass.other.sp.10_25mm = sum(dataset.other[dataset.other$dbhmm.1.3>=10 & dataset.other$dbhmm.1.3<25, "AGB"]) #between 10-25mm
  biomass.other.sp.25_100mm = sum(dataset.other[dataset.other$dbhmm.1.3>=25 & dataset.other$dbhmm.1.3<100, "AGB"]) #between 25-100mm
  biomass.other.sp.100_200mm = sum(dataset.other[dataset.other$dbhmm.1.3>=100 & dataset.other$dbhmm.1.3<200, "AGB"]) #between 100-200mm

  ## because we have different size structure (plots have different dbh at which they start measuring trees), we did a weighting 
  # biomass correction:
  # for plots that did not have a size class we multiplied by 
  # (biomass of the type of forests known for other forests that have the size class)*(total area without that size class)/(area in which the numerator biomass was calculated from)
  # eg. biomass 1-20cm in the young area of the 50ha plot, weighted for other sites as:
  # (biomass 1-20cm)*(38.54ha/1.92ha) which is 
  # all area without this size class(25ha, 10ha, knight-lang, and AVA)/1.92ha is the total of the 50ha that is young and has this size 1-19cm included
  #10-25mm 25ha(27.04), 10hayoung(6ha),  ava(4), knight(1.5ha)
  #25-100mm 25ha, 10hayoung, ava
  #100-200mm 10hayoung
  #biomass weighted by area (see above the hectares used) for less than 50cm, we do not weight 200-500mm because we have all plots with this data
  biomass.w.other.sp.10_25mm = biomass.other.sp.10_25mm * (27.04+6+4+1.5)/(1.92)
  biomass.w.other.sp.25_100mm = biomass.other.sp.25_100mm * (27.04+6+4)/(1.92)
  biomass.w.other.sp.100_200mm = biomass.other.sp.100_200mm * (6)/(1.92) #10ha starts measuring at 200mm
  biomass.other.sp.extra.less50.young = biomass.w.other.sp.10_25mm + biomass.w.other.sp.25_100mm + biomass.w.other.sp.100_200mm
  biomass.other.sp.less50.young = sum(na.omit(dataset[!dataset$sp %in% species.to.use & dataset$dbhmm.1.3<500, "AGB"]))
  
  #biomass weighted by area (see above the hectares used) for more than 50cm, we do not weight 200-500mm because we have all plots with this data
  biomass.39sp.less50$biomass.w.39sp.10_25mm = biomass.39sp.less50$biomass.39sp.10_25mm * (27.04+6+4+1.5)/(1.92)
  biomass.39sp.less50$biomass.w.39sp.25_100mm = biomass.39sp.less50$biomass.39sp.25_100mm * (27.04+6+4)/(1.92)
  biomass.39sp.less50$biomass.w.39sp.100_200mm = biomass.39sp.less50$biomass.39sp.100_200mm * (6)/(1.92) #10ha starts measuring at 200mm
  biomass.39sp.extra.less50.young = sum(biomass.39sp.less50[,5:7])
  biomass.39sp.less50.young = sum(na.omit(dataset[dataset$sp %in% species.to.use & dataset$dbhmm.1.3<500, "AGB"]))
  
  ## WOODY RESIDENCE TIMES of the species on the whole community ##
  
  #above-ground biomass for the whole community/ha
  #total hectares: 25ha 27.04ha; 10ha.young 6ha; AVA 4ha; knight-lang 1.5ha; 1.92 young from 50ha
  hectares.used.young = 40.46
  
  #total biomass
  young.other.more50 = dataset[!dataset$sp%in%species.to.use & dataset$dbhmm.1.3>=500,]
  young.other.less50 = dataset[!dataset$sp%in%species.to.use & dataset$dbhmm.1.3<500,]
  
  young.sp39.more50 = dataset[dataset$sp%in%species.to.use & dataset$dbhmm.1.3>=500,]
  young.sp39.less50 = dataset[dataset$sp%in%species.to.use & dataset$dbhmm.1.3<500,]
  young.sp39.less50$ML_mortality = sp39.less50.m.rate
  
  agb.other.more50.young = sum(na.omit(young.other.more50$AGB));agb.other.more50.young
  agb.39sp.more50.young = sum(na.omit(young.sp39.more50$AGB));agb.39sp.more50.young
  
  #all biomass + extras in hectares
  agb.all.young.ha = (agb.other.more50.young + agb.39sp.more50.young + biomass.other.sp.less50.young + biomass.39sp.less50.young +
                        biomass.other.sp.extra.less50.young + biomass.39sp.extra.less50.young)/hectares.used.young; agb.all.young.ha 
  
  #add the real mortality rate per species, and size 
  funk = function(x){ 
    use = match(toString(x), data$sp)
    mr = data[use, "ML_mortality"]
    return(mr)
  }
  young.sp39.more50$ML_mortality = unlist(lapply(young.sp39.more50$sp, funk))
  
  #woody residence time of all the community that includes all trees, with their corresponding mrates
  sp39.wrt.less50.young = c()
  sp39.wrt.more50.young = c()
  
 #calcualte woody residence times
  species.to.use = species.to.use
  funk.wrt = function(x, y){ #for the 39 species we have

    use1.young = young.sp39.less50[young.sp39.less50$sp == x,] #species x and individuals with dbh<50cm from our species
    use2.young = y[y$sp == x,] #species x and individuals with dbh>50cm from our species
    use3.young = biomass.39sp.less50[biomass.39sp.less50$sp == x,] #extra biomass for smaller size classes, use the weighted biomass per species
    agb.species.less50.young = (sum(na.omit(use1.young$AGB)) + sum(na.omit(use3.young[,5:7])))/hectares.used.young #Mg/ha
    agb.species.more50.young = sum(na.omit(use2.young$AGB))/hectares.used.young #Mg/ha
    
    if(length(agb.species.less50.young * unique(use1.young$ML_mortality))==0 & length(agb.species.more50.young * unique(use2.young$ML_mortality))==0){
      num.wrt.less50.young = 0 
      num.wrt.more50.young = 0
    } else if (length(agb.species.more50.young * unique(use2.young$ML_mortality))==0){
      num.wrt.less50.young = (agb.species.less50.young * unique(use1.young$ML_mortality)) 
      num.wrt.more50.young = 0 
    } else if (length(agb.species.less50.young * unique(use1.young$ML_mortality))==0){
      num.wrt.less50.young = 0 
      num.wrt.more50.young = (agb.species.more50.young * unique(use2.young$ML_mortality)) 
    }else{
      num.wrt.less50.young = (agb.species.less50.young * unique(use1.young$ML_mortality))
      num.wrt.more50.young = (agb.species.more50.young * unique(use2.young$ML_mortality)) #sp-specific giant mortality rate 
    }
    sp39.wrt.less50.young[x] <<- num.wrt.less50.young
    sp39.wrt.more50.young[x] <<- num.wrt.more50.young
  }
  t = lapply(species.to.use,funk.wrt, y = young.sp39.more50)
  
  agb.others.less50.ha = (biomass.other.sp.extra.less50.young+biomass.other.sp.less50.young)/hectares.used.young 
  agb.other.more50.young.ha = agb.other.more50.young/hectares.used.young 
  wrt2.young = (agb.all.young.ha)/(((sum(sp39.wrt.less50.young)+sum(sp39.wrt.more50.young))+(agb.others.less50.ha*spother.less50.m.rate)+(agb.other.more50.young.ha*spother.more50.m.rate))) #sp-specific wrt.less.50+wrt.more.50+other species biomass/community rate
  wrt2.young 
  wrts.young[k] = wrt2.young #k = 1 is the original
  print(k)
}

#### mature (old) forest ####
rm(biomass.39sp.less50)
dataset = mature
wrts.mature = rep(0,1000) #for CIs

#pre-split row indices by quadrat once, so each bootstrap draw is a fast lookup
mature.df = as.data.frame(mature)
mature.quads = unique(mature.df$quadrat_plot)
mature.quad.idx = split(seq_len(nrow(mature.df)), mature.df$quadrat_plot)

for(m in 1:1000){
  dataset = mature
  if(m == 1){
    dataset = dataset
  }else{
    #resample whole quadrats with replacement (not individual trees),
    #so trees within the same quadrat move together and spatial autocorrelation is preserved
    boot.quads = sample(mature.quads, size = length(mature.quads), replace = TRUE)
    dataset = mature.df[unlist(mature.quad.idx[boot.quads]), ]
  }
  species.to.use.mature = unique(ML_mortality$sp)
  extras.dataset = dataset[dataset$plot == "Gigante" | dataset$plot == "plot50ha",] 
  biomass.39sp.less50 = data.frame("species" = species.to.use.mature, "biomass.39sp.10_200mm" = rep(0, length(species.to.use.mature)))
  
  funk.wrt = function(x){ #39 species, this calculates the biomass in each class on the 50ha_old and Gigante peninsula (they have data from 10-200mm)
    #then use this value to extrapolate to other plots (extra biomass)
    
    use = extras.dataset[extras.dataset$sp == x,] #species x in the young are of the 50ha plot
    use2 = sum(na.omit(use[use$dbhmm.1.3>=10 & use$dbhmm.1.3<200, "AGB"])) #between 10-200mm
    biomass.39sp.less50[biomass.39sp.less50$species == x, "biomass.39sp.10_200mm"] <<- use2
  }
  t = lapply(species.to.use.mature, funk.wrt) #fun(x,  funk)
  
  #now, calculate biomass for other species for the same class sizes for the 50ha plot young forest
  dataset.other = extras.dataset[-(which(extras.dataset$sp %in% species.to.use.mature)),]; length(unique(dataset.other$sp)) #139,543rows and 298spp
  biomass.other.sp.10_200mm = sum(na.omit(dataset.other[dataset.other$dbhmm.1.3>=10 & dataset.other$dbhmm.1.3<200, "AGB"])) #between 10-25mm
  
  ## because we have different size structure (plots have different dbh at which they start measuring trees), we did a weighting 
  # biomass correction:
  # for plots that did not have a size class we multiplied by 
  # (biomass of the type of forests known for other forests that have the size class)*(total area without that size class)/(area in which the numerator biomass was calculated from)
  # eg. biomass 1-20cm in the old area of the 50ha plot and Gigante Peninsula (36.24+(50-1.92young-0.44mix)), weighted for other sites as:
  # (biomass 1-20cm)*(sites area/84.32) which is 
  # all area without this size class(25ha, 10ha, knight-lang, and AVA)/1.92ha is the total of the 50ha that is young and has this size 1-19cm included
  #10-200mm 10ha old(4ha), Pearson (6ha), Zetek (6ha) and Drayton (6ha) 
  #biomass weighted by area (see above the hectares used) for less than 50cm, we do not weight 200-500mm because we have all plots with this data
  biomass.w.other.sp.10_200mm = biomass.other.sp.10_200mm * (4+6+6+6)/(36.24+47.06)
  biomass.other.sp.extra.less50.mature = biomass.w.other.sp.10_200mm
  biomass.other.sp.less50.mature = sum(na.omit(dataset[!dataset$sp %in% species.to.use.mature & dataset$dbhmm.1.3<500, "AGB"]))
  
  #biomass weighted by area (see above the hectares used) for more than 50cm, we do not weight 200-500mm because we have all plots with this data
  biomass.39sp.less50$biomass.w.39sp.10_200mm = biomass.39sp.less50$biomass.39sp.10_200mm * (4+6+6+6)/(36.24+47.06)
  biomass.39sp.extra.less50.mature = sum(biomass.39sp.less50$biomass.w.39sp.10_200mm)
  biomass.39sp.less50.mature = sum(na.omit(dataset[dataset$sp %in% species.to.use.mature & dataset$dbhmm.1.3<500, "AGB"]))
  
  ## WOODY RESIDENCE TIMES of the species on the whole community ##
  
  #above-ground biomass for the whole community/ha
  #biomass is 47% C (IPCC, 2006); this study mature biomass: 161.58 Mg/ha-1
  #total hectares: 50ha.old 47.64 (50-1.92young-0.44mixed)ha; 10ha.old 4ha; Pearson 6ha; zetek 6ha; drayton 6ha, gigante peninsula 36.24
  hectares.used.mature = 105.88
  
  #total biomass
  mature.other.more50 = dataset[!dataset$sp%in%species.to.use.mature & dataset$dbhmm.1.3>=500,]
  mature.other.less50 = dataset[!dataset$sp%in%species.to.use.mature & dataset$dbhmm.1.3<500,]
  
  mature.sp39.more50 = dataset[dataset$sp%in%species.to.use.mature & dataset$dbhmm.1.3>=500,]
  mature.sp39.less50 = dataset[dataset$sp%in%species.to.use.mature & dataset$dbhmm.1.3<500,]
  mature.sp39.less50$ML_mortality = sp39.less50.m.rate #0.0300364
  
  agb.other.more50.mature = sum(na.omit(mature.other.more50$AGB));agb.other.more50.mature
  agb.39sp.more50.mature = sum(na.omit(mature.sp39.more50$AGB));agb.39sp.more50.mature
  
  #all biomass + extras in hectares
  agb.all.mature.ha = (agb.other.more50.mature + agb.39sp.more50.mature + biomass.other.sp.less50.mature + biomass.39sp.less50.mature +
                         biomass.other.sp.extra.less50.mature + biomass.39sp.extra.less50.mature)/hectares.used.mature; agb.all.mature.ha #161.58
    #add the real mortality rate per species, and size 
  funk = function(x){ 
    use = match(toString(x), data$sp)
    mr = data[use, "ML_mortality"]
    return(mr)
  }
  mature.sp39.more50$ML_mortality = unlist(lapply(mature.sp39.more50$sp, funk))
  
  #woody residence time of all the community that includes all trees, with their corresponding mrates
  sp39.wrt.less50.mature= c()
  sp39.wrt.more50.mature = c()
  
  #calculate the woody residence times
  funk.wrt2 = function(x, y){ #for the 39 species we have
    use1.mature = mature.sp39.less50[mature.sp39.less50$sp == x,] #species x and individuals with dbh<50cm from our species
    use2.mature = y[y$sp == x,] #species x and individuals with dbh>50cm from our species
    use3.mature = biomass.39sp.less50[biomass.39sp.less50$species == x,] #extra biomass for smaller size classes, use the weighted biomass per species
    agb.species.less50.mature = (sum(na.omit(use1.mature$AGB)) + (use3.mature$biomass.w.39sp.10_200mm))/hectares.used.mature #Mg/ha
    agb.species.more50.mature = sum(na.omit(use2.mature$AGB))/hectares.used.mature #Mg/ha
    
    if(length(agb.species.less50.mature * unique(use1.mature$ML_mortality))==0 & length(agb.species.more50.mature * unique(use2.mature$ML_mortality))==0){
      num.wrt.less50.mature = 0 
      num.wrt.more50.mature = 0
    } else if (length(agb.species.more50.mature * unique(use2.mature$ML_mortality))==0){
      num.wrt.less50.mature = (agb.species.less50.mature * unique(use1.mature$ML_mortality)) #community mortality rates for trees <50cm of our species(~3%)
      num.wrt.more50.mature = 0 
    } else if (length(agb.species.less50.mature * unique(use1.mature$ML_mortality))==0){
      num.wrt.less50.mature = 0 
      num.wrt.more50.mature = (agb.species.more50.mature * unique(use2.mature$ML_mortality)) 
    }else{
      num.wrt.less50.mature = (agb.species.less50.mature * unique(use1.mature$ML_mortality)) #community mortality rates for trees <50cm of our species(~3%)
      num.wrt.more50.mature = (agb.species.more50.mature * unique(use2.mature$ML_mortality)) #sp-specific giant mortality rate 
    }
    sp39.wrt.less50.mature[x] <<- num.wrt.less50.mature
    sp39.wrt.more50.mature[x] <<- num.wrt.more50.mature
  }
  t = lapply(species.to.use.mature, funk.wrt2, y = mature.sp39.more50) #fun(x,  funk, y)
  
  agb.others.less50.ha = (biomass.other.sp.extra.less50.mature+biomass.other.sp.less50.mature)/hectares.used.mature 
  agb.other.more50.mature.ha = agb.other.more50.mature/hectares.used.mature 
  wrt2.mature = (agb.all.mature.ha)/(((sum(sp39.wrt.less50.mature)+sum(sp39.wrt.more50.mature))+(agb.others.less50.ha*spother.less50.m.rate)+(agb.other.more50.mature.ha*spother.more50.m.rate))) #sp-specific wrt.less.50+wrt.more.50+other species biomass/community rate
  wrt2.mature 
  wrts.mature[m] = wrt2.mature #k = 1 is the original
  print(m)
}
m = 1; wrts.mature[m]

#how much bigger is the wrt the mature is over the secondary [(big - small)/small]*100
((wrts.mature[1]-wrts.young[1])/wrts.young[1])*100 #5.12% bigger in old-growth forest

### *Figure 3d: Woody residence times and old-growth and secondary forest
to.plot = data.frame("Type" = c("Old-growth", "Mature Secondary"), "wrt" = c(wrts.mature[1], wrts.young[1]), 
                     "CI.up" = c(quantile(wrts.mature[-1],0.975), quantile(wrts.young[-1],0.975)),
                     "CI.low" = c(quantile(wrts.mature[-1],0.025), quantile(wrts.young[-1],0.025)),
                     "Age" = c("120", "500+"))
breaks <- seq(44, 51, by = 1)
labels <- ifelse(breaks %% 5 == 0, as.character(breaks), "")

wrt.forest = ggplot(to.plot, aes(x = Type, y = wrt)) +
  geom_point(aes(color = "blue"), shape = 19)+ theme_basis +
  scale_x_discrete(labels = c("ca.80-135",">400"))+
  scale_y_continuous(limits = c(43,51), breaks = breaks, labels = labels)+
  labs(x = "Forest Age" , y = "Woody residence time (yr)")+
  theme(legend.position = "none") +
  geom_errorbar(aes(x = Type, y = wrt, ymax = CI.up ,ymin = CI.low, color= "blue"), 
                width = 0.2)
