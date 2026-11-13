########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Woody residence time predicted change when giant tree species start dying at an average community giant tree mortality rate
## Controls for the size class distribution
## Produces: Figures S2b, 2
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggplot2); library(dplyr)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12, face ="bold")) #, face ="bold"

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                       simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                       !is.na(simple_combined$dbhmm.1.3),]

##select all data from the census intervals: all plots
temp.all.plots.all.sp = simple_combined_sp[simple_combined_sp$openyear == 2005 & simple_combined_sp$endyear == 2010 |
                                             simple_combined_sp$openyear == 2003 & simple_combined_sp$endyear == 2008 |
                                             simple_combined_sp$openyear == 2004 & simple_combined_sp$endyear == 2014 |
                                             simple_combined_sp$openyear == 2000 & simple_combined_sp$endyear == 2014,]

#select only the alive trees
temp.all.plots.all.sp = temp.all.plots.all.sp[temp.all.plots.all.sp$mortality == 0,]

#split every plot into "young/secondary" (baseline) vs "old-growth" (old)
#SECONDARY FORESTS
#lang.k, plot25ha, AVA, part of 10ha
plotlang = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plotKnight",]
plot25ha = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plot25ha",]
plotava = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plotAVA",]

#the 10ha plot has one part that is mature forest (remove it and use it below)
plot10ha.temp = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plot10ha",]

#around 420 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 420m there is mature forest, to the right of 420m it is secondary.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 780m there is young forest, to the right of 780m it has both secondary/mature.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and above 50m in the y-axis it is young forest.
load("all2.rdata")
plot10ha.1 = plot10ha.temp[plot10ha.temp$gx>250 & plot10ha.temp$gx<750,] 
plot10ha.2 = plot10ha.temp[plot10ha.temp$gx>800 & plot10ha.temp$gy>50,] 
plot10ha.old = plot10ha.temp[plot10ha.temp$gx<200,] 
plot10ha.young = rbind(plot10ha.1, plot10ha.2)

##where is the young forest areas in the 50ha plot?
#loads the topographic types for each quadrat from Harms et al., 2001 J.of Ecology
harms = read.table("harms.txt",h=T)
harms$quadrat = paste(harms$GX20/20, (harms$GY20/20), sep = "") #code to name the quadrats on the 50haplot
young.forest.50ha = harms[harms$Habitat == "young",] #48 quadrats
mixed.forest.50ha = c("2123","2222","2422","2522","2622","2722","3320","3520","3620","3722","3824") #this is mixed in vicinity of the young forest #11quads
#Add young forest data to the secondary forest data from above (young 50ha is 1.92ha = 48quadrats*400m2)
plot50ha = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plot50ha",]
plot50ha.young = plot50ha[which(plot50ha$quadrat %in% young.forest.50ha$quadrat),]
plot50ha.old = plot50ha[-which(plot50ha$quadrat %in% young.forest.50ha$quadrat),]
plot50ha.old = plot50ha.old[-which(plot50ha.old$quadrat %in% mixed.forest.50ha),]
#bind the young/secondary datasets
baseline = rbind(plotlang, plot25ha, plotava, plot10ha.young, plot50ha.young)

#OLD-GROWTH FORESTS
#add data from Pearson, Zetek and Drayton as mature forests
plotPearson = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plotPearson",]
plotZetek = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plotZetek",]
plotDrayton = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "plotDrayton",]
#gigante peninsula plot >600 years
plotgigante = temp.all.plots.all.sp[temp.all.plots.all.sp$plot == "Gigante",]
plotgigante = plotgigante[plotgigante$gigante_plot >= 0,] #trees that are in plots that measure trees >=1cm
plotgigante = plotgigante[plotgigante$dbhmm.1.3 > 0,]
#bind the old-growth datasets
old = rbind(plot10ha.old, plot50ha.old, plotPearson, plotZetek, plotDrayton, plotgigante)

#add AGB, wood density, and taxon
add.taxon.wd.agb = function(dataset){
  load("bci.spptable.rdata")
  bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")
  
  funk = function(species){
    use = match(toString(species), bci.spptable$sp)
    bci.spptable[use, "taxon"]
  }
  dataset$taxon = unlist(lapply(dataset$sp, funk))
  
  funk = function(species){
    use = match(toString(species), bci.spptable$sp)
    bci.spptable[use, "Family"]
  }
  dataset$family = unlist(lapply(dataset$sp, funk))
  
  jw_data = read.table("WD_species.txt", h = T)
  jw_data = data.frame(jw_data)
  jw_data$SP. = tolower(jw_data$sp6)
  
  funk = function(species){
    use = match(toString(species), jw_data$SP.)
    jw_data[use, "wd100.mean"]
  }
  dataset$wd = unlist(lapply(dataset$sp, funk))
  dataset = as.data.frame(dataset)
  
  #extras wood from ctfs data when the species did not have reported wood density in the previous dataset
  load("wsg.ctfs.Rdata")
  ctfs = wsg.ctfs3
  ctfs$taxon = paste(ctfs$genus, ctfs$species, sep=" ")
  with.nas = dataset[is.na(dataset$wd),]
  unique.with.nas = unique(with.nas$taxon)
  ctfs[which(ctfs$taxon == "Acacia melanoceras"), "taxon"] = "Vachellia melanoceras"
  
  for (tax in unique.with.nas){
    use = match(toString(tax), ctfs$taxon)
    wd2 = ctfs[use,"wsg"]
    dataset[which(dataset$taxon == tax), "wd"] = wd2
  }
  
  #for species with no wood density value: the mean of all other species in the site was used
  with.nas = dataset[is.na(dataset$wd),]
  unique.with.nas = unique(with.nas$sp)
  dataset[which(dataset$sp %in% unique.with.nas), "wd"] = mean(na.omit(dataset$wd))
  
  #calculate AGB - Chave et al., 2015 (Eq. 7), moist tropical forest, D in cm, wsg in g/cm3, Kg/1000 -> Mg
  palms = dataset[which(dataset$family == "Arecaceae"),]
  palms$AGB = (exp(-3.3488 + 2.7484*log(palms$dbhmm.1.3/10))/1000) * (1.188714)
  
  trees = dataset[which(dataset$family != "Arecaceae"),]
  E = 0.051763982 #for BCI
  trees$AGB = (exp(-1.803-0.976*E+(0.976*log(trees$wd))+
                     (2.673*log(trees$dbhmm.1.3/10))-
                     (0.0299*(log(trees$dbhmm.1.3/10))^2))/1000)
  
  rbind(trees, palms)
}

baseline = add.taxon.wd.agb(baseline)
old = add.taxon.wd.agb(old)
#bind the data
temp.all.plots.all.sp = rbind(baseline, old)

### analysis ### 
#add the species we have mortality data, and do the analysis only for that= 39 species
load("ML_mortality_50CI.Rdata")
species.to.use = unique(ML_mortality$sp)
data = ML_mortality

all.comm.giant.m.rate = 0.0206576

#assigning community mortality rate to the trees lower than 50cm from our species
load("ML_mortality_sp39.less50.Rdata") #from Other.spp.mortality.R
sp39.less50.m.rate = ML_mortality_sp39.less50$ML_mortality
#mortality rate of all individuals that are not part of our 39species that have dbh>=50cm
load("ML_mortality_other.more50.Rdata")
spother.more50.m.rate = ML_mortality_other.more50$ML_mortality
#mortality rate of all individuals that are not part of our 39species that have dbh<50cm
load("ML_mortality_other.less50.Rdata")
spother.less50.m.rate = ML_mortality_other.less50$ML_mortality

### Extra biomass: this is becuase some plots do not have the same size structure and trees are not measured starting at the same dbh (missing small size classes)
# we extrapolated from the plots that measure small size classe weighted by the area of the plots missing
# for young/secondary the reference is: young quadrats of plot50ha
dataset = baseline
biomass.39sp.less50.young.df = data.frame(species = species.to.use,
                                          biomass.39sp.10_25mm = 0,
                                          biomass.39sp.25_100mm = 0,
                                          biomass.39sp.100_200mm = 0)

funk.extra.young = function(x){
  use = dataset[dataset$sp == x & dataset$plot == "plot50ha",] #species x in the young quadrats of the 50ha plot
  use2 = sum(use[use$dbhmm.1.3>=10 & use$dbhmm.1.3<25, "AGB"]) #10-25mm
  use3 = sum(use[use$dbhmm.1.3>=25 & use$dbhmm.1.3<100, "AGB"]) #25-100mm
  use4 = sum(use[use$dbhmm.1.3>=100 & use$dbhmm.1.3<200, "AGB"]) #100-200mm
  biomass.39sp.less50.young.df[biomass.39sp.less50.young.df$species == x, "biomass.39sp.10_25mm"] <<- use2
  biomass.39sp.less50.young.df[biomass.39sp.less50.young.df$species == x, "biomass.39sp.25_100mm"] <<- use3
  biomass.39sp.less50.young.df[biomass.39sp.less50.young.df$species == x, "biomass.39sp.100_200mm"] <<- use4
}
invisible(lapply(species.to.use, funk.extra.young))

ref.50ha.young = dataset[dataset$plot == "plot50ha",]
other.50ha.young = ref.50ha.young[-(which(ref.50ha.young$sp %in% species.to.use)),]
biomass.other.sp.10_25mm = sum(other.50ha.young[other.50ha.young$dbhmm.1.3>=10 & other.50ha.young$dbhmm.1.3<25, "AGB"])
biomass.other.sp.25_100mm = sum(other.50ha.young[other.50ha.young$dbhmm.1.3>=25 & other.50ha.young$dbhmm.1.3<100, "AGB"])
biomass.other.sp.100_200mm = sum(other.50ha.young[other.50ha.young$dbhmm.1.3>=100 & other.50ha.young$dbhmm.1.3<200, "AGB"])

#area missing each size class: 10-25mm -> plot25ha(27.04), 10ha.young(6), AVA(4), Knight(1.5)
#                              25-100mm -> plot25ha, 10ha.young, AVA
#                              100-200mm -> 10ha.young only (10ha starts measuring at 200mm)
biomass.other.sp.extra.less50.young = biomass.other.sp.10_25mm * (27.04+6+4+1.5)/1.92 +
  biomass.other.sp.25_100mm * (27.04+6+4)/1.92 +
  biomass.other.sp.100_200mm * (6)/1.92

biomass.39sp.less50.young.df$biomass.w.39sp.10_25mm = biomass.39sp.less50.young.df$biomass.39sp.10_25mm * (27.04+6+4+1.5)/1.92
biomass.39sp.less50.young.df$biomass.w.39sp.25_100mm = biomass.39sp.less50.young.df$biomass.39sp.25_100mm * (27.04+6+4)/1.92
biomass.39sp.less50.young.df$biomass.w.39sp.100_200mm = biomass.39sp.less50.young.df$biomass.39sp.100_200mm * (6)/1.92

biomass.39sp.extra.less50.young = sum(biomass.39sp.less50.young.df[,c("biomass.w.39sp.10_25mm","biomass.w.39sp.25_100mm","biomass.w.39sp.100_200mm")])

hectares.used.young = 40.46 #plot25ha(27.04)+10ha.young(6)+AVA(4)+Knight(1.5)+50ha.young(1.92)

#for old-growth/mature the reference: Gigante + old quadrats of plot50ha
dataset = old
biomass.39sp.less50.mature.df = data.frame(species = species.to.use, biomass.39sp.10_200mm = 0)
ref.mature = dataset[dataset$plot == "Gigante" | dataset$plot == "plot50ha",]

funk.extra.mature = function(x){
  use = ref.mature[ref.mature$sp == x,]
  use2 = sum(na.omit(use[use$dbhmm.1.3>=10 & use$dbhmm.1.3<200, "AGB"])) #10-200mm
  biomass.39sp.less50.mature.df[biomass.39sp.less50.mature.df$species == x, "biomass.39sp.10_200mm"] <<- use2
}
invisible(lapply(species.to.use, funk.extra.mature))

other.mature.ref = ref.mature[-(which(ref.mature$sp %in% species.to.use)),]
biomass.other.sp.10_200mm = sum(na.omit(other.mature.ref[other.mature.ref$dbhmm.1.3>=10 & other.mature.ref$dbhmm.1.3<200, "AGB"]))

#area missing the 10-200mm class: 10ha.old(4), Pearson(6), Zetek(6), Drayton(6)
#reference area: old quadrats of 50ha (47.64) + Gigante (36.24), the mature plots that do measure it
biomass.other.sp.extra.less50.mature = biomass.other.sp.10_200mm * (4+6+6+6)/(47.64+36.24)

biomass.39sp.less50.mature.df$biomass.w.39sp.10_200mm = biomass.39sp.less50.mature.df$biomass.39sp.10_200mm * (4+6+6+6)/(47.64+36.24)
biomass.39sp.extra.less50.mature = sum(biomass.39sp.less50.mature.df$biomass.w.39sp.10_200mm)

hectares.used.mature = 105.88 #50ha.old(47.64)+10ha.old(4)+Pearson(6)+Zetek(6)+Drayton(6)+Gigante(36.24)
hectares.used = hectares.used.young + hectares.used.mature #146.34

#remove the species that we have mortality rates of giants from the whole dataset to use them with their own mortality rates
sp39 = temp.all.plots.all.sp[temp.all.plots.all.sp$sp%in%species.to.use,]; length(unique(sp39$sp)); nrow(sp39) #39 species
spother = temp.all.plots.all.sp[!temp.all.plots.all.sp$sp%in%species.to.use,]; length(unique(spother$sp)) + length(unique(sp39$sp)); nrow(spother)
spother.more50 = spother[spother$dbhmm.1.3 >= 500,]; nrow(spother.more50); length(unique(spother.more50$sp))
spother.less50 = spother[spother$dbhmm.1.3 < 500,]; nrow(spother.less50)

#divide the dataset in below 50cm and above 50cm dbh for the species with giant mortality rates
sp39$dbhmm.1.3 = as.numeric(sp39$dbhmm.1.3)
sp39.more.50 = sp39[sp39$dbhmm.1.3>=500,];length(unique(sp39.more.50$sp))
sp39.less.50 = sp39[sp39$dbhmm.1.3<500,];length(unique(sp39.less.50$sp))
sp39.more.50 = sp39.more.50[sp39.more.50$sp %in% unique(sp39.less.50$sp),];length(unique(sp39.more.50$sp))

#assigning community mortality rate to the trees lower than 50cm from our species
sp39.less.50$ML_mortality= sp39.less50.m.rate

#add mortality rates for giant trees
funk = function(x){
  use = match(toString(x), data$sp)
  mr = data[use, "ML_mortality"]
  return(mr)
}
sp39.more.50$ML_mortality = unlist(lapply(sp39.more.50$sp, funk))

#WOODY RESIDENCE TIMES of the species on the whole community

#above-ground biomass for the whole community/ha and the extra biomass terms
agb.all = (sum(na.omit(temp.all.plots.all.sp$AGB)) + biomass.39sp.extra.less50.young + biomass.39sp.extra.less50.mature +
             biomass.other.sp.extra.less50.young + biomass.other.sp.extra.less50.mature)/hectares.used; agb.all
agb.other = (sum(na.omit(spother$AGB)) + biomass.other.sp.extra.less50.young + biomass.other.sp.extra.less50.mature)/hectares.used; agb.other
agb.other.more50 = sum(na.omit(spother.more50$AGB))/hectares.used; agb.other.more50
agb.other.less50 = (sum(na.omit(spother.less50$AGB)) + biomass.other.sp.extra.less50.young + biomass.other.sp.extra.less50.mature)/hectares.used; agb.other.less50

#species 39 giant over total biomass
(sum(na.omit(sp39.more.50$AGB))/hectares.used)/agb.all*100

#all species giant over total biomass
temp.all.plots.all.sp_subset_50 = temp.all.plots.all.sp[temp.all.plots.all.sp$dbhmm.1.3>=500,]
total.agb_50 = sum(na.omit(temp.all.plots.all.sp_subset_50$AGB))/hectares.used 
total.carbon_50 = 0.47*total.agb_50;total.carbon_50 

temp.all.plots.all.sp_subset_less50 = temp.all.plots.all.sp[temp.all.plots.all.sp$dbhmm.1.3<500,]
total.agb_less50 = sum(na.omit(temp.all.plots.all.sp_subset_less50$AGB))/hectares.used 

#species 39 giants over biomass of giants
(sum(na.omit(sp39.more.50$AGB))/hectares.used)/total.agb_50*100

### *Figure S2b: bar plot of biomass stored in different size classes
df <- data.frame(
  size_class = c("1-50", ">50"),
  carbon = c(total.agb_less50, total.agb_50)
)
df <- df %>%
  mutate(size_class = recode(size_class, ">50" = "≥50"))
df$size_class <- factor(df$size_class, levels = c("1-50", "≥50"))

biomass_stock = ggplot(df, aes(x = size_class, y = carbon)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(x = expression(bold("Size class (cm)"))) +
  labs(y = expression(bold("Aboveground biomass" ~('Mg'~ha^{-1})))) +
  ylim(0, 100) + theme_basis

#woody residence time of all the community that includes all trees, with their corresponding mrates
sp39.wrt.less50 = c()
sp39.wrt.more50 = c()

#species with at least one individual >=50cm dbh - the ones we run the giant-mortality experiment on
species.to.use = unique(sp39.more.50$sp)
funk.wrt = function(x, y){
  
  use1 = sp39.less.50[sp39.less.50$sp == x,] #species x and individuals with dbh<50cm from our species
  use2 = y[y$sp == x,] #species x and individuals with dbh>50cm from our species
  #extra biomass for species x's smaller size classes, extrapolated above (young + mature)
  extra.young.x = biomass.39sp.less50.young.df[biomass.39sp.less50.young.df$species == x,
                                               c("biomass.w.39sp.10_25mm","biomass.w.39sp.25_100mm","biomass.w.39sp.100_200mm")]
  extra.mature.x = biomass.39sp.less50.mature.df[biomass.39sp.less50.mature.df$species == x, "biomass.w.39sp.10_200mm"]
  extra.x = sum(unlist(extra.young.x)) + sum(extra.mature.x)
  agb.species.less50 = (sum(na.omit(use1$AGB)) + extra.x)/hectares.used #Mg/ha
  agb.species.more50 = sum(na.omit(use2$AGB))/hectares.used #Mg/ha
  
  #woody residence time = sum(biomass of a speciesi)/ (total biomass of a species * mrate)
  den.wrt.less50 = (agb.species.less50 * unique(use1$ML_mortality)) #community mortality rates for trees <50cm of our species(~3%)
  den.wrt.more50 = (agb.species.more50 * unique(use2$ML_mortality)) #sp-specific giant mortality rate
  
  sp39.wrt.less50[x] <<- den.wrt.less50
  sp39.wrt.more50[x] <<- den.wrt.more50
}
t = lapply(species.to.use,funk.wrt, y = sp39.more.50) #fun(x,  funk, y)

#wood residence time for the community, all individuals
#woody residence time =  (total biomass of all species)/ (sum(biomass of a speciesi) * (1 / (mortality rate of species)))
#baseline (observed mortality rates) - kept under its own name so the loop below does not overwrite it
wrt.normal = (agb.all)/((sum(sp39.wrt.less50)+sum(sp39.wrt.more50))+(agb.other.less50*spother.less50.m.rate)+(agb.other.more50*spother.more50.m.rate))

#table for values
wrt.table = data.frame("species" = species.to.use, "wrt" = NA)

#now we change the LONG-LIVED GIANT mortality rate to mean community rate
#change for each species at the time
for(i in species.to.use){
  use = sp39.more.50
  use[use$sp == i, "ML_mortality"] = all.comm.giant.m.rate
  sp39.wrt.less50 = c()
  sp39.wrt.more50 = c()
  t = lapply(species.to.use, funk.wrt, y = use) 
  
  wrt = (agb.all)/((sum(sp39.wrt.less50)+sum(sp39.wrt.more50))+(agb.other.less50*spother.less50.m.rate)+(agb.other.more50*spother.more50.m.rate))
  wrt.table[wrt.table$species == i,"wrt"] = as.numeric(wrt)
}

wrt.table[nrow(wrt.table)+1,] = c("normal", wrt.normal)

#add mortality rate for graph
funk = function(x){
  use = match(toString(x), data$sp)
  mr = data[use, "ML_mortality"]
  return(mr)
}
wrt.table$ML_mortality = unlist(lapply(wrt.table$species, funk))
wrt.table$wrt = as.numeric(wrt.table$wrt)

#add taxon for graph
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk = function(x){
  use = match(toString(x), bci.spptable$sp)
  mr = bci.spptable[use, "taxon"]
  return(mr)
}
wrt.table$taxon = unlist(lapply(wrt.table$species, funk))
wrt.table[wrt.table$species == "normal", "taxon"] = "Normal"

#plots
#add a taxon 2 to get first letter of genus and species name
wrt.table$taxon2 = sub("^(\\w)\\w*\\s(\\w+)$", "\\1. \\2", wrt.table$taxon)

#if woody residence time (AWRT) is changed while woody
#productivity (AWP) stays the same, then AGB will change as follows:
#(new AGB) = (old AGB) *(new AWRT) / (old AWRT).
#an X% decrease in AWRT will lead to an X% decrease in AGB, all else equal.
newAGB = agb.all*(wrt.table$wrt)/wrt.table[wrt.table$species == "normal","wrt"]

#if the focal species dies at a giant tree community wide mortality rate I calculated a new WRT, then
#this is the new AWRT, then calculated the newAGB. Then calculated the percentage change.
perc_AGB =  100 - (100*newAGB/agb.all)
wrt.table$perc_AGB = perc_AGB

#reorder wrt.table by mortality rates
wrt.table = wrt.table[order(wrt.table$ML_mortality),]
sum(wrt.table$perc_AGB [1:5])

### *Figure 2: expected change in agb when each giant tree starts dying at a mean giant tree rate
agb.plot.perc = ggplot(wrt.table[1:nrow(wrt.table)-1,], aes(x = -perc_AGB, y = reorder(taxon2, ML_mortality)))+
  theme(axis.text.y = element_text(face = "italic")) +
  geom_col(fill = "royalblue", width = 0.8, alpha = 0.5) +
  theme(legend.position = "none") + scale_x_continuous(limits =c(-5,2.5))+
  theme_basis + labs(x = expression(bold("Expected change in AGB (%)")), y = "Species")+
  geom_vline(xintercept = 0, linewidth = 0.2, color = "darkred")
