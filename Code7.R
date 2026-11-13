########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Biomass mortality rates in old-growth and secondary forest weighting by hectare-years
########################################################################################################################

### this code calculates the biomass mortality rate for secondary forests and compare it with
## biomass mortality rate of old-growth forests. 
# Biomass mortality rate = equation 4 in Kohyama et al. (2019), and 
# Biomass mortality rate = (log(AGB0/AGBS0)/T)*100
# AGB0 represents the total biomass at time 0
# AGBS0 represents the initial biomass of survivors at time 0
#T represents the length of the census interval (years). 

#clear workspace
rm(list=ls()); library(ggplot2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12))+ #, face ="bold"
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                       simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                       !is.na(simple_combined$dbhmm.1.3),]
#SECONDARY FORESTS
# we cannot use Dent's data since it does not have any demographics
#lang.k, plot25ha, AVA, part of 10ha 
plotlang = simple_combined_sp[simple_combined_sp$plot == "plotKnight",] 
plot25ha = simple_combined_sp[simple_combined_sp$plot == "plot25ha",] 
plotava = simple_combined_sp[simple_combined_sp$plot == "plotAVA",] 

#the 10ha plot has one part that is mature forest (remove it and use it below)
plot10ha.temp = simple_combined_sp[simple_combined_sp$plot == "plot10ha",]

#around 420 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 420m there is mature forest, to the right of 420m it is secondary.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and cuts the left size of the 10ha plot, then to the left of 780m there is young forest, to the right of 780m it has both secondary/mature.
#around 780 meter on x-axis of the 50haplot starts the mixed/young forest, and above 50m in the y-axis it is young forest. 
load("all2.rdata")
#View(all2) Use all2 from Maekem 100yearsBCI to select the quadrats
plot10ha.1 = plot10ha.temp[plot10ha.temp$gx>250 & plot10ha.temp$gx<750,] #young 10ha area 1 #5ha
plot10ha.2 = plot10ha.temp[plot10ha.temp$gx>800 & plot10ha.temp$gy>50,] #young 10ha area 2 #1ha
plot10ha.old = plot10ha.temp[plot10ha.temp$gx<200,] #old 10haplot ~4ha
plot10ha.young = rbind(plot10ha.1, plot10ha.2)  #6ha

##where is the young forest areas in the 50ha plot?
#loads the topographic types for each quadrat
harms = read.table("harms.txt",h=T)
harms$quadrat = paste(harms$GX20/20, (harms$GY20/20), sep = "") #code to name the quadrats on the 50haplot
young.forest.50ha = harms[harms$Habitat == "young",] #48 quadrats 
mixed.forest.50ha = c("2123","2222","2422","2522","2622","2722","3320","3520","3620","3722","3824") #this is mixed in vicinity of the young forest #11quads

#Add young forest data to the secondary forest data from above (young 50ha is 1.92ha = 48quadrats*400m2)
plot50ha = simple_combined_sp[simple_combined_sp$plot == "plot50ha",]
plot50ha.young = plot50ha[which(plot50ha$quadrat %in% young.forest.50ha$quadrat),] #63,247 rows
plot50ha.old = plot50ha[-which(plot50ha$quadrat %in% young.forest.50ha$quadrat),] #1,516,665 rows
plot50ha.old = plot50ha.old[-which(plot50ha.old$quadrat %in% mixed.forest.50ha),] #1,502,438 rows

#bind the 4 young/secondary datasets: Lang-Knight- Plot25ha, PlotAVA, Plot10ha
baseline = rbind(plotlang, plot25ha, plotava, plot10ha.young, plot50ha.young) #75,354
unique(baseline$plot)

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

#adding wood density to the dataframe for later biomass calculations
jw_data = read.table("WD_species.txt", h = T)
jw_data = data.frame(jw_data)
jw_data$SP. = tolower(jw_data$sp6)

#adding wd in g/cm3
funk_wd = function(species){ 
  use = match(toString(species), jw_data$SP.)
  wd = jw_data[use, "wd100.mean"]
  return(wd)
}
baseline$wd = unlist(lapply(baseline$sp, funk_wd))
baseline = as.data.frame(baseline)

#extras wd from ctfs data 
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

baseline[which(baseline$sp %in% unique.with.nas), "wd"] =
  mean(na.omit(baseline$wd))

#calculate their AGB
#when height is unknown - Chave et al., 2015 (Equation 7) for moist tropical forest
#D in "cm"; wsg in "g/cm3"
#Kg/1000 to convert to Mg
#dbhmm.1.3 is already corrected with taper equations when the data was cleaned

#For palm species -> ln(ABG) = -3.3488 + 2.7484*ln(dbh) * CF (in Kilograms)
#CF correction factor by log-transforming the data see Chave et al., 2005 (CF = exp((RSE^2)/2) = 1.188714)
palms = baseline[baseline$family == "Arecaceae",]
palms$AGB = (exp(-3.3488 + 2.7484*log(palms$dbhmm.1.3/10))/1000) * (1.188714)

trees = baseline[which(!baseline$family == "Arecaceae"),]
E = 0.051763982 #for BCI
trees$AGB = (exp(-1.803-0.976*E+(0.976*log(trees$wd))+
                   (2.673*log(trees$dbhmm.1.3/10))-
                   (0.0299*(log(trees$dbh/10))^2))/1000)

baseline = rbind(trees, palms)
length(unique(baseline$sp)) # 268 species
#Basal Area
baseline$ba = (pi/4)*(baseline$dbhmm.1.3)^2
#Convert to m2 instead of mm2
baseline$ba = baseline$ba/1000000

# First we calculate the extra biomass (computed once from the real data, not resampled
# only the main per plot tree data below is quadrat-bootstrapped)
#plot50ha
tmp.50ha = baseline[which(baseline$plot == "plot50ha"),] #baseline is the young forests area
ha.50ha.young = 1.92

# extra biomass loss for small size classes
##10-25mm 25ha(27.04), 10hayoung(6ha),  ava(4), knight(1.5ha)
#25-100mm 25ha, 10hayoung, ava
#100-200mm 10hayoung
#10-25mm
tmp.50ha.10_25 = tmp.50ha[tmp.50ha$dbhmm.1.3>=10 & tmp.50ha$dbhmm.1.3<25,]; unique(tmp.50ha.10_25$openyear); unique(tmp.50ha.10_25$endyear)

funk2B0 = function(x){
  census = tmp.50ha.10_25[tmp.50ha.10_25$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.10_25 = unlist(lapply(unique(tmp.50ha.10_25$openyear), funk2B0)) #0.85 #biomass/hectare

funk2BS0 = function(x){
  census = tmp.50ha.10_25[tmp.50ha.10_25$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.10_25 = unlist(lapply(unique(tmp.50ha.10_25$openyear), funk2BS0));BS0.50ha.10_25 #biomass/hectare

#25-100mm
tmp.50ha.25_100 = tmp.50ha[tmp.50ha$dbhmm.1.3>=25 & tmp.50ha$dbhmm.1.3<100,]; unique(tmp.50ha.25_100$openyear); unique(tmp.50ha.25_100$endyear)

funk3B0 = function(x){
  census = tmp.50ha.25_100[tmp.50ha.25_100$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.25_100 = unlist(lapply(unique(tmp.50ha.25_100$openyear), funk3B0)) #0.85

funk3BS0 = function(x){
  census = tmp.50ha.25_100[tmp.50ha.25_100$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.25_100 = unlist(lapply(unique(tmp.50ha.25_100$openyear), funk3BS0)); BS0.50ha.25_100 #206

#100-200mm
tmp.50ha.100_200 = tmp.50ha[tmp.50ha$dbhmm.1.3>=100 & tmp.50ha$dbhmm.1.3<200,]; unique(tmp.50ha.100_200$openyear); unique(tmp.50ha.100_200$endyear)

funk4B0 = function(x){
  census = tmp.50ha.100_200[tmp.50ha.100_200$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.100_200 = unlist(lapply(unique(tmp.50ha.100_200$openyear), funk4B0))

funk4BS0 = function(x){
  census = tmp.50ha.100_200[tmp.50ha.100_200$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.100_200 = unlist(lapply(unique(tmp.50ha.100_200$openyear), funk4BS0)); BS0.50ha.100_200 
al = B0.50ha.100_200

safe.mean = function(x) if (length(x) == 0) 0 else mean(x)
B0.extra.y.10_25 = safe.mean(B0.50ha.10_25)
BS0.extra.y.10_25 = safe.mean(BS0.50ha.10_25)
B0.extra.y.25_100 = safe.mean(B0.50ha.25_100)
BS0.extra.y.25_100 = safe.mean(BS0.50ha.25_100)
B0.extra.y.100_200 = safe.mean(B0.50ha.100_200)
BS0.extra.y.100_200 = safe.mean(BS0.50ha.100_200)

#compute one biomass-loss-rate row per census interval present for a given plot subset of the
#current bootstrap draw. 
census_interval_rates = function(plot_df, plot_name, ha, extra_B0 = 0, extra_BS0 = 0){
  results = data.frame()
  for(x in sort(unique(plot_df$openyear))){
    use = plot_df[plot_df$openyear == x,]
    
    interval_length = mean(use$endyear - use$openyear, na.rm = TRUE)
    B0 = sum(na.omit(use$AGB))/ha + extra_B0
    BS0 = sum(na.omit(use[use$mortality == 0,"AGB"]))/ha + extra_BS0
    
    results = bind_rows(results, data.frame(plot = plot_name, census = x, interval_length = interval_length,
                                            B0 = B0, BS0 = BS0, ha = ha,
                                            rate = log(B0/BS0)/interval_length,
                                            weight = interval_length*ha))
  }
  results
}

#weighted mean rate (%) across all (plot, census interval) rows: weight = interval length x plot
#area (ha), so larger/longer-monitored plots contribute more to the mean regardless of how much
#biomass they happened to start with
weighted_loss_rate = function(all_rates){
  ok = !is.na(all_rates$rate) & !is.na(all_rates$weight) & all_rates$weight > 0
  sum(all_rates$rate[ok]*all_rates$weight[ok]) / sum(all_rates$weight[ok]) * 100
}
#collapse once to one row per (plot, plot_treeID, openyear), keeping the row with the highest
#endcensus 
baseline.collapsed = baseline %>%
  group_by(plot, plot_treeID, openyear) %>% arrange(desc(endcensus)) %>% slice(1) %>%
  ungroup() %>% as.data.frame()

#pre-split row indices by quadrat once, so each bootstrap draw is a fast lookup
baseline.quads = unique(baseline.collapsed$quadrat_plot)
baseline.quad.idx = split(seq_len(nrow(baseline.collapsed)), baseline.collapsed$quadrat_plot)
biomass.loss.rate.y.boot = rep(0, 1000) #for CIs

for(k in 1:1000){
  if (k == 1){
    dataset = baseline.collapsed
  }else{
    #resample whole quadrats with replacement (not individual trees)
    boot.quads = sample(baseline.quads, size = length(baseline.quads), replace = TRUE)
    dataset = baseline.collapsed[unlist(baseline.quad.idx[boot.quads]),]
  }
  ha.50ha.young = 1.92
  ha_lang = 1.5
  ha_25ha = 27.04
  ha_AVA = 4
  ha_10ha = 6
  
  #each plot is missing a different subset of the small size classes, different minimum-dbh
  #so the correction added to each plot is not the same total:
  #10-25mm:  plot25ha, plot10ha.young, plotAVA, plotKnight (all four)
  #25-100mm: plot25ha, plot10ha.young, plotAVA (not plotKnight)
  #100-200mm: plot10ha.young only
  extra_B0.knight = B0.extra.y.10_25
  extra_BS0.knight = BS0.extra.y.10_25
  extra_B0.25ha = B0.extra.y.10_25 + B0.extra.y.25_100
  extra_BS0.25ha = BS0.extra.y.10_25 + BS0.extra.y.25_100
  extra_B0.AVA = B0.extra.y.10_25 + B0.extra.y.25_100
  extra_BS0.AVA = BS0.extra.y.10_25 + BS0.extra.y.25_100
  extra_B0.10ha = B0.extra.y.10_25 + B0.extra.y.25_100 + B0.extra.y.100_200
  extra_BS0.10ha = BS0.extra.y.10_25 + BS0.extra.y.25_100 + BS0.extra.y.100_200
  
  #plot50ha:no extra small-size-class correction
  #drop the 2015 census (few trees recorded that year)
  tmp.50ha = dataset[dataset$plot == "plot50ha" & floor(dataset$openyear) != 2015,]
  rates.50ha.y = census_interval_rates(tmp.50ha, "plot50ha.young", ha.50ha.young)
  
  #small-size-class correction (see above) estimated from plot50ha above
  rates.lang = census_interval_rates(dataset[dataset$plot == "plotKnight",], "plotKnight", ha_lang, extra_B0.knight, extra_BS0.knight)
  rates.25ha = census_interval_rates(dataset[dataset$plot == "plot25ha",], "plot25ha", ha_25ha, extra_B0.25ha, extra_BS0.25ha)
  rates.AVA = census_interval_rates(dataset[dataset$plot == "plotAVA",], "plotAVA", ha_AVA, extra_B0.AVA, extra_BS0.AVA)
  rates.10ha = census_interval_rates(dataset[dataset$plot == "plot10ha",], "plot10ha.young", ha_10ha, extra_B0.10ha, extra_BS0.10ha)
  
  all.rates.y = bind_rows(rates.50ha.y, rates.lang, rates.25ha, rates.AVA, rates.10ha)
  
  #weighted biomass loss rate (%) across every plot x census-interval combination
  biomass.loss.rate.y = weighted_loss_rate(all.rates.y); biomass.loss.rate.y
  
  biomass.loss.rate.y.boot[k] = biomass.loss.rate.y
  if (k %% 100 == 0) print(k)
}

#OLD GROWTH FORESTS
#add data from Pearson, Zetek and Drayton as mature forests
plotPearson = simple_combined_sp[simple_combined_sp$plot == "plotPearson",] #866 rows
plotZetek = simple_combined_sp[simple_combined_sp$plot == "plotZetek",] #866 rows
plotDrayton = simple_combined_sp[simple_combined_sp$plot == "plotDrayton",] #843 rows
#gigante peninsula plot >600 years
plotgigante = simple_combined_sp[simple_combined_sp$plot == "Gigante",];nrow(plotgigante) #12964 rows
plotgigante = plotgigante[plotgigante$gigante_plot >= 0,] #trees that are in plots that measure trees >=1cm
plotgigante = plotgigante[plotgigante$dbhmm.1.3 > 0,] #trees that are in plots that measure trees >=1cm

old = rbind(plot10ha.old, plot50ha.old, plotPearson, plotZetek, plotDrayton, plotgigante) #1,505,317

## Add taxon to the data
old$taxon = unlist(lapply(old$sp, funk_taxon))
old$family = unlist(lapply(old$sp, funk_family))

#adding wood density in g/cm3 to the dataframe for later biomass calculations
old$wd = unlist(lapply(old$sp, funk_wd))
old = as.data.frame(old)

#extras wood from ctfs data when the species did not have reported wood density in the previous dataset
q = unlist(lapply(unique.with.nas, funk_ctfs))

#for species with No wood density value: the mean of all other species in the site was used
with.nas = old[is.na(old$wd),]
unique.with.nas = unique(with.nas$sp)

old[which(old$sp %in% unique.with.nas), "wd"] =
  mean(na.omit(old$wd))

#calculate their AGB
#when height is unknown - Chave et al., 2015 (Equation 7) for moist tropical forest
#D in "cm"; wsg in "g/cm3"
#Kg/1000 to convert to Mg
#dbhmm.1.3 is already corrected with taper equations when the data was cleaned

#For palm species -> ln(ABG) = -3.3488 + 2.7484*ln(dbh) * CF (in Kilograms)
#CF correction factor by log-transforming the data see Chave et al., 2005 (CF = exp((RSE^2)/2) = 1.188714)
palms = old[old$family == "Arecaceae",]
palms$AGB = (exp(-3.3488 + 2.7484*log(palms$dbhmm.1.3/10))/1000) * (1.188714)

trees = old[which(!old$family == "Arecaceae"),]
E = 0.051763982 #for BCI
trees$AGB = (exp(-1.803-0.976*E+(0.976*log(trees$wd))+
                   (2.673*log(trees$dbhmm.1.3/10))-
                   (0.0299*(log(trees$dbhmm.1.3/10))^2))/1000)

old = rbind(trees, palms)
length(unique(old$sp)) # 335 species
#Basal Area
old$ba = (pi/4)*(old$dbhmm.1.3)^2
#Convert to m2 instead of mm2
old$ba = old$ba/1000000

# First we calculate the extra biomass to add it to the calculations (computed once from the
# real data, not resampled - only the main per-plot tree data below is quadrat-bootstrapped)
#plot50ha - this plot has many census intervals
tmp.50ha = old[which(old$plot == "plot50ha"),] #old is the old growth forests area

# extra biomass loss for small size classes
##10-25mm 25ha(27.04), 10hayoung(6ha),  ava(4), knight(1.5ha)
#25-100mm 25ha, 10hayoung, ava
#100-200mm 10hayoung
#10-25mm
tmp.50ha.10_25 = tmp.50ha[tmp.50ha$dbhmm.1.3>=10 & tmp.50ha$dbhmm.1.3<25,]; unique(tmp.50ha.10_25$openyear); unique(tmp.50ha.10_25$endyear)

funk2B0 = function(x){
  census = tmp.50ha.10_25[tmp.50ha.10_25$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.10_25 = unlist(lapply(unique(tmp.50ha.10_25$openyear), funk2B0)) #0.85 #biomass/hectare

funk2BS0 = function(x){
  census = tmp.50ha.10_25[tmp.50ha.10_25$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.10_25 = unlist(lapply(unique(tmp.50ha.10_25$openyear), funk2BS0));BS0.50ha.10_25 #biomass/hectare

#25-100mm
tmp.50ha.25_100 = tmp.50ha[tmp.50ha$dbhmm.1.3>=25 & tmp.50ha$dbhmm.1.3<100,]; unique(tmp.50ha.25_100$openyear); unique(tmp.50ha.25_100$endyear)

funk3B0 = function(x){
  census = tmp.50ha.25_100[tmp.50ha.25_100$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.25_100 = unlist(lapply(unique(tmp.50ha.25_100$openyear), funk3B0)) #0.85
if(length(B0.50ha.25_100)>7){B0.50ha.25_100 = B0.50ha.25_100[1:7]}

funk3BS0 = function(x){
  census = tmp.50ha.25_100[tmp.50ha.25_100$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.25_100 = unlist(lapply(unique(tmp.50ha.25_100$openyear), funk3BS0)); BS0.50ha.25_100 #206
if(length(BS0.50ha.25_100)>7){BS0.50ha.25_100 = BS0.50ha.25_100[1:7]}

#100-200mm
tmp.50ha.100_200 = tmp.50ha[tmp.50ha$dbhmm.1.3>=100 & tmp.50ha$dbhmm.1.3<200,]; unique(tmp.50ha.100_200$openyear); unique(tmp.50ha.100_200$endyear)

funk4B0 = function(x){
  census = tmp.50ha.100_200[tmp.50ha.100_200$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use.total.biomass = sum(na.omit(use$AGB))/ha.50ha.young
  return(use.total.biomass)
}
B0.50ha.100_200 = unlist(lapply(unique(tmp.50ha.100_200$openyear), funk4B0)); B0.50ha.100_200
if(length(B0.50ha.100_200)>7){B0.50ha.100_200 = B0.50ha.100_200[1:7]}

funk4BS0 = function(x){
  census = tmp.50ha.100_200[tmp.50ha.100_200$openyear == x,]
  use =  census %>% #this includes mortality = 1,0 (all trees in that census)
    group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)
  use2 = use[use$mortality == 0,]
  use.total.biomass = sum(na.omit(use2$AGB))/ha.50ha.young
  return(use.total.biomass)
}
BS0.50ha.100_200 = unlist(lapply(unique(tmp.50ha.100_200$openyear), funk4BS0)); BS0.50ha.100_200 
if(length(BS0.50ha.100_200)>7){BS0.50ha.100_200 = BS0.50ha.100_200[1:7]}

B0.extra.o.10_25 = safe.mean(B0.50ha.10_25)
BS0.extra.o.10_25 = safe.mean(BS0.50ha.10_25)
B0.extra.o.25_100 = safe.mean(B0.50ha.25_100)
BS0.extra.o.25_100 = safe.mean(BS0.50ha.25_100)
B0.extra.o.100_200 = safe.mean(B0.50ha.100_200)
BS0.extra.o.100_200 = safe.mean(BS0.50ha.100_200)

#collapse once to one row per (plot, plot_treeID, openyear)
old.collapsed = old %>%
  group_by(plot, plot_treeID, openyear) %>% arrange(desc(endcensus)) %>% slice(1) %>%
  ungroup() %>% as.data.frame()

#pre split row indices by quadrat once, so each bootstrap draw is a fast lookup
old.quads = unique(old.collapsed$quadrat_plot)
old.quad.idx = split(seq_len(nrow(old.collapsed)), old.collapsed$quadrat_plot)
biomass.loss.rate.o.boot = rep(0, 1000) #for CIs

for(m in 1:1000){
  if (m == 1){
    dataset = old.collapsed
  }else{
    #resample whole quadrats with replacement (not individual trees)
    boot.quads = sample(old.quads, size = length(old.quads), replace = TRUE)
    dataset = old.collapsed[unlist(old.quad.idx[boot.quads]),]
  }
  ha_pearson = 6
  ha_zetek = 6
  ha_drayton = 6
  ha_gigante = 36.24 #already removed the area with not small size class
  ha_10ha_old = 4
  ha.50ha.o = 47.64
  
  extra_B0.o = B0.extra.o.10_25 + B0.extra.o.25_100 + B0.extra.o.100_200
  extra_BS0.o = BS0.extra.o.10_25 + BS0.extra.o.25_100 + BS0.extra.o.100_200
  
  rates.pearson = census_interval_rates(dataset[dataset$plot == "plotPearson" & dataset$openyear == "2004",],
                                        "plotPearson", ha_pearson, extra_B0.o, extra_BS0.o)
  rates.zetek = census_interval_rates(dataset[dataset$plot == "plotZetek",], "plotZetek", ha_zetek, extra_B0.o, extra_BS0.o)
  rates.drayton = census_interval_rates(dataset[dataset$plot == "plotDrayton",], "plotDrayton", ha_drayton, extra_B0.o, extra_BS0.o)
  rates.10ha.old = census_interval_rates(dataset[dataset$plot == "plot10ha",], "plot10ha.old", ha_10ha_old, extra_B0.o, extra_BS0.o)
  
  tmp.gigante = dataset[dataset$plot == "Gigante" & dataset$gigante_plot >= 0,]
  rates.gigante = census_interval_rates(tmp.gigante, "Gigante", ha_gigante)

  tmp.50ha.o = dataset[dataset$plot == "plot50ha" & floor(dataset$openyear) != 2015,]
  rates.50ha.o = census_interval_rates(tmp.50ha.o, "plot50ha.old", ha.50ha.o)
  
  all.rates.o = bind_rows(rates.pearson, rates.zetek, rates.drayton, rates.10ha.old, rates.gigante, rates.50ha.o)
  
  #weighted biomass loss rate (%) across every plot x census-interval combination
  biomass.loss.rate.o = weighted_loss_rate(all.rates.o); biomass.loss.rate.o
  
  biomass.loss.rate.o.boot[m] = biomass.loss.rate.o
  if(m %% 100 == 0) print(m)
}

#point estimates (iteration 1 = real, unresampled data)
biomass.loss.rate.y = biomass.loss.rate.y.boot[1]; biomass.loss.rate.y #1.72%
biomass.loss.rate.o = biomass.loss.rate.o.boot[1]; biomass.loss.rate.o #1.77%

#95% CI from the quadrat cluster bootstrap (999 resampled draws, excluding the real-data draw)
biomass.loss.rate.y.CI = quantile(biomass.loss.rate.y.boot[-1], c(0.025, 0.975), na.rm = TRUE); biomass.loss.rate.y.CI
biomass.loss.rate.o.CI = quantile(biomass.loss.rate.o.boot[-1], c(0.025, 0.975), na.rm = TRUE); biomass.loss.rate.o.CI

