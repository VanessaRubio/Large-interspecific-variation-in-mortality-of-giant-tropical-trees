########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Are the hectares sampled in the old-growth forest influencing the observed max diameters?
## Produces: Figure S6
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggplot2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12)) #, face ="bold"

setwd("/Users/vanessarubio/Library/CloudStorage/GoogleDrive-rubiov@caryinstitute.org/My Drive/Cary/A.Longevity/")

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
#select only the trees with dbh in 1968
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
#View(all2) Use all2 from Maekem 100yearsBCI to select the quadrats
plot10ha.1 = plot10ha.temp[plot10ha.temp$gx>250 & plot10ha.temp$gx<750,] #young 10ha area 1 #5ha
plot10ha.2 = plot10ha.temp[plot10ha.temp$gx>800 & plot10ha.temp$gy>50,] #young 10ha area 2 #1ha
plot10ha.old = plot10ha.temp[plot10ha.temp$gx<200,] #old 10haplot ~4ha

plot10ha = rbind(plot10ha.1, plot10ha.2)
plot10ha$site = "P10"; plot10ha$sp_code = plot10ha$sp; plot10ha$dbh = plot10ha$dbhmm.1.3; plot10ha$site_age = 140; #site ages in plot25ha are 144-100 years old, mostly 100
plot10ha = plot10ha[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")]

plot10ha.old$site = "P10"; plot10ha.old$sp_code = plot10ha.old$sp; plot10ha.old$dbh = plot10ha.old$dbhmm.1.3; plot10ha.old$site_age = 140; #site ages in plot25ha are 144-100 years old, mostly 100
plot10ha.old = plot10ha.old[, c("site", "site_age", "sp_code", "dbh", "quadrat_plot")] #2ha

#bind the 5 datasets: Dent, Lang-Knight- Plot25ha, PlotAVA, Plot10ha
lang.k$quadrat_plot = paste(lang.k$q20, "LAN", sep = "_") #q20 is the 20x20m quadrat id in the Lang-Knight treefile
lang.k = lang.k[,c("site", "site_age", "sp6", "dbh1968", "quadrat_plot")]
colnames(lang.k) = c("site", "site_age", "sp_code", "dbh", "quadrat_plot")
baseline = rbind(baseline, lang.k, plot25ha, plotava, plot10ha)
baseline$sp_code = tolower(baseline$sp_code)

#retrieve the species to use based on all of our analysis for the paper
load("ML_mortality_50CIRdata")
data = ML_mortality
species.to.use = tolower(unique(data$sp))

## young forest and mixed vicinity forest are total = 2.36 ha (with 1.92ha of young forest and 0.44ha of mixed forest)

# 50-ha plot
#load("data/updated_50haplotdata_2019/bci.tree/bci.tree8.rdata") #census2015
##where is the young forest areas in the 50ha plot?
#loads the topographic types for each quadrat
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

#Add to the baseline: mature secondary forest
baseline = rbind(baseline, plot50ha.young)

#old growth forest 

#where is the mixed forest in the vicinity of the above mixed forest?
#remove the quadrats that have young forest
#plot50ha.old = plot50ha.old[-which(plot50ha.old$quad %in% mixed.forest.50ha),]; nrow(plot50ha.old) #169457
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

#bind the with 50ha datase datasets of old-growth forest
mature = rbind(plot50ha.old, plotpearson, plotzetek, plotdrayton, plot10ha.old, plotGigante)

#select only the species to use here:
baseline_giants = baseline[baseline$sp_code %in% species.to.use,]
mature_giants = mature[mature$sp_code %in% species.to.use,]

# ##graphing for ONLY giant trees (500 mm)
young = baseline_giants[baseline_giants$dbh >= 500,]
young = young[!young$sp_code== "vataer",] #remove "vataer" as there are not individuals in the old forest, they are mostly in Gigante peninsula

species = (unique(young$sp_code)); length(unique(young$sp_code)) #37 species
old = mature_giants[mature_giants$dbh >= 500,]
old = old[old$sp_code %in% species,]; length(unique(old$sp_code)) #37 species

##graphing for ONLY giant trees species
#create a data.frame that has each species, m_rate, young forest values, old growth forest, so it is easier to graph
data.plots = data.frame("species" = species, "old.diam" = rep(0, length(species)),
                        "young.diam" = rep(0, length(species)))

new_col_names = paste0("dev", 1:1000) #deviation from the 1:1 line
data.plots[new_col_names] = NA

#areas in ha #check this
y.ha = 4+1.5+27.04+4+6+1.92 #(dent+lang+25ha+ava+10ha_young+50ha_young) = 44.46ha
o.ha = 105.88 #4+47.64+6+6+6+36.24 #(10ha_old+50ha_old+pearson+zetek+drayton+gigante_peninsula) = 105.88ha

#pre-split row indices by quadrat_plot once, so each bootstrap draw is a fast lookup
young.quads = unique(baseline$quadrat_plot)
young.quad.idx = split(seq_len(nrow(baseline)), baseline$quadrat_plot)
old.quads = unique(mature$quadrat_plot)
old.quad.idx = split(seq_len(nrow(mature)), mature$quadrat_plot)

funk =  function(x){ #observed 90% quantile diameter
  print(x)
  y.use = young[young$sp_code == x,]
  o.use = old[old$sp_code == x,]
  #use the 90th percentile of giant tree diameter
  data.plots[data.plots$species == x, "young.diam"] <<- quantile(y.use$dbh, probs = 0.90)
  data.plots[data.plots$species == x, "old.diam"] <<- quantile(o.use$dbh, probs = 0.90)
}

a = unlist(lapply(species, funk))

### here: select hectares of the old-growth as the total number of hectares in the mature secondary
## the data is sorted by quadrats 20x20m so hectares in mature secondary = 44.46 ha or 1,111.5 quadrats or 444,600 m2
funk_deviation =  function(x){ #deviation of the old-growth from the 1:1 line (90% quantile diameter)
  print(i)
  y.use = young[young$sp_code == x,]
  o.use = old[old$sp_code == x,]
  #use the 90th percentile of giant tree diameter deviation
  data.plots[data.plots$species == x, paste0("dev",i)] <<- quantile(o.use$dbh, probs = 0.90) - quantile(y.use$dbh, probs = 0.90)
}

for(i in 1:1000){
  quadrats.draw = sample(old.quads, 1111, replace = F) #draws the same number of quadrats from the old-growth
  old = mature[mature$quadrat_plot %in% quadrats.draw,] #selects those quadrats
  old = old[old$dbh >= 500,]
  old = old[old$sp_code %in% species,]; length(unique(old$sp_code))
  tmp = unlist(lapply(species, funk_deviation))
}

## Add taxon to the data
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "taxon"]
  return(tx)
}
data.plots$taxon = unlist(lapply(data.plots$species, funk))
data.plots$taxon2 = sub("^(\\w)\\w*\\s(\\w+)$", "\\1. \\2", data.plots$taxon)

#what percentage of the time old-growth is larger? basically how many times is the deviation positive/total values 
data.plots$positive = rowSums(data.plots[,4:1003]>0, na.rm = T)
data.plots$percentage = (data.plots$positive/1000)*100

### *Figure S6: plot violin plots for deviation
violin = ggplot(data.frame(taxon2 = rep(data.plots$taxon2, times = 1000),
                           dev = as.vector(as.matrix(data.plots[, paste0("dev", 1:1000)])),
                           mort = rep(data$ML_mortality[match(data.plots$species,
                                                              tolower(data$sp))], times = 1000)),
                aes(x = dev, y = reorder(taxon2, mort))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +   #the 1:1 line
  geom_violin(fill = "#ffc425", colour = NA, scale = "width", na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", size = 1.2, na.rm = TRUE) +
  geom_text(data = data.frame(taxon2 = data.plots$taxon2,
                              mort = data$ML_mortality[match(data.plots$species,
                                                             tolower(data$sp))],
                              percentage = data.plots$percentage),
            aes(x = Inf, label = sprintf("%.1f%%", percentage)),
            hjust = 1, size = 2.8, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_basis +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  labs(x = "Deviation in maximum diameter (mm)", y = NULL)
