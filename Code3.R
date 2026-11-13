########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Calculates AGB and diameter of giant trees
## Produces: Figures S2a, 1b
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggplot2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12, face ="bold"))

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                       simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                       !is.na(simple_combined$dbhmm.1.3),]

#total forest carbon stocks in a census period
simple_combined_sp_subset = simple_combined_sp[simple_combined_sp$openyear == 2005 & simple_combined_sp$endyear == 2010 |
                                                 simple_combined_sp$openyear == 2003 & simple_combined_sp$endyear == 2008 |
                                                 simple_combined_sp$openyear == 2004 & simple_combined_sp$endyear == 2014 |
                                                 simple_combined_sp$openyear == 2000 & simple_combined_sp$endyear == 2014,]

simple_combined_sp_subset = simple_combined_sp_subset[simple_combined_sp_subset$mortality == 0,] #alive trees

#remove the trees from the quadrats that only have trees >=20cm in the gigante peninsula, "gigante_plot" = -1
simple_combined_sp_subset = simple_combined_sp_subset[simple_combined_sp_subset$gigante_plot >= 0,] #trees that are in plots that measure trees >=1cm

## Add taxon to the data
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "taxon"]
  return(tx)
}
simple_combined_sp_subset$taxon = unlist(lapply(simple_combined_sp_subset$sp, funk))

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  f = bci.spptable[use, "Family"]
  return(f)
}
simple_combined_sp_subset$family = unlist(lapply(simple_combined_sp_subset$sp, funk))

#adding wood density to the dataframe for later biomass calculations
jw_data = read.table("WD_species.txt", h = T)
jw_data = data.frame(jw_data)
jw_data$SP. = tolower(jw_data$sp6)

#adding wood density in g/cm3
funk= function(species){ 
  use = match(toString(species), jw_data$SP.)
  wd = jw_data[use, "wd100.mean"]
  return(wd)
}
simple_combined_sp_subset$wd = unlist(lapply(simple_combined_sp_subset$sp, funk))
simple_combined_sp_subset = as.data.frame(simple_combined_sp_subset)

#extras wood from ctfs data when the species did not have reported wood density in the previous dataset
load("wsg.ctfs.Rdata")
ctfs = wsg.ctfs3
ctfs$taxon = paste(ctfs$genus, ctfs$species, sep=" ")
with.nas = simple_combined_sp_subset[is.na(simple_combined_sp_subset$wd),]
unique.with.nas = unique(with.nas$taxon)
ctfs[which(ctfs$taxon == "Acacia melanoceras"), "taxon"] =  "Vachellia melanoceras"

funk= function(tax){ 
  use = match(toString(tax), ctfs$taxon)
  wd2 = ctfs[use,"wsg"]
  simple_combined_sp_subset[which(simple_combined_sp_subset$taxon == tax), "wd"] <<- wd2
}
q = unlist(lapply(unique.with.nas, funk))

#for species with No wood density value: the mean of all other species in the site was used
with.nas = simple_combined_sp_subset[is.na(simple_combined_sp_subset$wd),]; nrow(with.nas)
(nrow(with.nas)/nrow(simple_combined_sp_subset))*100 #0.08% of individuals dont have wd
unique.with.nas = unique(with.nas$sp)

simple_combined_sp_subset[which(simple_combined_sp_subset$sp %in% unique.with.nas), "wd"] =
  mean(na.omit(simple_combined_sp_subset$wd))

#calculate their AGB
#when height is unknown - Chave et al., 2015 (Equation 7) for moist tropical forest
#D in "cm"; wsg in "g/cm3"
#Kg/1000 to convert to Mg
#dbhmm.1.3 is already corrected with taper equations when the data was cleaned

#For palm species -> ln(ABG) = -3.3488 + 2.7484*ln(dbh) * CF (in Kilograms)
#CF correction factor by log-transforming the data see Chave et al., 2005 (CF = exp((RSE^2)/2) = 1.188714)
palms = subset(simple_combined_sp_subset, family == "Arecaceae")
#palms = simple_combined_sp_subset[simple_combined_sp_subset$family == "Arecaceae",]
palms$AGB = (exp(-3.3488 + 2.7484*log(palms$dbhmm.1.3/10))/1000) * (1.188714)

trees = subset(simple_combined_sp_subset, family != "Arecaceae")
#trees = simple_combined_sp_subset[!simple_combined_sp_subset$family == "Arecaceae",]
E = 0.051763982 #for BCI
trees$AGB = (exp(-1.803-0.976*E+(0.976*log(trees$wd))+
                   (2.673*log(trees$dbhmm.1.3/10))-
                   (0.0299*(log(trees$dbhmm.1.3/10))^2))/1000)

simple_combined_sp_subset = rbind(trees, palms)

#total hectares: 50ha; gigante 36.24 ha; 25ha 27.04ha; 10ha 10ha; AVA 4ha; pearson 6ha; zetek 6ha; drayton 6ha; knight-lang 1.5ha
#we use all the data because we only need individuals higher than 50cm
hectares.used = 146.34

#relationship between diameter and giant tree mortality
funk = function(x, tmp){
  use = tmp[tmp$sp == x, ]
  mean_dbh = mean(use$dbhmm.1.3, na.rm = TRUE)
  return(mean_dbh)
}

#mean DBH (mm)
data$mean.DBH = sapply(species.to.use, funk, tmp = simple_combined_sp_subset2.giants)
#CIs
dbhs = list()
for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  idx = sample(nrow(simple_combined_sp_subset2.giants), 
               nrow(simple_combined_sp_subset2.giants), 
               replace = TRUE)
  tmp = simple_combined_sp_subset2.giants[idx,]
  dbhs[[i]] = unlist(lapply(species.to.use, funk, tmp))
}

dbh_mat = do.call(rbind, dbhs)  # rows = bootstrap reps, cols = species
colnames(dbh_mat) = species.to.use

CIs = apply(dbh_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_df = as.data.frame(t(CIs))
ci_df$sp = species.to.use 
colnames(ci_df) = c("diam.loCI", "diam.hiCI", "sp")

data = merge(data, ci_df, by.x = "sp", by.y = "sp", all.x = TRUE)

#above ground biomass (Mg)
funk = function(x, tmp){
  use = tmp[tmp$sp == x,]
  mean_agb = mean(use$AGB, na.rm = TRUE)
  return(mean_agb)
}
#mean AGB 
data$mean.AGB = sapply(species.to.use, funk, tmp = simple_combined_sp_subset2.giants)
#CIs
agbs = list()
for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  idx = sample(nrow(simple_combined_sp_subset2.giants), 
               nrow(simple_combined_sp_subset2.giants), 
               replace = TRUE)
  tmp = simple_combined_sp_subset2.giants[idx,]
  agbs[[i]] = unlist(lapply(species.to.use, funk, tmp))
}

agb_mat = do.call(rbind, agbs)  # rows = bootstrap reps, cols = species
colnames(agb_mat) = species.to.use

CIs = apply(agb_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_df = as.data.frame(t(CIs))
ci_df$sp = species.to.use 
colnames(ci_df) = c("agb.loCI", "agb.hiCI", "sp")

data = merge(data, ci_df, by.x = "sp", by.y = "sp", all.x = TRUE)
 
### *Figure S2a: plot AGB and giant tree mortality
model = (lm(log(data$mean.AGB) ~ log(data$ML_mortality)))
summary(model)
cross = ggplot(data = data, aes(x = ML_mortality, y = mean.AGB)) + theme_basis +
  geom_point(col = "lightcoral", cex = 1.5) +
  geom_errorbar(aes(ymin = agb.loCI, ymax = agb.hiCI), cex = 0.3, col = "lightcoral")+
  geom_errorbar(aes(xmin = ML_rate_loCI, xmax = ML_rate_hiCI), cex = 0.3, col = "lightcoral")+
  scale_x_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20"))+
  scale_y_log10(expression(bold("Mean AGB of giant trees (Mg)")),
                breaks = c(seq(0,9,1), seq(10,60,10)),
                labels = c("0",rep("",9),"10","20","30","40","50","60"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10")

### *Figure 1b: plot DBH and giant tree mortality
model2dbh = (lm(log(data$mean.DBH) ~ log(data$ML_mortality)))
summary(model2dbh)
cross2 = ggplot(data = data, aes(x = ML_mortality, y = mean.DBH)) + theme_basis +
  geom_point(col = "lightcoral", cex = 1.5) +
  geom_errorbar(aes(ymin = diam.loCI, ymax = diam.hiCI), cex = 0.3, col = "lightcoral")+
  geom_errorbar(aes(xmin = ML_rate_loCI, xmax = ML_rate_hiCI), cex = 0.3, col = "lightcoral")+
  scale_x_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20"))+
  scale_y_log10(expression(bold("Mean DBH of giant trees (mm)")),
                limits = c(500, 2000))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") 
