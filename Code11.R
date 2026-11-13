########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Sapling recruiment and giant tree mortality
## Sapling recruitment is calculates per capita (50ha plot)
## Produces: Figures 4d, S11
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggplot2); library(lmodel2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=17))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 14, face = "bold")) #, face ="bold"

#census-interval data
load("survival_data_assembled") #from CleanData.R
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp <- simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                        simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                        !is.na(simple_combined$dbhmm.1.3),]

#load the species.to.use data
load("ML_mortality_50CI.Rdata")
species.to.use = ML_mortality$sp

#calcualte the abundance (# trees) of the species.to.use in each of the census intervals of
#the 50ha plot only
sc_50ha = simple_combined_sp[simple_combined_sp$plot == "plot50ha" & simple_combined_sp$sp %in% species.to.use,]

#the 50ha plot was censused in 1982,1985,1990,1995,2000,2005,2010,2015; openyear records the
#exact measurement date within a census, so floor() assigns each row to its census year
sc_50ha$census_year = floor(sc_50ha$openyear)

#total abundance (# trees/ha, since the plot is 50 ha) per species per census
species_N_by_census = aggregate(dbhmm.1.3 ~ sp + census_year, data = sc_50ha, FUN = length)
names(species_N_by_census)[names(species_N_by_census) == "dbhmm.1.3"] = "n_trees"
species_N_by_census$n_trees = species_N_by_census$n_trees/50 #trees/ha per census

#one row per species, one column per census year
species_N_wide = reshape(species_N_by_census, idvar = "sp", timevar = "census_year",
                         direction = "wide")
names(species_N_wide) = gsub("n_trees\\.", "N_", names(species_N_wide))

species_N_wide = species_N_wide[,-9] #removes 2015 where only diptpa were visited

#add data from Ruger 2020, Science
ruger = read.table("ruger_data_s1.txt", h =T)

#add each species' sapling recruitment rate (F_2, recruits ha^-1 yr^-1) from Ruger et al. 2020
species_N_wide$F_2 = ruger$F_2[match(species_N_wide$sp, ruger$sp2)]
#add top layer mortality rates
species_N_wide$mort.top = ruger$mu1_2[match(species_N_wide$sp, ruger$sp2)]

#recruitment per unit abundance, for each census: F_2 (recruits ha^-1 yr^-1) / N_<year> (which is in trees ha^-1)
census_cols = grep("^N_", names(species_N_wide), value = TRUE)
for (col in census_cols) {
  species_N_wide[[paste0("rate_", col)]] <- species_N_wide$F_2 / species_N_wide[[col]]
}

#mean recruitment per abundance per species, averaged across all censuses
rate_cols = paste0("rate_", census_cols)
species_N_wide$mean_rate_N = rowMeans(species_N_wide[,rate_cols], na.rm = TRUE)

#now add the ML_mortality rates and make graph
species_N_wide$ML_mortality = ML_mortality$ML_mortality[match(species_N_wide$sp, ML_mortality$sp)]

### *Figure 4d: giant tree mortality and per capita sapling recruitment 
model = lmodel2(log(species_N_wide$ML_mortality)~log(species_N_wide$mean_rate_N), nperm = 99)
plot = ggplot(data = species_N_wide, aes(y = ML_mortality,x = mean_rate_N))+
  geom_point(size = 2.5) + theme_basis +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") +
  scale_x_log10(expression(bold("Sapling recruitment"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01),0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10","20")) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10" )

### *Figure S11: mortality as a tree in the top layer and per capita sapling recruitment 
model2 = lmodel2(log(species_N_wide$mort.top)~log(species_N_wide$mean_rate_N), nperm = 99)
plot2 = ggplot(data = species_N_wide, aes(y = mort.top,x = mean_rate_N))+ geom_point(size = 3) +
  theme_basis +
  scale_y_log10(expression(bold("Top layer mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.04,.01)),
                labels = c("0.1",rep("",8),"1.0","2.0", "3.0", "4.0")) +
  theme(legend.position = "none") +
  scale_x_log10(expression(bold("Sapling recruitment"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01),0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10","20")) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10" )
