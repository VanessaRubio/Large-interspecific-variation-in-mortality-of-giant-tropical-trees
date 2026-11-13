########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Estimates mortality rates of small size class trees of the same species as the giant trees
## Produces: Figures 1c, 1d, 1e, S8
########################################################################################################################

#clear workspace
rm(list=ls()); library(dplyr)

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp <- simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                        simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                        !is.na(simple_combined$dbhmm.1.3),]

### select only the species in this study (39 species)
load("ML_mortality_50CI.Rdata")
species.to.use = unique(ML_mortality$sp)

simple_combined_sp <- simple_combined_sp[simple_combined_sp$sp %in% species.to.use,]

#select trees with DBH between 1-10cm
simple_combined_1_10 = simple_combined_sp[simple_combined_sp$dbhmm.1.3>=10 & simple_combined_sp$dbhmm.1.3<100,]
#select trees with DBH between 10-50cm
simple_combined_10_50 = simple_combined_sp[simple_combined_sp$dbhmm.1.3>=100 & simple_combined_sp$dbhmm.1.3<500,]
#select trees with DBH between 1-50cm
simple_combined_1_50 = simple_combined_sp[simple_combined_sp$dbhmm.1.3>=10 & simple_combined_sp$dbhmm.1.3<500,]

#recalculate opencensus and endcensus based for trees >50cm
simple_combined_1_10 <- simple_combined_1_10 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus = opencensus-min(opencensus),
         endcensus = opencensus + (endcensus-opencensus))

simple_combined_10_50 <- simple_combined_10_50 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus = opencensus-min(opencensus),
         endcensus = opencensus + (endcensus-opencensus))

simple_combined_1_50 <- simple_combined_1_50 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus = opencensus-min(opencensus),
         endcensus = opencensus + (endcensus-opencensus))

#this is for data that are not time dependent, which includes these data
pertree_1_10 <- simple_combined_1_10 %>%
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)

pertree_10_50 <- simple_combined_10_50 %>%
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)

pertree_1_50 <- simple_combined_1_50 %>%
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)

#select species with more than 14-30 individuals, this way the rates can be calculated for a bigger number of species
names_1_10_15indivs <- names(table(pertree_1_10$sp)[table(pertree_1_10$sp)>=14])
pertree_1_10_15indivs <- pertree_1_10[pertree_1_10$sp %in% names_1_10_15indivs,]

names_10_50_15indivs <- names(table(pertree_10_50$sp)[table(pertree_10_50$sp)>=14]) 
pertree_10_50_15indivs <- pertree_10_50[pertree_10_50$sp %in% names_10_50_15indivs,]

names_1_50_30indivs <- names(table(pertree_1_50$sp)[table(pertree_1_50$sp)>=30])
pertree_1_50_30indivs <- pertree_1_50[pertree_1_50$sp %in% names_1_50_30indivs,]

#calculate average survival using iterative method, Maximum likelihood method (Kubo et al. 2000) -- ###
##################################################################################
# repeat the following code for all the size classes
allintervals_1_10_15indivs <- simple_combined_1_10[simple_combined_1_10$sp %in% names_1_10_15indivs,]

#calculate census interval
allintervals_1_10_15indivs$interval <- allintervals_1_10_15indivs$endcensus-allintervals_1_10_15indivs$opencensus

#create for-loop
for(j in 1:length(names_1_10_15indivs)){
  sp_df <- allintervals_1_10_15indivs[allintervals_1_10_15indivs$sp==names_1_10_15indivs[j],]
  
  #create important metrics for contextualizing the data
  sp_name <- names_1_10_15indivs[j]
  total_years_monitored <- sum(sp_df$interval)
  trees <- length(unique(sp_df$treeID))
  trees_killed <- sum(sp_df$mortality)
  
  #create vector with intervals for all dead and all surviving trees
  dead_intervals <- sp_df$interval[sp_df$mortality == 1]
  living_intervals <- sp_df$interval[sp_df$mortality == 0]
  
  #perform manual minimization
  minimizations<-rep(NA,100001)
  potential_rates <- seq(0,1,1/100000)
  for(i in 1:100001) {
    minimizations[i] <- abs(-sum(living_intervals) + sum((dead_intervals*exp(-potential_rates[i]*dead_intervals))/(1-exp(-potential_rates[i]*dead_intervals))))
  }
  numerical_approach <- data.frame("minimization_attempts" = minimizations,
                                   "mortality_rates" = potential_rates)
  #select the minimum of the absolute value to get the closest value to zero
  numerical_solution <- numerical_approach[which.min(numerical_approach$minimization_attempts),]
  numerical_solution
  
  #calculate confidence intervals for mortality here
  resampled_mortality <- rep(NA,1000)
  resampled_minimizations<-rep(NA,100001)
  potential_rates <- seq(0,1,1/100000)
  
  for(k in 1:1000){
    #create vector with intervals for all dead and all surviving trees
    resampled_sp_df <- sample_n(as.data.frame(sp_df),size=nrow(sp_df),replace=TRUE)
    resampled_dead_intervals <- resampled_sp_df$interval[resampled_sp_df$mortality == 1]
    resampled_living_intervals <- resampled_sp_df$interval[resampled_sp_df$mortality == 0]
    
    for(l in 1:100001) {
      resampled_minimizations[l] <- abs(-sum(resampled_living_intervals) + sum((resampled_dead_intervals*exp(-potential_rates[l]*resampled_dead_intervals))/(1-exp(-potential_rates[l]*resampled_dead_intervals))))
    }
    resampled_mortality[k] <- potential_rates[which.min(resampled_minimizations)]
  }
  ML_rate_loCI <- quantile(resampled_mortality,0.025)
  ML_rate_hiCI <- quantile(resampled_mortality,0.975)
  
  #combine results into single data frame for the maximum likelihood mortality approach (Kubo et al. 2000)
  if(j == 1){
    ML_mortality<-data.frame("sp"=sp_name,"ML_mortality" = numerical_solution$mortality_rates,
                             "ML_rate_loCI" = ML_rate_loCI, "ML_rate_hiCI" = ML_rate_hiCI,
                             "total_years_monitored" = total_years_monitored,
                             "sample_size" = trees, "trees_killed" = trees_killed)
    
  }else{
    ML_mortality_inter<-data.frame("sp"=sp_name,"ML_mortality" = numerical_solution$mortality_rates, 
                                   "ML_rate_loCI" = ML_rate_loCI, "ML_rate_hiCI" = ML_rate_hiCI,
                                   "total_years_monitored" = total_years_monitored,
                                   "sample_size" = trees, "trees_killed" = trees_killed)
    ML_mortality <- bind_rows(ML_mortality,ML_mortality_inter)
  }
}

######################################################### GRAPHS ############################################################### 
library(lmodel2); library(ggplot2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12, face ="bold")) #, face ="bold"

#giant tree mortality rates
load("ML_mortality_50CI.Rdata"); m.giant = ML_mortality
#1-10cm dbh tree mortality rates
load("ML_mortality_1_10.Rdata"); m.1.10 = ML_mortality
#10-50cm dbh tree mortality rates
load("ML_mortality_10_50.Rdata"); m.10.50 = ML_mortality
#1-50cm dbh tree mortality rates
load("ML_mortality_1_50.Rdata"); m.1.50 = ML_mortality

#bind datasets together
funk= function(species, dataset){ 
  use = match(toString(species), dataset$sp)
  ml = dataset[use,c("ML_mortality", "ML_rate_loCI", "ML_rate_hiCI")]
  return(ml)
}

allml = do.call("rbind", lapply(m.giant$sp, funk, dataset = m.1.50)) 
colnames(allml) = c("ML1.50", "ML1.50.loCI", "ML1.50.hiCI")
m.giant = cbind(m.giant,allml)

allml = do.call("rbind", lapply(m.giant$sp, funk, dataset = m.1.10)) 
colnames(allml) = c("ML1.10", "ML1.10.loCI", "ML1.10.hiCI")
m.giant = cbind(m.giant,allml)

allml = do.call("rbind", lapply(m.giant$sp, funk, dataset = m.10.50)) 
colnames(allml) = c("ML10.50", "ML10.50.loCI", "ML10.50.hiCI")
m.giant = cbind(m.giant,allml)

#here I add taxon for graphing purposes
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "taxon"]
  return(tx)
}
m.giant$taxon = unlist(lapply(m.giant$sp, funk))

#adding the data from Rüger 2020, Science
ruger = read.table("ruger_data_s1.txt", h =T)

#adding species-specific PCAs values to the data
funk= function(species){ 
  use = match(toString(species), ruger$sp2)
  pcs = ruger[use,c("PC1score", "PC2score")]
  return(pcs)
}
all.pcs = do.call("rbind", lapply(m.giant$sp, funk))
m.giant = cbind(m.giant, all.pcs)
m.giant$taxon2 = sub("^(\\w)\\w*\\s(\\w+)$", "\\1. \\2", m.giant$taxon)

#here I only run linear regressions for species with m rates 2%%
species.to.use2 = m.giant[m.giant$ML_mortality<0.0200,]$sp

#organize the data so the datasets share the same species names as different size classes can have different species
m.1.10b = m.1.10[m.1.10$sp %in% species.to.use2,]
temp.giant = m.giant[m.giant$sp %in% m.1.10b$sp,]
m.1.10b$ML50 = temp.giant$ML_mortality; temp.giant$sp%in%m.1.10b$sp
m.1.10b$ML50_loCI = temp.giant$ML_rate_loCI
m.1.10b$ML50_hiCI = temp.giant$ML_rate_hiCI

m.10.50b = m.10.50[m.10.50$sp %in% species.to.use2,]
temp.giant = m.giant[m.giant$sp %in% m.10.50b$sp,]
m.10.50b$ML50 = temp.giant$ML_mortality
m.10.50b$ML50_loCI = temp.giant$ML_rate_loCI
m.10.50b$ML50_hiCI = temp.giant$ML_rate_hiCI

m.1.50b = m.1.50[m.1.50$sp %in% species.to.use2,]
temp.giant = m.giant[m.giant$sp %in% m.1.50b$sp,]
m.1.50b$ML50 = temp.giant$ML_mortality
m.1.50b$ML50_loCI = temp.giant$ML_rate_loCI
m.1.50b$ML50_hiCI = temp.giant$ML_rate_hiCI

### *Figure 1c: giant tree mortality vs. 1-50cm tree mortality
model.1.50 = lmodel2(log(m.1.50$ML50)~ log(m.1.50$ML_mortality), nperm = 99)
model.1.50
model.1.50b = lmodel2(log(m.1.50b$ML50)~ log(m.1.50b$ML_mortality), nperm = 99)
model.1.50b
plot.1.50 = 
  ggplot(data = m.1.50, aes(y = ML50, x = ML_mortality)) + 
  geom_point(color = "grey50") +
  theme_basis + geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "black")+ 
  geom_point(data = m.1.50[m.1.50$ML50 <= 0.02,], shape = 21, fill = NA, color = "red", 
             size = 2, stroke = 1.5)+
  scale_x_log10(expression(bold("1-50cm tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.20),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "solid", linewidth = 0.2)+
  geom_smooth(data = m.1.50b, method = "lm", se = FALSE, linetype = "22", col = "red")

### *Figure 1d: giant tree mortality vs. 1-10cm tree mortality
model.1.10 = lmodel2(log(m.1.10$ML50)~ log(m.1.10$ML_mortality), nperm = 99)
model.1.10
model.1.10b = lmodel2(log(m.1.10b$ML50)~ log(m.1.10b$ML_mortality), nperm = 99)
model.1.10b
plot.1.10 = 
  ggplot(data = m.1.10, aes(y = ML50, x = ML_mortality)) + 
  geom_point(color = "grey50") +
  theme_basis + geom_smooth(method = "lm", se = FALSE, linetype = "22", color = "black")+ 
  geom_point(data = m.1.10[m.1.10$ML50 <= 0.02,], shape = 21, fill = NA, color = "red", 
             size = 2, stroke = 1.5)+
  scale_x_log10(expression(bold("1-10cm tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.20),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "solid", size = 0.2)+
  geom_smooth(data = m.1.10b, method = "lm", se = FALSE, linetype = "22", col = "red")

### *Figure 1e: giant tree mortality vs. 10-50cm tree mortality
model.10.50 = lmodel2(log(m.10.50$ML50)~ log(m.10.50$ML_mortality), nperm = 99)
model.10.50
model.10.50b = lmodel2(log(m.10.50b$ML50)~ log(m.10.50b$ML_mortality), nperm = 99)
model.10.50b
plot.10.50 = 
  ggplot(data = m.10.50, aes(y = ML50, x = ML_mortality)) + 
  geom_point(color = "grey50") +
  geom_point(data = m.10.50[m.10.50$ML50 <= 0.02,], shape = 21, fill = NA, color = "red", 
             size = 2, stroke = 1.5)+
  theme_basis + geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "black")+
  scale_x_log10(expression(bold("10-50cm tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.20),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.18),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20")) +
  geom_abline(slope = 1, intercept = 0, color = "gray10", linetype = "solid", size = 0.2)+
  geom_smooth(data = m.10.50b, method = "lm", se = FALSE, linetype = "22", col = "red")

##graphs of the 1-10 and 10-50 small size class mortality and the growth-survival and stature-recruitment tradeoff
data = m.giant

### *Figure S8a: mortality 1-10cm vs PC1 (growth-survival tradeoff)
model.pcm2 = lm(log(data$ML1.10)~data$PC1score)
species_PC1_m2 = 
  ggplot(data, aes(y = ML1.10, x = PC1score)) +
  geom_point(aes(color = "blue"))+
  theme_basis + labs(x = "Growth-survival tradeoff")+
  scale_y_log10(expression(bold("1-10cm tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"),
                limits = c(0.001, 0.2)) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC1score, y = ML1.10, ymax = ML1.10.hiCI, ymin = ML1.10.loCI, color= "blue"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = -4.3, y = 0.001 ,label = paste("R^2==", round(summary(model.pcm2)$adj.r.squared, 
                                                                            digits = 3)), parse = T)+
  annotate(geom = "text", x = -4.5, y = 0.0016, label = ifelse(summary(model.pcm2)$coefficients[2, "Pr(>|t|)"]<0.01,"P < 0.01",
                                                               paste("P = ", round(summary(model.pcm2)$coefficients[2, "Pr(>|t|)"], digits = 3))))
### *Figure S8c: mortality 10-50cm vs PC1 (growth-survival tradeoff)
model.pcm1 = lm(log(data$ML10.50)~data$PC1score)
species_PC1_m1 =
  ggplot(data, aes(y = ML10.50, x = PC1score)) +
  geom_point(aes(color = "blue"))+
  theme_basis + labs(x = "Growth-survival tradeoff")+
  scale_y_log10(expression(bold("10-50cm tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"),
                limits = c(0.001, 0.12)) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC1score, y = ML10.50, ymax = ML10.50.hiCI, ymin = ML10.50.loCI, color= "blue"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = -4.7, y = 0.0012 ,label = paste("R^2==", round(summary(model.pcm1)$adj.r.squared, 
                                                                             digits = 3)), parse = T)+
  annotate(geom = "text", x = -4.7, y = 0.0016, label = ifelse(summary(model.pcm1)$coefficients[2, "Pr(>|t|)"]<0.01,"P < 0.01",
                                                               paste("P = ", round(summary(model.pcm1)$coefficients[2, "Pr(>|t|)"], digits = 3))))
### *Figure S8b: mortality 1-10cm vs PC2 (stature-recruitment tradeoff)
model.pc2m2 = lm(log(data$ML1.10)~data$PC2score)
species_PC2_m2 =
  ggplot(data, aes(y = ML1.10, x = PC2score)) +
  geom_point(aes(color = "blue"))+
  theme_basis + labs(x = "Stature-recruitment tradeoff")+
  scale_y_log10(expression(bold("1-10cm tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"),
                limits = c(0.001, 0.2)) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC2score, y = ML1.10, ymax = ML1.10.hiCI, ymin = ML1.10.loCI, color= "blue"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = 0.3, y = 0.0010 ,label = paste("R^2==", round(summary(model.pc2m2)$adj.r.squared, 
                                                                            digits = 3)), parse = T)+
  annotate(geom = "text", x = 0.3, y = 0.0016, label = ifelse(summary(model.pc2m2)$coefficients[2, "Pr(>|t|)"]<0.01,"P < 0.01",
                                                              paste("P = ", round(summary(model.pc2m2)$coefficients[2, "Pr(>|t|)"], digits = 3))))
### *Figure S8d: mortality 10-50cm vs PC2 (stature-recruitment tradeoff)
model.pc2m1 = lm(log(data$ML10.50)~data$PC2score)
species_PC2_m1 =
  ggplot(data, aes(y = ML10.50, x = PC2score)) +
  geom_point(aes(color = "blue"))+
  theme_basis + labs(x = "Stature-recruitment tradeoff")+
  scale_y_log10(expression(bold("10-50cm tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"),
                limits = c(0.001, 0.12)) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC2score, y = ML10.50, ymax = ML10.50.hiCI, ymin = ML10.50.loCI, color= "blue"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = 0.3, y = 0.0012 ,label = paste("R^2==", round(summary(model.pc2m1)$adj.r.squared, 
                                                                            digits = 3)), parse = T)+
  annotate(geom = "text", x = 0.3, y = 0.0016, label = ifelse(summary(model.pc2m1)$coefficients[2, "Pr(>|t|)"]<0.01,"P < 0.01",
                                                              paste("P = ", round(summary(model.pc2m1)$coefficients[2, "Pr(>|t|)"], digits = 3))))
