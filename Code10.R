########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Giant tree mortality and functional traits/charcteristics related to size-dependent drivers of mortality 
## size-dependent drivers of mortality: drought, wind, and lightning
## Produces: Figures S12, 5i, 5j, 5k
########################################################################################################################

#clear workspace
rm(list=ls()); library(openxlsx); library(ggplot2); library(readxl)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", face ="bold", size = 12))

load("ML_mortality_50CI.Rdata")
species.to.use = ML_mortality$sp

##### Here loads the Hydraulics ##### 
hyd = read.csv(file = "ks species summary.csv")
output_df_plotting_comb = ML_mortality

#match the data - species and P50 (Stem water potential at which 50% of maximum stem hydraulic conductivity is lost)
funk = function(x){
  use = match(x, tolower(hyd$sp))
  h = hyd[use, "p50"]
  return(h)
}
output_df_plotting_comb$p50 = unlist(lapply(output_df_plotting_comb$sp, funk))
nrow(output_df_plotting_comb[!is.na(output_df_plotting_comb$p50),]) #only 12 species

species.hydraulics = output_df_plotting_comb[!is.na(output_df_plotting_comb$p50),]$sp

### *Figure S12a: giant tree mortality vs. P50 (Stem water potential at which 50% of maximum stem hydraulic conductivity is lost) 
model = lm(data = output_df_plotting_comb,log(ML_mortality) ~ p50)
h = ggplot(data = output_df_plotting_comb, aes(x = p50, y = ML_mortality)) + 
  labs(x = "P50 (MPa)")+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"))+
  geom_errorbar(aes(ymin = ML_rate_loCI, ymax = ML_rate_hiCI), cex = 0.3)+
  geom_point() + theme_basis 

#load turgor loss point data
turgor = read_excel("Kunert-et-al_2021_NewPhyt_SItable_HM.xlsx", sheet = "modified")
turgor = as.data.frame(turgor)

#match the data - species and TLP (turgor loss point)
funk = function(x){
  use = match(x, tolower(turgor$spec6))
  t = turgor[use, "TLPmean"]
  return(t)
}
output_df_plotting_comb$TLPmean = unlist(lapply(output_df_plotting_comb$sp, funk))
nrow(output_df_plotting_comb[!is.na(output_df_plotting_comb$TLPmean),]) #only 25 species
shapiro.test(output_df_plotting_comb$TLPmean) #0.2352 do not transform

### *Figure 5j: giant tree mortality vs. TLP 
model2 = lm(data = output_df_plotting_comb, log(ML_mortality) ~ TLPmean)
t = ggplot(data = output_df_plotting_comb2, aes(x = TLPmean, y = ML_mortality)) + 
  labs(x = "Turgor loss point (MPa)")+
  theme(plot.background = element_rect(fill = "grey98")) +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"))+
  geom_point(size = 3) + theme_basis +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10" )+
  geom_errorbar(aes(ymin = ML_rate_loCI, ymax = ML_rate_hiCI), cex = 0.5)

#run the analysis only with the 12 species that we have data for P50 "species.hydarulics"
subset.species = output_df_plotting_comb[output_df_plotting_comb$sp %in% species.hydarulics,]
nrow(output_df_plotting_comb[!is.na(output_df_plotting_comb$TLPmean),]) #only 25 species

### *Figure S12c: giant tree mortality vs. TLP for smaller sample size, same species that we have the data for P50
model3 = lm(data = subset.species, log(ML_mortality) ~ TLPmean)
t3 = ggplot(data = subset.species, aes(x = TLPmean, y = ML_mortality)) + 
  labs(x = "Turgor loss point (MPa)")+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"))+
  geom_point() + theme_basis +
  geom_errorbar(aes(ymin = ML_rate_loCI, ymax = ML_rate_hiCI), cex = 0.3)

##### Here loads the Deciduousness ##### 
#load deciduousness data
deci = read.table("deciduousness.txt", h = T)
deci$sp6 = tolower(deci$sp6)

#match the data
funk = function(x){
  use = match(x, deci$sp6)
  dec = deci[use, "deciduous"]
  return(dec)
}
output_df_plotting_comb$deciduous = unlist(lapply(output_df_plotting_comb$sp, funk))

## Add taxon to the data
load("bci.spptable.rdata")
bci.spptable$taxon = paste(bci.spptable$Genus, bci.spptable$Species, sep=" ")

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "taxon"]
  output_df_plotting_comb[output_df_plotting_comb$sp == species, "taxon"] <<- tx
}
q = unlist(lapply(output_df_plotting_comb$sp, funk))

#remove the species with unknown deciduous data, and divide the data into a new variable for deciduous (D) vs. evergreen (E)
output_df_plotting_comb = output_df_plotting_comb[-which(output_df_plotting_comb$deciduous == "UNK"),]
output_df_plotting_comb$new.deciduous = ifelse(output_df_plotting_comb$deciduous == "DB" | 
                                                 output_df_plotting_comb$deciduous == "DF" |
                                                 output_df_plotting_comb$deciduous == "DO", "D","E")
output_df_plotting_comb$deciduous = as.factor(output_df_plotting_comb$deciduous)

## plots
### *Figure S12b: giant tree mortality vs. leaf habit 
bx = ggplot(data = output_df_plotting_comb, aes(x = deciduous, y = ML_mortality, 
                                                fill = new.deciduous)) + 
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10"))+
  geom_boxplot(alpha = 0.3) + theme_basis + labs(x = "Leaf habit")+
  scale_fill_manual(values = c("D" = "tan1", "E" = "darkgreen")) +theme(legend.position = "none")

### *Figure 5i: giant tree mortality vs. leaf habit (new.deciduous, i.e. only D vs. E)
model.decid = aov(log(ML_mortality) ~ new.deciduous, data = output_df_plotting_comb)
summary(model.decid)

bx2 = ggplot(data = output_df_plotting_comb, aes(x = new.deciduous, y = ML_mortality, fill = new.deciduous)) + 
  theme(plot.background = element_rect(fill = "grey98")) +
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01),0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10",""))+
  geom_boxplot(alpha = 0.3) + theme_basis + labs(x = "Leaf habit")+ theme(legend.position = "none")+
  scale_fill_manual(values = c("D" = "tan1", "E" = "darkgreen"))

##### Here loads the Wind mortality risk ##### 

##this code uses Toby Jackson's equations from the paper "Wind Shapes the Growth Strategies of Trees in a Tropical Forest" (2024; Ecology letters)
#from Jackson's script: select their trees with dbh higher than 50cm
#wind.risk.37 = are mean values of risk of snapping at 37km/hper species (same order as the subset_paper) from Jackson's data for trees higher than 50cm
#wind.risk.rel = are mean relative values (100*risk37/maxrisk) per species (same order as the subset_paper) from Jackson's data for trees higher than 50cm 
load("dt_analysis_t50.Rdata"); load("ML_mortality_50CI.Rdata")

models <- lapply(unique(dt_analysis_t50$genus), function(g) {
  subset_data <- dt_analysis_t50[dt_analysis_t50$genus == g, ]
  lm(rel_risk ~ tree_dbh_m, data = subset_data)
})

#maximum risk for species with dbh>=50cm
all.50 = dt_analysis_t50[dt_analysis_t50$tree_dbh_m >=0.50,]
max.risk = max(all.50$risk_at_37)

#create a table to store the data
table.wind = data.frame("species" = c("tab1gu", "diptpa", "virosu", "anacex", "jac1co"), 
                        "mean.rel.risk" = rep(0,5), "low.CI" = rep(0,5), "high.CI" = rep(0,5))

Handroanthus = as.data.frame(dt_analysis_t50[dt_analysis_t50$genus == "Handroanthus",])
model.H = summary(lm(rel_risk ~ tree_dbh_m, data = Handroanthus))
H.50 = Handroanthus[Handroanthus$tree_dbh_m >=0.5,]
H.50 %<>% mutate(rel_risk.new = 100*(risk_at_37/max.risk))
table.wind[table.wind$species == "tab1gu", "mean.rel.risk"] = mean(H.50$rel_risk.new)

Dypterix = as.data.frame(dt_analysis_t50[dt_analysis_t50$genus == "Dypterix",])
model.D = summary(lm(rel_risk ~ tree_dbh_m, data = Dypterix))
D.50 = Dypterix[Dypterix$tree_dbh_m >=0.5,]
D.50 %<>% mutate(rel_risk.new = 100*(risk_at_37/max.risk))
table.wind[table.wind$species == "diptpa", "mean.rel.risk"] = mean(D.50$rel_risk.new)

Virola = as.data.frame(dt_analysis_t50[dt_analysis_t50$genus == "Virola",])
model.V = summary(lm(rel_risk ~ tree_dbh_m, data = Virola))
V.50 = Virola[Virola$tree_dbh_m >=0.5,]
V.50 %<>% mutate(rel_risk.new = 100*(risk_at_37/max.risk))
table.wind[table.wind$species == "virosu", "mean.rel.risk"] = mean(V.50$rel_risk.new)

Anacardium = as.data.frame(dt_analysis_t50[dt_analysis_t50$genus == "Anacardium",])
model.A = summary(lm(rel_risk ~ tree_dbh_m, data = Anacardium))
A.50 = Anacardium[Anacardium$tree_dbh_m >=0.5,]
A.50 %<>% mutate(rel_risk.new = 100*(risk_at_37/max.risk))
table.wind[table.wind$species == "anacex", "mean.rel.risk"] = mean(A.50$rel_risk.new)

Jacaranda = as.data.frame(dt_analysis_t50[dt_analysis_t50$genus == "Jacaranda",])
model.J = summary(lm(rel_risk ~ tree_dbh_m, data = Jacaranda))
J.50 = Jacaranda[Jacaranda$tree_dbh_m >=0.5,]
J.50 %<>% mutate(rel_risk.new = 100*(risk_at_37/max.risk))
table.wind[table.wind$species == "jac1co", "mean.rel.risk"] = mean(J.50$rel_risk.new)

datasets <- list(
  tab1gu = H.50, #Handroanthus
  diptpa = D.50, #Dypterix
  virosu = V.50, #Virola
  anacex = A.50, #Anacardium
  jac1co = J.50 #Jacaranda
)
#calculate confidence intervals of the wind mortality risk
for(sp in names(datasets)){
  to.use = datasets[[sp]]
  
  tmp = replicate(1000,{
    resampled = sample_n(to.use, size = nrow(to.use), replace = TRUE)
    mean(resampled$rel_risk.new)
  })
  to.use.loCI = quantile(tmp, 0.025)
  to.use.hiCI = quantile(tmp, 0.975)
  table.wind[table.wind$species == sp, "low.CI"]  = to.use.loCI
  table.wind[table.wind$species == sp, "high.CI"] = to.use.hiCI
}
#organize the data
funk= function(species){ 
  use = match(toString(species), ML_mortality$sp)
  ml = ML_mortality[use, c("ML_mortality","ML_rate_loCI","ML_rate_hiCI")]
  return(ml)
}
all.t = do.call("rbind", lapply(table.wind$species, funk))
table.wind = cbind(table.wind,all.t)
model.wind = lm(log(ML_mortality)~ mean.rel.risk, data = table.wind)

### *Figure 5k: giant tree mortality and wind mortality risk
wind = ggplot(data = table.wind , aes(x = mean.rel.risk, y = ML_mortality))+
  theme(plot.background = element_rect(fill = "grey98")) +
  theme_basis + geom_point(col = "black", data = table.wind, size = 3, aes(x = mean.rel.risk,y = ML_mortality))+
  xlab("Wind mortality risk (relative)") + xlim(-1,60) + 
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.03,.01)),
                labels = c("0.1",rep("",8),"1.0","2.0","3.0"))+
  geom_errorbar(aes(x = mean.rel.risk, y = ML_mortality, xmax = high.CI, xmin = low.CI))+
  geom_errorbar(aes(x = mean.rel.risk, y = ML_mortality, ymax = ML_rate_hiCI, ymin = ML_rate_loCI))

##### Here loads the lightning strikes ##### 
#Lightning model from Gora et al., 2020 JEcol
lightning = read.csv("final_frequency_ov50.csv") 

#species survival-mortality rates 
#select our 39 giant species
lightning = lightning[lightning$Mnemonic %in% species.to.use,]
output = ML_mortality

#add the new columns to the data
measures = c("prob.direct", "prob.indirect", "prob.aggregated",
              "freq.direct", "freq.aggregated")
new.cols = as.vector(t(outer(measures, c("", ".loCI", ".hiCI"), paste0)))
output[new.cols] = NA

#calculate averages and CIs for only 36/39 species
funk= function(x){ 
  use = lightning[lightning$Mnemonic == x,]
  output[output$sp == x,"prob.direct"] <<- mean(na.omit(use$prob.direct.strike))
  output[output$sp == x,"prob.direct.loCI"] <<- quantile(na.omit(use$prob.direct.strike),0.025)
  output[output$sp == x,"prob.direct.hiCI"] <<- quantile(na.omit(use$prob.direct.strike),0.975)
  output[output$sp == x,"prob.indirect"] <<- mean(na.omit(use$indirect_dmg_probability))
  output[output$sp == x,"prob.indirect.loCI"] <<- quantile(na.omit(use$indirect_dmg_probability),0.025)
  output[output$sp == x,"prob.indirect.hiCI"] <<- quantile(na.omit(use$indirect_dmg_probability),0.975)
  output[output$sp == x,"prob.aggregated"] <<- mean(na.omit(use$agg_dmg_probability))
  output[output$sp == x,"prob.aggregated.loCI"] <<- quantile(na.omit(use$agg_dmg_probability),0.025)
  output[output$sp == x,"prob.aggregated.hiCI"] <<- quantile(na.omit(use$agg_dmg_probability),0.975)
  output[output$sp == x,"freq.direct"] <<- mean(na.omit(use$freq.direct.strike))
  output[output$sp == x,"freq.direct.loCI"] <<- quantile(na.omit(use$freq.direct.strike),0.025)
  output[output$sp == x,"freq.direct.hiCI"] <<- quantile(na.omit(use$freq.direct.strike),0.975)
  output[output$sp == x,"freq.aggregated"] <<- mean(na.omit(use$agg_dmg_frequency))
  output[output$sp == x,"freq.aggregated.loCI"] <<- quantile(na.omit(use$agg_dmg_frequency),0.025)
  output[output$sp == x,"freq.aggregated.hiCI"] <<- quantile(na.omit(use$agg_dmg_frequency),0.975)
}
a = unlist(lapply(output$sp, funk))

output$species.names = output$sp
output$strikes = output$freq.direct/output$ML_mortality #mean freq direct strikes in their lifetime (1/m)

#removing Alchornea costaricensis (no direct probability), and Enterolobium sc and vataer (no data) 
output = output[output$sp != "alchco" & output$sp != "entesc"  & output$sp !="vataer",]

### *Figure S12d: giant tree mortality and direct lightning strikes
mod =  lm(log(output$ML_mortality)~(log(output$strikes)))
dir.strikes = ggplot(data = output, aes(y = ML_mortality, x = strikes)) +
  geom_point(data = output, size = 3, aes(x = strikes,y = ML_mortality), 
             stroke = 0.5, alpha = 0.8)+ geom_vline(xintercept = 1,linetype="dotted")+
  theme_basis + 
  scale_x_log10(expression(bold("Direct lightining strike")))+ 
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10","20")) +
  theme(legend.position = "none")
