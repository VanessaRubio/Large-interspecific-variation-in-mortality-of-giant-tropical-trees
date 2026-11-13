########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Estimates mortality rates of giant trees, and compares them against predictors
## Produces: Figures S3, 1a, 4a, s9, s10, 4b, 4c
########################################################################################################################

#clear workspace
rm(list=ls()); library(dplyr); library(ggplot2); library(lmodel2)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", face ="bold", size = 16))

#census-interval data
load("survival_data_assembled.Rdata") 
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp <- simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                        simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                        !is.na(simple_combined$dbhmm.1.3),]

#select only trees >50cm 
simple_combined_sp_ov50 <- simple_combined_sp[simple_combined_sp$dbhmm.1.3>=500,]

#recalculate opencensus and endcensus based for trees >50cm
simple_combined_sp_ov50 <- simple_combined_sp_ov50 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus_ov50 = opencensus-min(opencensus),
         endcensus_ov50 = opencensus_ov50 + (endcensus-opencensus))

#this is for data that are not time dependent, which includes these data
pertree_ov50 <- simple_combined_sp_ov50 %>%
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1)

#select species with more than 30 individuals
names_ov50_30indivs <- names(table(pertree_ov50$sp)[table(pertree_ov50$sp)>=30]) 
pertree_ov50_30indivs <- pertree_ov50[pertree_ov50$sp %in% names_ov50_30indivs,]

#convert table to data frame
df = as.data.frame(sort(table(pertree_ov50$sp), decreasing = T))
names(df) = c("sp", "count")

#sort
df$sp = factor(df$sp, levels = df$sp)
df = df[-which(df$sp == "dead"), ]

### *Figure S3: Abundance curve for giant tree (DBH 50cm) sample size per species*
figure_s3 = ggplot(df, aes(y = count, x = sp)) + theme_basis+
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(y = "Abundance", x = "Species")+ lims(y=c(0,400))+
  theme(axis.text.x = element_text(size = 5, angle  = 90, hjust = 1))

#### calculate average survival using iterative method, Maximum likelihood method (Kubo et al. 2000) 
#create a for loop to do this with all species with >30 individuals

#remove species that are not represented by at least 30 individuals
allintervals_ov50_30indivs <- simple_combined_sp_ov50[simple_combined_sp_ov50$sp %in% names_ov50_30indivs,]

#calculate census interval
allintervals_ov50_30indivs$interval <- allintervals_ov50_30indivs$endcensus-allintervals_ov50_30indivs$opencensus

#create for-loop
for(j in 1:length(names_ov50_30indivs)){ 
  sp_df <- allintervals_ov50_30indivs[allintervals_ov50_30indivs$sp==names_ov50_30indivs[j],]
  
  #create important metrics for contextualizing the data
  sp_name <- names_ov50_30indivs[j]
  total_years_monitored <- sum(sp_df$interval)
  trees <- length(unique(sp_df$plot_treeID))
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
    print(k)
  }

  ML_rate_loCI <- quantile(resampled_mortality,0.025)
  ML_rate_hiCI <- quantile(resampled_mortality,0.975)

    if(j ==1){
    ML_mortality<-data.frame("sp"=sp_name,"ML_mortality" = numerical_solution$mortality_rates,
                             "ML_rate_loCI" = NA, "ML_rate_hiCI" = NA,
                             "total_years_monitored" = total_years_monitored,
                             "sample_size" = trees, "trees_killed" = trees_killed)
  }else{
    ML_mortality_inter<-data.frame("sp"=sp_name,"ML_mortality" = numerical_solution$mortality_rates, 
                                   "ML_rate_loCI" = NA, "ML_rate_hiCI" = NA,
                                   "total_years_monitored" = total_years_monitored,
                                   "sample_size" = trees, "trees_killed" = trees_killed)
    ML_mortality <- bind_rows(ML_mortality,ML_mortality_inter)
  }
}

######################################################### GRAPHS ############################################################### 
#import species table to add mnemonics for matching species with other data
load("bci.spptable.rdata")

#add species mnemonics
bci.spptable_trim = bci.spptable[,c("sp","Family","Genus","Species")]
bci.spptable_trim$taxon = paste(bci.spptable_trim$Genus,bci.spptable_trim$Species,sep=" ")

#keep only one row for repeats
bci.spptable_trim = bci.spptable_trim %>% group_by(taxon) %>% slice(1)

#add species data to the residuals
data = left_join(ML_mortality,bci.spptable_trim, by = "sp")

#add species codes
data$sp_abrev = paste(substr(data$Genus,1,1),data$Species,sep = ". ")

### *Figure 1a: graph species-specific giant tree mortality rates*
figure_1a = ggplot(data = data, aes(x = ML_mortality, y = reorder(taxon, ML_mortality))) + 
  geom_point()+ theme(legend.position = "none") + theme_basis + labs(y = "Species") +
  theme(axis.text.y = element_text(face = "italic", angle = 0, size = 10, hjust = 1)) +
  scale_x_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.0008, 0.12),
                breaks = c(0.0008,seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.12),
                labels = c("","0.1",rep("",8),"1.0",rep("",8),"10", ""))+ 
  geom_errorbar(aes(xmin = ML_rate_loCI, xmax = ML_rate_hiCI), cex = 0.3)

### add data from Ruger 2020, Science
ruger = read.table("data/updated_50haplotdata_2019/bci.tree/ruger_data_s1.txt", h =T)

#adding PC2 from Ruger's paper
funk= function(species){ 
  use = match(toString(species), ruger$sp2)
  pcs = ruger[use,c(7,8,10,20:28)]
  return(pcs)
}
all.pcs = do.call("rbind", lapply(data$sp, funk)) #applies the function in all the data, unlist+as.numeric to create a vector instead of the list
data = cbind(data,all.pcs)

#selecting the species for which we have complete trait data
data = data[-which(data$sp == "vataer"),]

#### Giant tree mortality vs. PCA1 and PCA2 axes ####

## graphs without Cecropia insignis and Zanthoxylum ekmanii (see Figure 4a)
data2 = data[-which(data$sp == "cecrin"),]
data2 = data2[-which(data2$sp == "zantbe"),]

#type 2 regression
model.m1 = lmodel2(log(data$ML_mortality)~data$PC1score, nperm = 99)
model.m2 = lmodel2(log(data2$ML_mortality)~data2$PC1score, nperm = 99) #without Cecropia insignis and Zanthoxylum ekmanii 

### *Figure 4a: giant tree mortality and PC1 axis*
figure_4a =
  ggplot(data, aes(y = ML_mortality, x = PC1score)) +
  geom_point(aes(color = "blue"))+
  geom_point(data = subset(data, PC1score < -4), color = "#C80CFF", size = 1.1) +
  theme_basis + labs(x = "Growth-survival tradeoff")+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC1score, y = ML_mortality, ymax = ML_rate_hiCI, ymin = ML_rate_loCI, color= "blue"))+
  geom_errorbar(data = subset(data, PC1score < -4),
                aes(ymax = ML_rate_hiCI, ymin = ML_rate_loCI), color = "#C80CFF", width = 0) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  geom_smooth(data = data2, method = "lm", se = FALSE, size = 0.5, linetype = "dashed", col = "grey50")+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = -4.7, y = 0.00135 ,label = paste("R^2==", round(model.m1$rsquare, 
                                                                              digits = 2)), parse = T)+
  annotate(geom = "text", x = -4.7, y = 0.001 ,label = paste("R^2==", round(model.m2$rsquare, 
                                                                            digits = 2)), parse = T)
### *Figure S9: mortality rate vs PC1 scores - residuals*
m1 = lmodel2(log(data$ML_mortality)~data$PC1score, nperm = 99)
intercept.m1 = m1$regression.results$Intercept[2]  
slope.m1     = m1$regression.results$Slope[2]      
fitted.m1    = intercept.m1 + slope.m1 * data$PC1score
residuals.m1 = log(data$ML_mortality) - fitted.m1   
residuals.m1.abs = abs(residuals.m1)

data$fitted.m1    = fitted.m1
data$residuals.m1 = residuals.m1

figure_s9 = ggplot(data, aes(x = PC1score, y = residuals.m1.abs)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey50") +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE, col = "grey10", linewidth = 0.5) +
  theme_basis + labs(x = "PC1 score", y = "abs (residuals)")

### *Figure 4b: giant tree mortality and PC2 axis*
model.m = lmodel2(log(data$ML_mortality)~data$PC2score, nperm = 99)
model.m

figure_4b =
  ggplot(data, aes(y = ML_mortality, x = PC2score)) +
  geom_point(aes(color = "blue"))+
  theme_basis + labs(x = "Stature-recruitment tradeoff")+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") +
  geom_errorbar(aes(x = PC2score, y = ML_mortality, ymax = ML_rate_hiCI, ymin = ML_rate_loCI, color= "blue"))+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey10") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  annotate(geom = "text", x = 0.15, y = 0.001 ,label = paste("R^2==", round(model.m$rsquare, 
                                                                            digits = 2)), parse = T)

### *Figure 4c: giant tree mortality and demographic space*
figure_4c = ggplot(data = ruger, aes(x = PC1score,y =PC2score)) + 
  geom_point(color="grey80") + theme_basis +geom_point(data = data, 
                                                       aes(x = PC1score,y =PC2score,
                                                           color = ML_mortality*100,size = ML_mortality), 
                                                       stroke = 1, alpha = 0.8) +
  guides(alpha = "none", size = "none", color = guide_colorbar(
    direction = "horizontal",
    title.position = "top",     
    barwidth = 5,               
    barheight = 1))+
  labs(x = "Growth-survival tradeoff", y = "Stature-recruitment tradeoff")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"('%'~yr^{-1})))+
  geom_vline(xintercept = 0)+ geom_hline(yintercept = 0)+
  scale_y_continuous(limits =c(-5,5),expand=c(0,0),breaks = seq(-6,6,2),
                     labels = c("", -4, -2, 0, 2, 4, "") ) +
  scale_x_continuous(limits =c(-7.2,6),expand=c(0,0),breaks = seq(-6,6,2),
                     labels = c("", -4, -2, 0, 2, 4, "") )+
  theme(
    legend.background = element_blank(),
    legend.position = c(0.87, 0.89),    
    legend.key = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

# range of PCA captured in ruger vs. this study
sum(abs(range(data$PC1score)))/sum(abs(range(ruger$PC1score))) 
sum(abs(range(data$PC2score)))/sum(abs(range(ruger$PC2score))) 

#% of species above 0 in PC2 that we capture
nrow(data[data$PC2score>=0,])/nrow(ruger[ruger$PC2score>=0,]) 

### *Figure S10: Giant tree mortality vs. growth and mortality in each canopy layer (from Ruger et al. 2020)*
#growth in the top canopy layer vs. giant tree mortality
m1 = lmodel2(log(ML_mortality)~ log(G1_2), data = data, nperm = 99)
p1 = ggplot(data = data, aes(x = G1_2, y = ML_mortality))+
  geom_point(size = 3) + theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") + scale_x_log10(expression(bold("Top layer growth rate"~('mm'~yr^{-1})))) 

#mortality in the top layer vs. giant tree mortality
m2 = lmodel2(log(ML_mortality)~ log(mu1_2), data = data, nperm = 99)
p2 = ggplot(data = data, aes(x = mu1_2, y = ML_mortality))+
  geom_point(size = 3) + theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") + 
  scale_x_log10(expression(bold("Top layer mortality rate"~('%'~yr^{-1}))),
                limits = c(0.001,0.05),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.05,.01)),
                labels = c("0.1",rep("",8),"1.0","", "3.0", "", "5.0")) +
  theme(legend.position = "none") + geom_abline(intercept = 0, slope = 1, col = "gray10")+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed", col = "black" ) 

#growth second layer vs. giant mortality
m3 = lmodel2(log(ML_mortality)~ log(G2_2), data = data, nperm = 99)
p3 = ggplot(data = data, aes(x = G2_2, y = ML_mortality))+
  geom_point(size = 3) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  theme(legend.position = "none") + 
  scale_x_log10(expression(bold("Second layer growth rate"~('mm'~yr^{-1}))))

# mortality second layer vs. giant mortality
m4 = lmodel2(log(ML_mortality)~ log(mu2_2), data = data, nperm = 99)
p4 = ggplot(data = data, aes(x = mu2_2, y = ML_mortality))+
  geom_point(size = 3) + theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  scale_x_log10(expression(bold("Second layer mortality rate"~('%'~yr^{-1}))),
                limits = c(0.001,0.05),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.05,.01)),
                labels = c("0.1",rep("",8),"1.0","", "3.0", "", "5.0")) +
  geom_abline(intercept = 0, slope = 1, col = "gray10")+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed", col = "black" ) 

#growth third layer vs. giant mortality
m5 = lmodel2(log(ML_mortality) ~ log(G3_2), data = data, nperm = 99)
p5 = ggplot(data = data, aes(x = G3_2, y = ML_mortality))+
  geom_point(size = 3) + theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  scale_x_log10(expression(bold("Third layer growth rate"~('mm'~yr^{-1}))))

#mortality third layer vs. giant mortality
m6 = lmodel2(log(ML_mortality) ~ log(mu3_2), data = data, nperm = 99)
p6 = ggplot(data = data, aes(x = mu3_2, y = ML_mortality))+
  geom_point(size = 3) + theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  scale_x_log10(expression(bold("Third layer mortality rate"~('%'~yr^{-1}))),
                labels = function(x) x * 100) +
  geom_abline(intercept = 0, slope = 1, col = "gray10")+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed", col = "black" ) 

#growth four layer vs. giant mortality
m7 = lmodel2(log(ML_mortality) ~ log(G4_2), data = data, nperm = 99)
p7 = ggplot(data = data, aes(x = G4_2, y = ML_mortality))+
  geom_point(size = 3) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  scale_x_log10(expression(bold("Understory layer growth rate"~('mm'~yr^{-1}))))

#mortality four layer vs. giant mortality
m8 = lmodel2(log(ML_mortality) ~ log(mu4_2), data = data)
p8 = ggplot(data = data, aes(x = mu4_2, y = ML_mortality))+
  geom_point(size = 3) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) +
  scale_x_log10(expression(bold("Understory layer mortality rate"~('%'~yr^{-1}))),
                limits = c(0.001,0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0","", "3.0", "", "5.0", "", "7.0", "", "9.0", "")) +
  geom_abline(intercept = 0, slope = 1, col = "gray10")+
  geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed", col = "black" ) 
