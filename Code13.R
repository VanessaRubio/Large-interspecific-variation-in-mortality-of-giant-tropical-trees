########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Estimation of mean mortality rates by size class (1-10cm; 10-50cm; >50cm for all trees and species combined)
########################################################################################################################

#clear workspace
rm(list=ls()); library(dplyr); library(ggplot2)

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#Drop "*", "", and "unk?"
simple_combined_sp <- simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                        simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                        !is.na(simple_combined$dbhmm.1.3),]

###Calculate mortality using a single record for each tree category
#select only trees >50cm and other class sizes
simple_combined_sp_ov50 <- simple_combined_sp[simple_combined_sp$dbhmm.1.3>=500,]

simple_combined_sp_ov1_10 <- simple_combined_sp[simple_combined_sp$dbhmm.1.3>0 &
                                                  simple_combined_sp$dbhmm.1.3<100,]

simple_combined_sp_ov10_50 <- simple_combined_sp[simple_combined_sp$dbhmm.1.3>=100 &
                                                   simple_combined_sp$dbhmm.1.3<500,]

#recalculate opencensus and endcensus based for trees >50cm
simple_combined_sp50 <- simple_combined_sp_ov50 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus_ov50 = opencensus-min(opencensus),
         endcensus_ov50 = opencensus_ov50 + (endcensus-opencensus))

simple_combined_sp1_10 <- simple_combined_sp_ov1_10 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus_ov1_10 = opencensus-min(opencensus),
         endcensus_ov1_10 = opencensus_ov1_10 + (endcensus-opencensus))

simple_combined_sp10_50 <- simple_combined_sp_ov10_50 %>% group_by(plot_treeID) %>% 
  arrange(opencensus, .by_group = T) %>%
  mutate(opencensus_ov10_50 = opencensus-min(opencensus),
         endcensus_ov10_50 = opencensus_ov10_50 + (endcensus-opencensus))

#calculate census interval
simple_combined_sp_ov50$interval <- simple_combined_sp_ov50$endcensus-simple_combined_sp_ov50$opencensus
simple_combined_sp_ov1_10$interval <- simple_combined_sp_ov1_10$endcensus-simple_combined_sp_ov1_10$opencensus
simple_combined_sp_ov10_50$interval <- simple_combined_sp_ov10_50$endcensus-simple_combined_sp_ov10_50$opencensus

datasets = list(
  d.1.10 = simple_combined_sp_ov1_10, 
  d.10.50 = simple_combined_sp_ov10_50, 
  d.50 = simple_combined_sp_ov50
)

#estimates of size-class mortality rates of all species and individuals together
for(i in 1:length(datasets)) {
  class.size = datasets[[i]]
  total_years_monitored <- sum(class.size$interval)
  trees <- length(unique(class.size$treeID))
  trees_killed <- sum(class.size$mortality)
  
  #create vector with intervals for all dead and all surviving trees
  dead_intervals <- class.size$interval[class.size$mortality == 1]
  living_intervals <- class.size$interval[class.size$mortality == 0]
  
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
  resampled_mortality <- rep(NA,100000)
  resampled_minimizations<-rep(NA,100001)
  potential_rates <- seq(0,1,1/100000)
  
  for(k in 1:1000){
    #create vector with intervals for all dead and all surviving trees
    resampled_sp_df <- sample_n(as.data.frame(class.size),size=nrow(class.size),replace=TRUE)
    resampled_dead_intervals <- resampled_sp_df$interval[resampled_sp_df$mortality == 1]
    resampled_living_intervals <- resampled_sp_df$interval[resampled_sp_df$mortality == 0]
    
    for(l in 1:100001) {
      resampled_minimizations[l] <- abs(-sum(resampled_living_intervals) + sum((resampled_dead_intervals*exp(-potential_rates[l]*resampled_dead_intervals))/(1-exp(-potential_rates[l]*resampled_dead_intervals))))
    }
    resampled_mortality[k] <- potential_rates[which.min(resampled_minimizations)]
  }
  ML_rate_loCI <- quantile(na.omit(resampled_mortality),0.025)
  ML_rate_hiCI <- quantile(na.omit(resampled_mortality),0.975)
  
  ML_mortality1 = data.frame("sp"= size,"ML_mortality" = numerical_solution$mortality_rates,
                             "ML_rate_loCI" = ML_rate_loCI, "ML_rate_hiCI" = ML_rate_hiCI,
                             "total_years_monitored" = total_years_monitored,
                             "sample_size" = trees, "trees_killed" = trees_killed, "class" = size)
  
  ML_mortality = data.frame("sp"= size,"ML_mortality" = numerical_solution$mortality_rates,
                            "ML_rate_loCI" = ML_rate_loCI, "ML_rate_hiCI" = ML_rate_hiCI,
                            "total_years_monitored" = total_years_monitored,
                            "sample_size" = trees, "trees_killed" = trees_killed, "class" = size)
  #run after 2nd run
  ML_mortality <- bind_rows(ML_mortality1, ML_mortality)
  name = paste0("ML_mortality_size_class_", names(datasets)[i], ".Rdata")
  save(ML_mortality, file = name)
}

