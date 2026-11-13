########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Recruitment into the 50cm size class in old-growth and mature secondary forest 
## Produces: Figure S7
########################################################################################################################

#clear workspace
rm(list=ls()); library(openxlsx); library(ggplot2); library(dplyr)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", size = 12))

#census-interval data
load("survival_data_assembled.rec.Rdata") 
simple_combined = as.data.frame(simple_combined.rec2)

#Drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*"&simple_combined$sp!=""&
                                       simple_combined$sp!="unk?"&!is.na(simple_combined$sp)&
                                       !is.na(simple_combined$dbhmm.1.3),]

#plots with secondary forest = (dent+lang+25ha+ava+10ha_young+50ha_young)
plots_secondary = c("plot25ha", "plotKnight", "plotAVA", "plot10ha", "plot50ha")
load("young.forest.Rdata") #baseline
baseline = baseline

simple_combined_sp_young = simple_combined_sp[which(simple_combined_sp$plot %in% plots_secondary),]
secondary = simple_combined_sp_young[which(simple_combined_sp_young$plot_treeID %in% baseline$plot_treeID),]

###Calculate recruitment using a single record for each tree over 50cm
#recalculate opencensus and endcensus based for all trees once they became >=500mm
#this is when they start being giant trees
simple_rec_census = secondary %>%
  mutate(cat_dbh = dbhmm.1.3 >= 500) %>% #is their current dbh >=500? YES/true (was gigante already, "initial count"), NO/false
  group_by(plot_treeID, cat_dbh) %>% 
  arrange(opencensus, .by_group = T) %>% 
  mutate(opencensus_ov50 = opencensus-min(opencensus), #open census as a giant (only uses the rows with cat_dbh = TRUE)
         endcensus_ov50 = opencensus_ov50 + (endcensus-opencensus))%>% #years as giant (true above)
  mutate(opencensus_ov50 = case_when(dbhmm.1.3 < 500 ~ NA, #if it is not giant, do not change these columns
                                     T ~ opencensus_ov50),
         endcensus_ov50 = case_when(dbhmm.1.3 < 500 ~ NA,
                                    T ~ endcensus_ov50)) %>%
  data.frame()

#calculate open and close census for new individuals, since they just became 500mm then the time as giant is just
#the difference between the last census and the open census
#this includes also trees that are not giant, basically all trees that are not already a giant or just became 
simple_rec_census[simple_rec_census$cat_dbh == FALSE, "opencensus_ov50"] = 0
simple_rec_census[simple_rec_census$cat_dbh == FALSE, "endcensus_ov50"] = 
  simple_rec_census[simple_rec_census$cat_dbh == FALSE, "endcensus"]-
  simple_rec_census[simple_rec_census$cat_dbh == FALSE, "opencensus"]

#select the last census for each tree.
temporal = simple_rec_census %>% 
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1) #select the last row (last census where it was recorded)
temporal = as.data.frame(temporal)
sum(is.na(temporal$next.dbhmm.1.3)); length(which(temporal$mortality ==1)) #~same number of deaths as NAs in next.dbh

#which trees enter in the analysis: 
#we use the trees with nextdbh higher than 500 in the last census, new.recruits = new.dbh>500 and cat.dbh = False (but here use all becuase they were not giants already; all because we need false and true)
#the trees that dies (NA as nex.dbhmm.1.3) but had a dbh>500
to.use = which(temporal$next.dbhmm.1.3>=500|(is.na(temporal$next.dbhmm.1.3)&temporal$dbhmm.1.3>=500))
simple_combined_sp_ov50 = temporal[to.use,]

#select only the species in this analysis
#retrieve the species to use based on all of our analysis
load("ML_mortality_50CI.Rdata")
data = ML_mortality
species.to.use = tolower(unique(data$sp))

#trim baseline to use only our species
sp_indivs = simple_combined_sp_ov50[simple_combined_sp_ov50$sp %in% species.to.use,];nrow(sp_indivs) # 931 trees

#a couple trees have a NA as mortality, that means they are NEW recruits into the plot and NEW recruits into the 50cm
#change endcensus_ov50 to "0" (are just new-new recruits); and mortality to "0" (alive)
sp_indivs[is.na(sp_indivs$mortality), "mortality"] = 0
sp_indivs[is.na(sp_indivs$mortality), "endcensus_ov50"] = 0

#species to use
names_ov50_15indivs = unique(sp_indivs$sp) #17 species

###Recruitment rates per plot, per census interval, and per species

#restrict to the same giant-tree individuals used above (species.to.use, and
#trees that reach >=500mm), but keep EVERY cat_dbh==TRUE census row per tree (not
#just the last one, as "temporal"/sp_indivs does) so we can break rates out by census interval
giant_treeIDs = unique(simple_combined_sp_ov50$plot_treeID)

percensus_df = simple_rec_census[simple_rec_census$plot_treeID %in% giant_treeIDs &
                                   simple_rec_census$sp %in% species.to.use,]

#for cat_dbh==FALSE rows, keep ONLY the interval where the tree actually crosses
#into >=500mm (current dbh <500, next.dbhmm.1.3 >=500). Otherwise a tree would be
#counted as a "new recruit" (adding to NT) in every earlier sub-500mm census of its
#history, not just the one interval where it actually transitions into the giant class
is_recruit_interval = percensus_df$cat_dbh == FALSE &
  !is.na(percensus_df$next.dbhmm.1.3) & percensus_df$next.dbhmm.1.3 >= 500
percensus_df = percensus_df[percensus_df$cat_dbh == TRUE | is_recruit_interval,]

#a tree could in principle satisfy that crossing condition more than once if dbh
#measurements oscillate around 500mm; keep only the earliest crossing per tree
percensus_df = percensus_df %>%
  group_by(plot_treeID, cat_dbh) %>%
  arrange(opencensus, .by_group = TRUE) %>%
  filter(cat_dbh == TRUE | row_number() == 1) %>%
  ungroup() %>%
  data.frame()

#same NA-mortality handling as above: NA mortality = new-new recruit (alive, just entered)
percensus_df[is.na(percensus_df$mortality), "mortality"] = 0
percensus_df[is.na(percensus_df$mortality), "endcensus_ov50"] = 0

output_df_percensus_young = data.frame()

#for loop to calculate recruitment rates using Kohyama et al.2017, per plot and per census interval
for(p in seq_along(plots_secondary)){
  
  plot_df = percensus_df[percensus_df$plot == plots_secondary[p],]
  censuses = sort(unique(plot_df$openyear)) #open year is more informative
  
  for(c in seq_along(censuses)){
    
    census_df = plot_df[plot_df$openyear == censuses[c],]
    sp_list = unique(census_df$sp)
    
    for(i in seq_along(sp_list)){
      print(paste(plots_secondary[p], "openyear:",censuses[c], sp_list[i]))
      inter_df <- census_df[census_df$sp == sp_list[i],]
      
      N = nrow(inter_df[inter_df$cat_dbh == TRUE,]) ##initial count of individuals in time i (were more than 50cm at the census)
      S = N - sum(inter_df[inter_df$cat_dbh == TRUE, "mortality"]) #count of survivors (total giant trees (cat_dbh = T)- giant trees that died)
      NT = nrow(inter_df[inter_df$cat_dbh == FALSE,]) + S #final count in time i+1, what was already alive in the previous census + new recruits into the 50cm size class
      
      meant = mean(na.omit(inter_df$endcensus_ov50)) #mean time between censuses
            r_rate = (log(NT)-log(S))/meant
      
      #skip groups too small/degenerate to produce a usable rate (S=0, no time span, etc.)
      if(!is.na(r_rate) & is.finite(r_rate) & r_rate != 0){
        resampled_recruitment = rep(NA,1000)
        for(j in 1:1000){
          resampled_rec = sample_n(as.data.frame(inter_df),size=nrow(inter_df),replace=TRUE)
          S_j = nrow(resampled_rec[resampled_rec$cat_dbh == TRUE,]) - sum(resampled_rec[resampled_rec$cat_dbh == TRUE, "mortality"]) #count of survivors
          NT_j = nrow(resampled_rec[resampled_rec$cat_dbh == FALSE,]) + S_j
          
          resampled_recruitment[j] = (log(NT_j) - log(S_j))/mean(na.omit(resampled_rec$endcensus_ov50))
        }
        
        resampled_recruitment[which(is.infinite(resampled_recruitment))] = mean(resampled_recruitment[which(is.finite(resampled_recruitment))])/2 #if there are Inf in the calculation change it for the mean/2
        
        r_rate_loCI = quantile(resampled_recruitment,0.025, na.rm = TRUE)
        r_rate_hiCI = quantile(resampled_recruitment,0.975, na.rm = TRUE)
        inter_output_df = data.frame("plot" = plots_secondary[p], "census" = censuses[c],
                                     "sp" = sp_list[i], "trees" = N, "N_start" = N, "meant" = meant,
                                     "trees_dead" = sum(inter_df$mortality),
                                     "r_rate" = r_rate,"r_rate_loCI" = r_rate_loCI, "r_rate_hiCI" = r_rate_hiCI)
      }else{
        
        inter_output_df = data.frame("plot" = plots_secondary[p], "census" = censuses[c],
                                     "sp" = sp_list[i], "trees" = N, "N_start" = N, "meant" = meant,
                                     "trees_dead" = sum(inter_df$mortality),
                                     "r_rate" = r_rate,"r_rate_loCI" = r_rate, "r_rate_hiCI" = r_rate)
      }
      output_df_percensus_young <- bind_rows(output_df_percensus_young,inter_output_df)
    }
  }
}

###Weighted-average recruitment rate in secondary per species, pooling across plot

#weight each plot x census-interval rate by (interval length in years) x (N_start:
#individuals with cat_dbh==TRUE at the start of the interval, i.e. already >=500mm;
#this excludes new recruits, which have cat_dbh==FALSE and enter the count only in NT)
output_df_percensus_young$weight = output_df_percensus_young$meant * output_df_percensus_young$N_start

#drop interval-plot cases where r_rate could not be estimated (non-finite, e.g. S=0 or no time span)
#or that carry no weight (no giant trees at the start of the interval and/or zero-length interval)
weighting_df_young = output_df_percensus_young[is.finite(output_df_percensus_young$r_rate) &
                                                 is.finite(output_df_percensus_young$weight) &
                                                 output_df_percensus_young$weight > 0,]

names_percensus = unique(weighting_df_young$sp)
output_df_weighted_young = data.frame()

for(i in seq_along(names_percensus)){
  sp_df = weighting_df_young[weighting_df_young$sp == names_percensus[i],]
  
  r_rate_weighted = weighted.mean(sp_df$r_rate, w = sp_df$weight)
  r_rate_loCI_weighted = weighted.mean(sp_df$r_rate_loCI, w = sp_df$weight)
  r_rate_hiCI_weighted = weighted.mean(sp_df$r_rate_hiCI, w = sp_df$weight)
  
  inter_output_df = data.frame("sp" = names_percensus[i],
                               "n_intervals" = nrow(sp_df),
                               "N_start" = sum(sp_df$N_start),
                               "trees_dead" = sum(sp_df$trees_dead),
                               "r_rate" = r_rate_weighted,
                               "r_rate_loCI" = r_rate_loCI_weighted,
                               "r_rate_hiCI" = r_rate_hiCI_weighted)
  output_df_weighted_young <- bind_rows(output_df_weighted_young,inter_output_df)
}

###plots with old-growth forest

plots_old = c("plot50ha", "plotPearson", "plotZetek", "plotDrayton", "Gigante")
load("old.forest.Rdata") #baseline
old = old

simple_combined_sp_old = simple_combined_sp[which(simple_combined_sp$plot %in% plots_old),]
oldgrowth = simple_combined_sp_old[which(simple_combined_sp_old$plot_treeID %in% old$plot_treeID),]

###Calculate recruitment using a single record for each tree over 50cm 
#recalculate opencensus and endcensus based for all trees once they became >=500mm
#this is when they start being giant trees

simple_rec_census.o = oldgrowth %>%
  mutate(cat_dbh = dbhmm.1.3 >= 500) %>% #is their current dbh >=500? YES/true (was gigante already, "initial count"), NO/false
  group_by(plot_treeID, cat_dbh) %>% 
  arrange(opencensus, .by_group = T) %>% 
  mutate(opencensus_ov50 = opencensus-min(opencensus), #open census as a giant (only uses the rows with cat_dbh = TRUE)
         endcensus_ov50 = opencensus_ov50 + (endcensus-opencensus))%>% #years as giant (true above)
  mutate(opencensus_ov50 = case_when(dbhmm.1.3 < 500 ~ NA, #if it is not giant, do not change these columns
                                     T ~ opencensus_ov50),
         endcensus_ov50 = case_when(dbhmm.1.3 < 500 ~ NA,
                                    T ~ endcensus_ov50)) %>%
  data.frame()

#calculate open and close census for new individuals, since they just became 500mm then the time as giant is just
#the difference between the last census and the open census
#this includes also trees that are not giant, basically all trees that are not already a giant or just became 
simple_rec_census.o[simple_rec_census.o$cat_dbh == FALSE, "opencensus_ov50"] = 0
simple_rec_census.o[simple_rec_census.o$cat_dbh == FALSE, "endcensus_ov50"] = 
  simple_rec_census.o[simple_rec_census.o$cat_dbh == FALSE, "endcensus"]-
  simple_rec_census.o[simple_rec_census.o$cat_dbh == FALSE, "opencensus"]

#select the last census for each tree.
temporal.o = simple_rec_census.o %>% 
  group_by(plot_treeID) %>% arrange(desc(endcensus)) %>% slice(1) #select the last row (last census where it was recorded)
temporal.o = as.data.frame(temporal.o)
sum(is.na(temporal.o$next.dbhmm.1.3)); length(which(temporal.o$mortality ==1)) #~same number of deaths as NAs in next.dbh

#which trees enter in the analysis: 
#we use the trees with nextdbh higher than 500 in the last census, new.recruits = new.dbh>500 and cat.dbh = False (but here use all becuase they were not giants already; all because we need false and true)
#the trees that dies (NA as nex.dbhmm.1.3) but had a dbh>500
to.use.o = which(temporal.o$next.dbhmm.1.3>=500|(is.na(temporal.o$next.dbhmm.1.3)&temporal.o$dbhmm.1.3>=500))
simple_combined_sp_ov50.o = temporal.o[to.use.o,]

length(unique(simple_combined_sp_ov50.o$plot_treeID)); length(unique(simple_combined_sp_ov50.o$sp)); nrow(simple_combined_sp_ov50.o) #1218 unique plot_treeIDs of 99 species

pertree_ov50_15indivs.o <- simple_combined_sp_ov50.o %>% group_by(plot_treeID) %>%
  arrange(opencensus, .by_group = T)
pertree_ov50_15indivs.o = pertree_ov50_15indivs.o%>%arrange(plot_treeID)
length(unique(pertree_ov50_15indivs.o$sp))

#trim baseline to use only our species
pertree_ov50_15indivs.o = pertree_ov50_15indivs.o[pertree_ov50_15indivs.o$sp %in% species.to.use,] # 2150 trees
length(unique(pertree_ov50_15indivs.o$sp)) #39 species

#a couple trees have a NA as mortality, that means they are NEW recruits into the plot and NEW recruits into the 50cm
#change endcensus_ov50 to "0" (are just new-new recruits); and mortality to "0" (alive)
pertree_ov50_15indivs.o[is.na(pertree_ov50_15indivs.o$mortality), "mortality"] = 0
pertree_ov50_15indivs.o[is.na(pertree_ov50_15indivs.o$mortality), "endcensus_ov50"] = 0

#species to use
names_ov50_15indivs.o = unique(pertree_ov50_15indivs.o$sp) 

###Recruitment rates per plot, per census interval, and per species

#restrict to the same giant-tree individuals used above (species.to.use, and
#trees that reach >=500mm), but keep EVERY cat_dbh==TRUE census row per tree (not
#just the last one, as "temporal.o"/pertree_ov50_15indivs.o does) so we can break rates out by census interval
giant_treeIDs.o = unique(simple_combined_sp_ov50.o$plot_treeID)

percensus_df.o = simple_rec_census.o[simple_rec_census.o$plot_treeID %in% giant_treeIDs.o &
                                       simple_rec_census.o$sp %in% species.to.use,]

#for cat_dbh==FALSE rows, keep ONLY the interval where the tree actually crosses
#into >=500mm (current dbh <500, next.dbhmm.1.3 >=500). Otherwise a tree would be
#counted as a "new recruit" (adding to NT) in every earlier sub-500mm census of its
#history, not just the one interval where it actually transitions into the giant class
is_recruit_interval.o = percensus_df.o$cat_dbh == FALSE &
  !is.na(percensus_df.o$next.dbhmm.1.3) & percensus_df.o$next.dbhmm.1.3 >= 500
percensus_df.o = percensus_df.o[percensus_df.o$cat_dbh == TRUE | is_recruit_interval.o,]

#a tree could in principle satisfy that crossing condition more than once if dbh
#measurements oscillate around 500mm; keep only the earliest crossing per tree
percensus_df.o = percensus_df.o %>%
  group_by(plot_treeID, cat_dbh) %>%
  arrange(opencensus, .by_group = TRUE) %>%
  filter(cat_dbh == TRUE | row_number() == 1) %>%
  ungroup() %>%
  data.frame()

#same NA-mortality handling as above: NA mortality = new-new recruit (alive, just entered)
percensus_df.o[is.na(percensus_df.o$mortality), "mortality"] = 0
percensus_df.o[is.na(percensus_df.o$mortality), "endcensus_ov50"] = 0

output_df_percensus.o = data.frame()

#for loop to calculate recruitment rates using Kohyama et al.2017, per plot and per census interval
for(p in seq_along(plots_old)){
  
  plot_df.o = percensus_df.o[percensus_df.o$plot == plots_old[p],]
  censuses.o = sort(unique(plot_df.o$openyear)) #open year is more informative
  
  for(c in seq_along(censuses.o)){
    
    census_df.o = plot_df.o[plot_df.o$openyear == censuses.o[c],]
    sp_list.o = unique(census_df.o$sp)
    
    for(i in seq_along(sp_list.o)){
      inter_df <- census_df.o[census_df.o$sp == sp_list.o[i],]
      N = nrow(inter_df[inter_df$cat_dbh == TRUE,]) ##initial count of individuals in time i (were more than 50cm at the census)
      S = N - sum(inter_df[inter_df$cat_dbh == TRUE, "mortality"]) #count of survivors (total giant trees (cat_dbh = T)- giant trees that died)
      NT = nrow(inter_df[inter_df$cat_dbh == FALSE,]) + S #final count in time i+1, what was already alive in the previous census + new recruits into the 50cm size class
      
      meant = mean(na.omit(inter_df$endcensus_ov50)) #mean time between censuses
      
      #r_rate = (log(final count) - log(count of survivors))/time duration between two censuses (years)
      r_rate = (log(NT)-log(S))/meant
      
      #skip groups too small/degenerate to produce a usable rate (S=0, no time span, etc.)
      if(!is.na(r_rate) & is.finite(r_rate) & r_rate != 0){
        resampled_recruitment = rep(NA,1000)
        for(j in 1:1000){
          # r_rate = (log(final count) - log(count of survivors))/time duration between two censuses (years)
          resampled_rec = sample_n(as.data.frame(inter_df),size=nrow(inter_df),replace=TRUE)
          
          S_j = nrow(resampled_rec[resampled_rec$cat_dbh == TRUE,]) - sum(resampled_rec[resampled_rec$cat_dbh == TRUE, "mortality"]) #count of survivors
          NT_j = nrow(resampled_rec[resampled_rec$cat_dbh == FALSE,]) + S_j
          
          resampled_recruitment[j] = (log(NT_j) - log(S_j))/mean(na.omit(resampled_rec$endcensus_ov50))
        }
        
        resampled_recruitment[which(is.infinite(resampled_recruitment))] = mean(resampled_recruitment[which(is.finite(resampled_recruitment))])/2 #if there are Inf in the calculation change it for the mean/2
        
        r_rate_loCI = quantile(resampled_recruitment,0.025, na.rm = TRUE)
        r_rate_hiCI = quantile(resampled_recruitment,0.975, na.rm = TRUE)
        inter_output_df = data.frame("plot" = plots_old[p], "census" = censuses.o[c],
                                     "sp" = sp_list.o[i], "trees" = N, "N_start" = N, "meant" = meant,
                                     "trees_dead" = sum(inter_df$mortality),
                                     "r_rate" = r_rate,"r_rate_loCI" = r_rate_loCI, "r_rate_hiCI" = r_rate_hiCI)
      }else{
        
        inter_output_df = data.frame("plot" = plots_old[p], "census" = censuses.o[c],
                                     "sp" = sp_list.o[i], "trees" = N, "N_start" = N, "meant" = meant,
                                     "trees_dead" = sum(inter_df$mortality),
                                     "r_rate" = r_rate,"r_rate_loCI" = r_rate, "r_rate_hiCI" = r_rate)
      }
      output_df_percensus.o <- bind_rows(output_df_percensus.o,inter_output_df)
    }
  }
}

head(output_df_percensus.o)

###Weighted-average recruitment rate in old-growth per species, pooling across plot

#weight each plot x census-interval rate by (interval length in years) x (N_start:
#individuals with cat_dbh==TRUE at the start of the interval, i.e. already >=500mm;
#this EXCLUDES new recruits, which have cat_dbh==FALSE and enter the count only in NT)
output_df_percensus.o$weight = output_df_percensus.o$meant * output_df_percensus.o$N_start

#drop interval-plot cases where r_rate could not be estimated (non-finite, e.g. S=0 or no time span)
#or that carry no weight (no giant trees at the start of the interval and/or zero-length interval)
weighting_df.o = output_df_percensus.o[is.finite(output_df_percensus.o$r_rate) &
                                         is.finite(output_df_percensus.o$weight) &
                                         output_df_percensus.o$weight > 0,]

names_percensus.o = unique(weighting_df.o$sp)
output_df_weighted.o = data.frame()

for(i in seq_along(names_percensus.o)){
  sp_df = weighting_df.o[weighting_df.o$sp == names_percensus.o[i],]
  
  r_rate_weighted = weighted.mean(sp_df$r_rate, w = sp_df$weight)
  r_rate_loCI_weighted = weighted.mean(sp_df$r_rate_loCI, w = sp_df$weight)
  r_rate_hiCI_weighted = weighted.mean(sp_df$r_rate_hiCI, w = sp_df$weight)
  
  inter_output_df = data.frame("sp" = names_percensus.o[i],
                               "n_intervals" = nrow(sp_df),
                               "N_start" = sum(sp_df$N_start),
                               "trees_dead" = sum(sp_df$trees_dead),
                               "r_rate" = r_rate_weighted,
                               "r_rate_loCI" = r_rate_loCI_weighted,
                               "r_rate_hiCI" = r_rate_hiCI_weighted)
  output_df_weighted.o <- bind_rows(output_df_weighted.o,inter_output_df)
}

##### graphs ###
##add mortality rates to the recruitment rates for graphing purposes
funk = function(x){ 
  use = match(toString(x), data$sp)
  m = data[use, c("ML_mortality", "ML_rate_loCI", "ML_rate_hiCI")]
  return(m)
}
all = do.call("rbind", lapply(output_df_weighted_young$sp, funk))
output_df_weighted_young = cbind(output_df_weighted_young, all)
all = do.call("rbind", lapply(output_df_weighted.o$sp, funk))
output_df_weighted.o = cbind(output_df_weighted.o, all)

#add taxon
load("bci.spptable.rdata")

funk= function(species){ 
  use = match(toString(species), bci.spptable$sp)
  tx = bci.spptable[use, "Latin"]
  return(tx)
}
output_df_weighted_young$taxon = unlist(lapply(output_df_weighted_young$sp, funk))
output_df_weighted_young$taxon2 = sub("^(\\w)\\w*\\s(\\w+)$", "\\1. \\2", output_df_weighted_young$taxon)
output_df_weighted.o$taxon = unlist(lapply(output_df_weighted.o$sp, funk))
output_df_weighted.o$taxon2 = sub("^(\\w)\\w*\\s(\\w+)$", "\\1. \\2", output_df_weighted.o$taxon)

young_rates = output_df_weighted_young[, c("sp", "r_rate", "r_rate_loCI", "r_rate_hiCI")]
names(young_rates)[2:4] = paste0(names(young_rates)[2:4], "_young")

output_df_combined = merge(output_df_weighted.o, young_rates, by = "sp", all.x = TRUE)

##graph
#drop species with a zero or NA rate in each forest type
output_df_combined_toplot = output_df_combined[!is.na(output_df_combined$r_rate) &
                                                 !is.na(output_df_combined$r_rate_young) &
                                                 output_df_combined$r_rate != 0 &
                                                 output_df_combined$r_rate_young != 0,]
alp = 1

### *Figure S7: giant tree recruitment into old-growth and secondary forest
rec.old =
  ggplot(output_df_combined_toplot, aes(x = r_rate, y = r_rate_young))+
  guides(alpha = "none", size = "none")+
  scale_color_gradient(low= "blue2",high = "red1", name = expression("Giant tree \nmortality"~('%'~yr^{-1})))+
  theme_basis+ 
  scale_x_log10(expression(bold("Old-growth giant tree recruitment rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01), 0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20"),
                limits = c(0.001,0.2))+
  scale_y_log10(expression(bold("Secondary giant tree recruitment rate"~('%'~yr^{-1}))),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01),0.2),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10", "20"),
                limits = c(0.001,0.2))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey50" )+
  geom_errorbar(aes(x = r_rate, y = r_rate_young, ymax = r_rate_hiCI_young, ymin = r_rate_loCI_young), color = "grey50", size = 0.3)+
  geom_errorbar(aes(x = r_rate, y = r_rate_young, xmax = r_rate_hiCI, xmin = r_rate_loCI), color = "grey50", size = 0.3)+
  geom_point(data = output_df_combined_toplot, aes(x = r_rate, y = r_rate_young, color = ML_mortality*100),
             stroke = 0, alpha = alp, size = 3)

#how many species have higher values in secondary than old-growth
higher_secondary = output_df_combined_toplot[output_df_combined_toplot$r_rate_young > output_df_combined_toplot$r_rate,]
nrow(higher_secondary) #count
nrow(higher_secondary)/nrow(output_df_combined_toplot) #proportion of plotted species
higher_secondary$sp

#how many species are significantly different than the 1:1 line?
touches_x = output_df_combined_toplot$r_rate_young >= output_df_combined_toplot$r_rate_loCI &
  output_df_combined_toplot$r_rate_young <= output_df_combined_toplot$r_rate_hiCI
touches_y = output_df_combined_toplot$r_rate >= output_df_combined_toplot$r_rate_loCI_young &
  output_df_combined_toplot$r_rate <= output_df_combined_toplot$r_rate_hiCI_young

output_df_combined_toplot$touches_1to1 = touches_x | touches_y
output_df_combined_toplot$sig_diff = !output_df_combined_toplot$touches_1to1

sum(output_df_combined_toplot$sig_diff) #count significantly different
sum(output_df_combined_toplot$sig_diff)/nrow(output_df_combined_toplot) #proportion of plotted species
output_df_combined_toplot$sp[output_df_combined_toplot$sig_diff]


