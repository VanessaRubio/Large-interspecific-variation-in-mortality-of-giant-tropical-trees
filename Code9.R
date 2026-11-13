########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Giant tree mortality and commonly measured functional traits
## Produces: Figures 5a-h
########################################################################################################################

#clear workspace
rm(list=ls()); library(ggplot2); library(readxl); library(magrittr)

#this theme is for plotting
theme_basis = theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text = element_text(family = "Arial", colour="black", face ="bold", size = 16))

load("ML_mortality_50CI.Rdata")
data = ML_mortality

#add functional trait data from Wright et al. 2010
traits = read.table("BCITRAITS_20101220.txt", h = T)
traits$SP. = tolower(traits$SP.)
rownames(traits) = traits$SP.

funk= function(x){ 
  use = match(toString(x), traits$SP.)
  #seed dry mass, leaf thickness, leaf area, leaf mass per area, fracture toughness, maximum height, leaf dry matter content
  t= traits[use,c("SG100C_AVG","SEED_DRY","LEAFTHCK_AVD","LEAFAREA_AVD","LMALEAF_AVD","AVG_LAMTUF","HEIGHT_AVG", "LDMC_AVD")]
  return(t)
}

all.t = do.call("rbind", lapply(data$sp, funk))
data = cbind(data,all.t)

#add leaf chemical data - Wright
chem = as.data.frame(readxl::read_excel("leaf.chemical.xlsx", col_names = T))
chem$SPECIES = tolower(chem$SPECIES)
chem$LFP_SUN = as.numeric(chem[,"LFP_SUN"])
chem$LFN_SUN = as.numeric(chem[,"LFN_SUN"])

funk= function(x){ 
  use = match(toString(x), chem$SPECIES)
  #leaf phosphorus, leaf nitrogen
  c = chem[use,c("LFP_SUN","LFN_SUN")]
  return(c)
}

all.c = do.call("rbind", lapply(data$sp, funk))
data = cbind(data,all.c) #dim 39,19

#select only the species that have all the traits complete
#transform traits
shapiro.test(data$SG100C_AVG) #wood density, p-value = 0.3

shapiro.test(data$SEED_DRY) #seed mass, p-value = 0.00000000004
data$SEED_DRY=log(data$SEED_DRY)

shapiro.test(data$LEAFTHCK_AVD) #leaf thickness, p-value = 0.0001
data$LEAFTHCK_AVD=log(data$LEAFTHCK_AVD)

shapiro.test(data$LEAFAREA_AVD) #leaf area, p-value = 0.000000004
data$LEAFAREA_AVD=log(data$LEAFAREA_AVD)

shapiro.test(data$LMALEAF_AVD) #leaf mass per area, p-value = 0.00003
data$LMALEAF_AVD=log(data$LMALEAF_AVD)

shapiro.test(data$AVG_LAMTUF) #leaf toughness, p-value = 0.05825
shapiro.test(data$HEIGHT_AVG) #max.height, p-value = 0.3
shapiro.test(data$LDMC_AVD) #leaf dry matter content, p-value = 0.9
shapiro.test(as.numeric(data$LFP_SUN)) #leaf phosphorus concentration, p-value = 0.6
shapiro.test(as.numeric(data$LFN_SUN)) #leaf nitrogen concentration, p-value = 0.08

#graphs of individual traits
data.pca$ML.rate = data.pca$ML_mortality
ps = 3

### *Figure 5a: giant tree mortality and max height
m1 = lm(log(data$ML_mortality) ~ data$HEIGHT_AVG)
summary(m1)
maxh = ggplot(data = data, aes(x = HEIGHT_AVG, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  labs(x = "Maximum height (m)") + geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "solid", col = "grey50" )

### *Figure 5b: giant tree mortality and wood density
wd = ggplot(data = data, aes(x = SG100C_AVG, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Wood density"~('g'~cm^{-3}))))

### *Figure 5c: giant tree mortality and seed mass
sm = ggplot(data = data, aes(x = SEED_DRY, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  labs(x = "Seed mass (g)") 

### *Figure 5d: giant tree mortality and leaf mass per area
lma = ggplot(data = data, aes(x = LMALEAF_AVD, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Leaf mass per area"~('g'~m^{-2}))))

### *Figure 5e: giant tree mortality and leaf dry matter content
ldmc = ggplot(data = data, aes(x = LDMC_AVD, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Leaf dry matter content"~('g'~g^{-1}))))

### *Figure 5f: giant tree mortality and leaf phosphorus concentration
leafP = ggplot(data = data, aes(x = LFP_SUN, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Leaf P content"~('g'~g^{-1}))))

### *Figure 5g: giant tree mortality and leaf nitrogen concentration
leafN = ggplot(data = data, aes(x = LFN_SUN, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Leaf N content"~('g'~g^{-1}))))

### *Figure 5h: giant tree mortality and leaf area
la = ggplot(data = data, aes(x = LEAFAREA_AVD, y = ML_mortality))+ geom_point(size = ps) +  theme_basis+
  scale_y_log10(expression(bold("Giant tree mortality rate"~('%'~yr^{-1}))),
                limits = c(0.003, 0.1),
                breaks = c(seq(0.001,0.01,0.001),seq(0.02,0.1,.01)),
                labels = c("0.1",rep("",8),"1.0",rep("",8),"10")) + 
  scale_x_continuous(expression(bold("Leaf area"~(cm^{2}))))
