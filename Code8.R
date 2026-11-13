########################################################################################################################
### This code is part of the publication: Large interspecific variation in mortality of giant tropical trees: magnitude, importance, and ecology
## Expected R2 between giant-tree mortality and small size-class mortality if the relationship were perfect (R2=1), with sampling error on both axes.
## Produces: Figures S4
########################################################################################################################

#clear workspace
rm(list=ls())
library("dplyr")
library("lmodel2")

#observed rates giant tree mortality and size-class mortality
load("ML_mortality_50CI.Rdata");  m.giant = ML_mortality
load("ML_mortality_1_10.Rdata"); m.1.10 = ML_mortality
load("ML_mortality_10_50.Rdata"); m.10.50 = ML_mortality
load("ML_mortality_1_50.Rdata"); m.1.50 = ML_mortality

#add the giant-tree rate to each size-class table, matched by species
#(1-10cm has 38 of the 39 species, no anacex)
add.giant = function(small){
  idx = match(as.character(small$sp), as.character(m.giant$sp))
  stopifnot(!anyNA(idx))
  small$ML50 = m.giant$ML_mortality[idx]
  small
}
m.1.10 = add.giant(m.1.10)
m.10.50 = add.giant(m.10.50)
m.1.50 = add.giant(m.1.50)

#census-interval data
load("survival_data_assembled.Rdata")
simple_combined = simple_combined.vr

#drop "*", "", and "unk?"
simple_combined_sp = simple_combined[simple_combined$sp!="*" & simple_combined$sp!="" &
                                       simple_combined$sp!="unk?" & !is.na(simple_combined$sp) &
                                       !is.na(simple_combined$dbhmm.1.3),]
simple_combined_sp$interval = simple_combined_sp$endcensus - simple_combined_sp$opencensus

#intervals of one size class (dbh in mm, [min, max)), only for the species in spp
get.intervals = function(dbh.min, dbh.max, spp){
  d = simple_combined_sp[simple_combined_sp$dbhmm.1.3 >= dbh.min &
                           simple_combined_sp$dbhmm.1.3 <  dbh.max &
                           simple_combined_sp$sp %in% spp,
                         c("plot_treeID", "sp", "interval", "mortality")]
  d = as.data.frame(d)
  d = d[d$plot_treeID != "plotKnight_3338",] #unconsistence data
  n.bad = sum(!(d$interval > 0))
  if(n.bad > 0) cat(n.bad, "intervals <= 0 or NA removed\n")
  d[d$interval > 0 & !is.na(d$interval),]
}

int.giant = get.intervals(500, Inf, m.giant$sp)   #>=50cm
int.1.10 = get.intervals(10, 100, m.1.10$sp)    #1-10cm
int.10.50 = get.intervals(100, 500, m.10.50$sp)   #10-50cm
int.1.50 = get.intervals(10, 500, m.1.50$sp)    #1-50cm

#ML mortality rate (Kubo et al. 2000), using uniroot (faster) intead of minimizations as in Code1.R
ml.rate = function(dead, living){
  if(length(dead) == 0) return(NA_real_)
  S = sum(living)
  f = function(r) -S + sum(dead / expm1(r * dead))
  if(f(1) >= 0) return(1)
  uniroot(f, c(1e-12, 1), tol = 1e-10)$root
}

#rates for all species from a vector of 0/1 deaths
rates.by.sp = function(int.list, death.list){
  mapply(function(iv, dd) ml.rate(iv[dd == 1], iv[dd == 0]), int.list, death.list)
}

#simulation
#small is size-class table (m.1.10, m.10.50 or m.1.50)
#interval.size.class is the census intervals of that size class (int.1.10, int.10.50 or int.1.50)
#main title for the histogram
#draw.x is TRUE to simulate both axes 
simulate.R2.both = function(small, interval.size.class, main, n.draws = 1000, draw.x = TRUE, seed = 1){
  spp = as.character(small$sp)
  
  #observed fit and the giant-tree rates lying on the MA line (R2 = 1)
  obs = suppressMessages(lmodel2(log(small$ML50) ~ log(small$ML_mortality), nperm = 99))
  intercept = obs$regression.results$Intercept[2]  
  slope = obs$regression.results$Slope[2]      
  y.true = exp(intercept + slope * log(small$ML_mortality))
  x.true = small$ML_mortality
  
  #intervals split by species, in the same order as small size data
  iv.y = split(int.giant$interval, factor(int.giant$sp, levels = spp))
  iv.x = split(interval.size.class$interval, factor(interval.size.class$sp, levels = spp))
  stopifnot(all(lengths(iv.y) > 0), all(lengths(iv.x) > 0))
  
  #probability of dying in each interval
  p.y = mapply(function(iv, r) -expm1(-r * iv), iv.y, y.true, SIMPLIFY = FALSE)
  p.x = mapply(function(iv, r) -expm1(-r * iv), iv.x, x.true, SIMPLIFY = FALSE)
  
  set.seed(seed)
  sim.y = matrix(NA_real_, nrow = length(spp), ncol = n.draws,
                 dimnames = list(spp, paste0("draw_", 1:n.draws)))
  sim.x = sim.y
  r2    = rep(NA_real_, n.draws)
  n.sp  = rep(NA_integer_, n.draws)
  
  for(k in 1:n.draws){
    if(k %% 50 == 0) cat(main, "draw", k, "\n")
    
    #simulated deaths (0/1) for every interval, then re-estimated rates
    sim.y[, k] = rates.by.sp(iv.y, lapply(p.y, function(p) rbinom(length(p), 1, p)))
    sim.x[, k] = if(draw.x) rates.by.sp(iv.x, lapply(p.x, function(p) rbinom(length(p), 1, p))) else x.true
    
    #species with no simulated deaths (NA) or a 0 rate are left out of this draw
    ok = is.finite(sim.y[, k]) & sim.y[, k] > 0 & is.finite(sim.x[, k]) & sim.x[, k] > 0
    fit = suppressMessages(suppressWarnings(
      lmodel2(log(sim.y[ok, k]) ~ log(sim.x[ok, k]), nperm = 0)))
    r2[k] = fit$rsquare
    n.sp[k] = sum(ok)
  }
  
  r2.obs = obs$rsquare
  ci = quantile(r2, c(0.025, 0.975), na.rm = TRUE)

  hist(r2, breaks = seq(0, 1, by = 0.015), col = "#00BFC4", border = "white",
       main = main, xlab = expression(R^2), ylab = "Frequency", las = 1,
       xlim = c(0, 1), xaxs = "i", xaxt = "n")   
  axis(1, at = seq(0, 1, by = 0.2))              
  abline(v = r2.obs,                 col = "#FF68A1", lwd = 2,   lty = 3)  
  abline(v = mean(r2, na.rm = TRUE), col = "red",     lwd = 2,   lty = 1)  
}

### *Figure S4: simulated R2 for all size-classes
par(mfrow = c(2, 3))
both.1.50  = simulate.R2.both(m.1.50,  int.1.50,  "(>=50cm ~ 1-50cm)")
both.1.10  = simulate.R2.both(m.1.10,  int.1.10,  "(>=50cm ~ 1-10cm)")
both.10.50 = simulate.R2.both(m.10.50, int.10.50, "(>=50cm ~ 10-50cm)")
