#### FUNCTIONS TO RUN SIMULATION ####
library(dplyr)
library(INLA)
library(VGAM)
library(dplyr)

source("Simulated Data Analysis/0_1_functions_data_generation.R")
source("Simulated Data Analysis/0_2_functions_for_data_analysis.R")

#### NECESSARY FUNCTIONS / OBJECTS ####
# compute period based on phi1 (x1) and phi2(x2)
periodAR <- function(x1,x2) 2*pi/(acos(x1/(2*sqrt(-x2))))

# when the simulation crashes, replace with
failframe <- data.frame(ar1.mean = rep(NA,8),ar2.mean = rep(NA,8),
                        ar1.025 = rep(NA,8),ar1.50 = rep(NA,8),ar1.975 = rep(NA,8),
                        ar2.025 = rep(NA,8),
                        ar2.50 = rep(NA,8),
                        ar2.975 = rep(NA,8), 
                        ar.var=rep(NA,8))

mean.na <- function(x) 
  # mean function to use in apply setups
{
  m <- mean(as.numeric(x),na.rm = TRUE)
  return(m)
}


#### SINGLE SIMULATION FUNCTION ####
sim_models <- function(realparams = c(-1,-.4,1), 
                       ts.length = 20, ar_scale=25,
                       w.avg.mn = 30, w.avg.sd=8,inn.sd = 1.2, p.d1=.55,p.d2=.75, mod = "Mt")
{
  
  proceed = FALSE;  failures = 0
  while(proceed == FALSE){
    tryCatch({
      failures=failures+1
      sets <- generate.cr.data(params_ar = realparams, ar.scale = ar_scale,
                           timepoints = ts.length,model=mod,
                           w.avg.mn = w.avg.mn, w.avg.sd=w.avg.sd,inn.sd = inn.sd, p.d1=p.d1,p.d2=p.d2, plot.d=FALSE)
      
      
      tryset <- sets[[1]] # counts
      tryset2 <- sets[[2]] # cr data
      
      methods <- c("BASELINE A","CR-INLA A","CR-VGAM A","ObsCount A","BASELINE P","CR-INLA P", "CR-VGAM P", "ObsCount P")
      
      print(methods[1])
      rowdf <- fit.ar(tryset, truelog=TRUE)
      
      print(methods[2])
      rowdf  <- rbind(rowdf,fit.ar(ht.estimates(tryset2, timepoints=ts.length)))
      
      print(methods[3])
      rowdf  <- rbind(rowdf,fit.ar(vgam.ht(vgam.fitted(tryset2, model = mod), timepoints=ts.length)))
      
      print(methods[4])
      rowdf  <- rbind(rowdf,fit.ar(tryset)) 
      
      print(methods[5])
      rowdf  <- rbind(rowdf,ar2.poisson.rate(tryset,Lhood = "poisson", trueN=TRUE))
      
      print(methods[6])
      rowdf  <- rbind(rowdf,ar2.poisson.rate(ht.estimates(tryset2), Lhood = "poisson"))
      
      print(methods[7])
      rowdf  <- rbind(rowdf,ar2.poisson.rate(vgam.ht(vgam.fitted(tryset2, model = mod), timepoints=ts.length),Lhood = "poisson"))
      
      print(methods[8])
      
      
      rowdf  <- rbind(rowdf,ar2.poisson.rate(tryset,Lhood = "poisson"))
      
      proceed = TRUE},error=function(e){
        cat(paste("ERROR ",failures," \n"))})
    
    if(failures==8) {proceed = TRUE; rowdf <- failframe; print("break loop")}
  }
  
  #### if failure happens ####
  
  # if an error occurred, invalidate entire line
  # inefficient but should do for now
  # if (failures==8) rowdf <- failframe
  
  rowdf$method <- 1:8
  rowdf$truear <- paste(realparams[1],realparams[2],realparams[3])
  rownames(rowdf) <- NULL
  return(rowdf)
}
