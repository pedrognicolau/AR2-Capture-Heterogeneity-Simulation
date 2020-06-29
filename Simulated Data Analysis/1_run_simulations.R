source("Simulated Data Analysis/0_3_functions_for_run_sims_inparallel.R")

#### run in parallel function ####
#setup parallel backend to use many processors
library(doRNG)
library(doParallel)

loopparal <- function(sims=5,params = c(-1,-.5,0.08),ar_scale=20,
                      mn = 30, wsd=8,
                      p1=.55,p2=.75, md = "Mt", seed=123){
  #runs in parallel the function sim_models and returns dataframe
  start_time <- Sys.time()
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload computer
  registerDoParallel(cl)
  finalMatrix <- foreach(i=1:sims, .combine=rbind,.options.RNG=seed) %dorng% { 
    # change path # needs to be full path
    source("~/OneDrive - UiT Office 365/MusData/SCRIPT_FOR_PUBLICATION/Simulated Data Analysis/0_3_functions_for_run_sims_inparallel.R")
    tempMatrix = sim_models(realparams = params, p.d1=p1,p.d2=p2,
                            w.avg.mn=mn,w.avg.sd = wsd, mod=md) #calling a function
    tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  #stop cluster
  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)
  finalMatrix$nsims <- rep(1:sims,each=8)
  return(finalMatrix)
}

try1 <- loopparal(sims=7,params = c(-1,-.5,1),ar_scale=20,
          mn = 30, wsd=8,
          p1=.55,p2=.75, md = "Mt", seed=123)

##### RUN SIMULATIONS IN PARALEL #####
#### run simulations ####
# simloop <- function(phi1=c(-1,-.5,0,.5,1),phi2=c(-.2,-.5,-.8),vars=c(.08,0.2),
#                     pd1 = .55, pd2=.75, mod = "Mt", tsleng = 20, sims = 100,
#                     w.avg.sd=5,inn.sd=1.2,w.avg.mn=30, ar.scale=20){
sims=200
phi1=c(-1,-.5,0,.5,1)
phi2=c(-.2,-.5,-.8)
vars=c(1)
counter=1
sim_var <- list()
  start_time <- Sys.time()
  for(v in vars){  
      print(paste0("var ",v))
      for(ph1 in phi1){
          print(paste0("phi1 ",ph1))
          for(ph2 in phi2){
            print(paste0("phi2 ",ph2))
            print(paste("setting number",counter))
            # sim <- runsims(simulations = sims, realparams = c(ph1,ph2,v),
            #                ts.length = tsleng, ar_scale=ar.scale,
            #                w.avg.mn = w.avg.mn, w.avg.sd=w.avg.sd,inn.sd = inn.sd, p.d1=pd1,
            #                p.d2=pd2, mod = mod)
            sim <- try(loopparal(sims=sims,params = c(ph1,ph2,v),ar_scale=20,
                                  mn = 30, wsd=5,
                                  p1=.55,p2=.75, md = "Mt", seed=counter^1.234))
            # if(typeof(a)!="character")
            sim_var[[counter]] <- sim
            saveRDS(sim_var, "Results/results_list_200_var1.rds")
            counter=counter+1
          }
        }
      }
  end_time <- Sys.time()
  print(paste0("Total running time ",end_time - start_time))

proclist <- summze(sim_var)
check1 <- proclist
str(sim_var)
#### summarize simulations function ####
Bias.2 <- function(est,true) return((as.numeric(est) - true)^2)
true.var <- function(phi1,phi2,inn.var) return(((1-phi2)/(1+phi2))*((inn.var/((1-phi2)^2-phi1^2))))

summze <- function(simloops)
{
  for(s in 1:length(simloops))
  {
    sim.s <- simloops[[s]]
    realparams <- as.numeric(strsplit(sim.s$truear[1], " ")[[1]])
    methods <- unique(sim.s$method)
    
    quantsmin1 <- apply(sim.s[,3:5],1,min); quantsmax1 <- apply(sim.s[,3:5],1,max)
    quantsmin2 <- apply(sim.s[,6:8],1,min); quantsmax2 <- apply(sim.s[,6:8],1,max)
    ar1cond <- realparams[1] >= quantsmin1 & realparams[1] <= quantsmax1
    ar2cond <- realparams[2] >= quantsmin2 & realparams[2] <= quantsmax2
    sim.s$ar1cov <- ifelse(ar1cond,1,0) ; sim.s$ar2cov <- ifelse(ar2cond,1,0)
    sim.s$ar12cov <- ifelse(ar1cond & ar2cond,1,0)
    for(m in methods)
    { 
      sim.m <- sim.s[sim.s$method==m,]
      
      dset <- apply(sim.m[,c(1:9,12:ncol(sim.m))],2,mean.na)
      # RMSE1
      dset$RMSE1 <- sqrt(mean(Bias.2(sim.m$ar1.mean,realparams[1]),na.rm=TRUE))
      # RMSE2
      dset$RMSE2 <- sqrt(mean(Bias.2(sim.m$ar2.mean,realparams[2]),na.rm=TRUE))
      # RMSE12
      dset$RMSE12 <- sqrt(mean(Bias.2(sim.m$ar1.mean,realparams[1])+
                               Bias.2(sim.m$ar2.mean,realparams[2]),na.rm=TRUE))
      # number of simulations
      dset$nsims <- as.numeric(nrow(sim.m[is.na(sim.m$ar1.mean)==FALSE,]))
      dset$truear1 <- realparams[1]; dset$truear2 <- realparams[2]; dset$inn.var <- realparams[3]
      dset$truevar <- true.var(realparams[1],realparams[2],realparams[3])
      dset$method <- m
      
      dset <- as.data.frame(dset)
      if(m==1) chunk <- dset
      else chunk <- rbind(chunk,dset)
    }
    print(s)
    if(s==1) fset <- chunk
    else fset <- rbind(fset,chunk)
  }
  return(fset)
}
check1 <- summze(sim_var)

v04 <- filter(check1,inn.var==.64&method>=5)[,c(1,2,10:15,16:19,21)]
summary(v04)
1.18+2.3+2.4+3+2+2.1+2+1.6+1.8+1.76+1.9+1.78+1.54+1.51+1.4+1.23
30*5*2
300/24
