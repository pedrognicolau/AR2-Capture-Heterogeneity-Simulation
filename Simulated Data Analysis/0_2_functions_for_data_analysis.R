### FUNCTIONS TO FIT MODELS IN SIMULATION SETUP ###
require(dplyr)
require(INLA)
require(VGAM)

logna <- function(x) ifelse(x==0,0,log(x))

# Standardize inla fitted values to add up to one
stdize <- function(matrix){
  #set counter
  fit2 <- matrix
  j=1
  for(i in seq(1,length(fit2),3)){ #matrix needs to be made of triplets
    seqn <- i:(i+2)
    total <- sum(fit2[seqn])
    
    for(j in seqn) fit2[j]<- fit2[j]/total
    j=j+3
  }
  return(fit2)
}

# Transform triplet sequence into data frame with 3 columns
# creates data frame from sequence of probability values
makedf <- function(df){
  for(i in seq(1,length(df),3))
  {
    if (i == 1) {A <- c(df[i]); B <- c(df[i+1]); C <- c(df[i+2])}
    else {A <- c(A, df[i]); B <- c(B,df[i+1]); C <- c(C,df[i+2])}
  }
  dataf <- data.frame("0,1"=A,"1,0"=B,"1,1"=C)
  return(dataf)
}

# Horvitz-Thompson Estimator
HT <- function(x) sum(1/(1-x)) #x is the probability of 0,0

partialN <- function(p00.st.tp){
  df <- p00.st.tp
  colnames(df) <- c("p00","st","tp")
  
  newdata <- distinct(p00.st.tp[,2:3]) ; colnames(newdata) <- c("station","timepoint")
  newdata$N.est <- NA
  
  sts <- sort(unique(df$st))
  c <- 1
  for(s in 1:length(sts)){
    # print(sts[s])
    df_st <- filter(df, st == sts[s])
    
    for (timep in sort(unique(df_st$tp)))
    {
      df_tp <- filter(df_st, tp == timep)
      nt <- HT(df_tp$p00)
      
      newdata$N.est[newdata$station==sts[s]&newdata$timepoint==timep] <- nt
    }
  }
  
  return(newdata)
}

# Computes capture probability using weight as a single covariate and gammas (see paper for definition)
prob.comp <- function(gama1,gama2,weight) (exp((gama1-gama2)*weight))

ht.estimates <- function(ch.dataset, p.equation=TRUE, model = "Mt", timepoints = 20)
  # Fits inla model to capture history data and return N
  # formula includes id, alt.id, weight.s, alt.id2, sex.new, station.i, timepoint.i
  #
  # allows two ways of computing probabilies
  # one through computing gammas described in paper (p.equation=TRUE)
  # another solving the equation system using the estimates of the conditional capture history probabilities (p.equation=FALSE)
  #
{
  #formula for inla model
  formula <- 
    count ~ f(id, initial = -10, fixed=T) + f(alt.id, weight, fixed = T, constr = T) + 
    f(timepoint.i, model = "rw1", scale.model = TRUE, hyper = list(prec=list(prior="pc.prec", 
                                                                             param = c(U=3,alpha=0.01))))
  
  #fit model
  mod <- inla(formula, family = "Poisson", data = ch.dataset, control.predictor = list(compute=T),
              control.compute = list(dic=T))

  # try Monte Carlo estimates 
  # retrieve fitted values
  fitval <- mod$summary.fitted.values$mean
  #turn them into a data frame
  prob_df <- makedf(stdize(fitval))
  
  rawdata <- filter(ch.dataset, count == 1)
  
  if(p.equation==TRUE){
    
    # RETRIEVE THE GAMMAS
    g01 <- mod$summary.random$alt.id$mean[1]
    g10 <- mod$summary.random$alt.id$mean[2]
    g11 <- mod$summary.random$alt.id$mean[3]
    
    prob_df$p2_inla <- prob.comp(g11,g10,rawdata$weight)/(1+prob.comp(g11,g10,rawdata$weight))
    
    # MODEL WITH TIME
    if (model == "Mt"){prob_df$p1_inla <- prob.comp(g11,g01,rawdata$weight)/(1+prob.comp(g11,g01,rawdata$weight))}
    
    # MODEL WITH BEHAVIOUR
    if (model == "Mb"){prob_df$p1_inla <- 1- prob_df$p2_inla*prob.comp(g01,g11,rawdata$weight)}
  }
  #calculate p1 and p2 based on Mt
  
  if(p.equation==FALSE){
    alpha01 <- prob_df$X1.1/prob_df$X0.1
    alpha10 <- prob_df$X1.1/prob_df$X1.0
    prob_df$p2_inla <- (alpha10)/(1+alpha10)
    prob_df$p1_inla <- (alpha01-prob_df$p2_inla)/(alpha01)
  }
  
  #NON-CAPTURE PROBABILITY
  prob_df$p00_inla <- (1-prob_df$p1_inla)*(1-prob_df$p2_inla)
  
  #lump datasets
  totaldf <- data.frame(rawdata,p1.est=prob_df$p1_inla,
                        p2.est=prob_df$p2_inla, p00.est=prob_df$p00_inla,
                        p00true = (1-rawdata$p1)*(1-rawdata$p2))
  
  #obtaining time series
  ts <- arrange(partialN(p00.st.tp=data.frame(totaldf$p00.est,1,totaldf$timepoint)),
                timepoint)
  # add zeros
  ts$timepoint <- factor(ts$timepoint, levels=c(1:timepoints))
  adjustedcounts <- as.data.frame.table(xtabs(N.est~timepoint, data=ts))
  
  colnames(adjustedcounts)[2] <- "count"
  return(adjustedcounts)
  
}

fit.inla.ar <- function(tsdata = NA)
  # fits ar model to variable logcount 
  # needs to have vars logcount and timepoint
{
  
  modi <- inla(logcount ~ f(timepoint, model = "ar", order = 2,
                            hyper = list(prec = list(prior = "pc.prec",
                                                     param = c(U = 3, alpha = 0.01))
                                         ,
                                         pacf1 = list(prior = "pc.cor0",
                                                      param = c(U=0.5,alpha=0.5)),
                                         pacf2 = list(prior = "pc.cor0",
                                                      param = c(U=0.5,alpha=0.5)))),
               control.inla = list(h=.001),
               control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
               data = tsdata)
  return (modi)
  
}

fit.ar <- function(timeseries, trueN=FALSE, truelog=FALSE)
  # receives timeseries dataframe and uses fit.inla.ar function, to retrieve a row containing all necessary information
{
  # no true log
  if(truelog == FALSE){
    #if true N
    if(trueN==TRUE) {timeseries$count <- timeseries$Ntotal}
    timeseries$logcount <- log(timeseries$count+1)
  }
  else {colnames(timeseries)[4] <- "logcount"}
  
    timeseries$id <- 1:nrow(timeseries)
  # fit AR model 
    modi <- fit.inla.ar(timeseries)

    # MC estimation
    modi.hsample <- as.data.frame(inla.hyperpar.sample(2000,modi)) #inla
    
    modi.hsample$phi1 <- modi.hsample$`PACF1 for timepoint`*(1-modi.hsample$`PACF2 for timepoint`) #inla
    ar1q <- quantile(modi.hsample$phi1, c(0.025,0.5,0.975)) #inla
    
    ars <- data.frame(ar1.mean = mean(modi.hsample$phi1),ar2.mean = modi$summary.hyperpar$mean[3],
                      ar1.025 = ar1q[1],ar1.50 = ar1q[2],ar1.975 = ar1q[3],
                      ar2.025 = modi$summary.hyperpar$`0.025quant`[3],
                      ar2.50 = modi$summary.hyperpar$`0.5quant`[3],
                      ar2.975 = modi$summary.hyperpar$`0.975quant`[3], 
                      ar.var=1/modi$summary.hyperpar$mean[1])
  return(ars)
  
}

ar2.poisson.rate <- function(countdf, Lhood = "poisson", trueN=FALSE) # "nbinomial"
  # receives count data.frame with count and timepoint and returns ar coefficients for poisson rate
  #an ar2 in a data.frame
  # if inla.ar == TRUE it fits an inla ar2 and returns credible intervals 
{

  if(trueN==TRUE) {countdf$count<-countdf$Ntotal} 
  
  form.pois <- count ~ f(timepoint, model = "iid") #}
  
  #scale.model makes sure the smoothing doesnt depend on the length of the time series
  modpois <- inla(form.pois, family = Lhood, data = countdf, control.predictor = list(compute=T))
  
  count.df <- data.frame(logcount=modpois$summary.linear.predictor$mean, timepoint=1:nrow(countdf))

  # fit ar model to poisson rate
  modi <- fit.inla.ar(count.df) #
    
    # fixes gaussian precision
    
    modi.hsample <- as.data.frame(inla.hyperpar.sample(2000,modi)) #inla
    modi.hsample$phi1 <- modi.hsample$`PACF1 for timepoint`*(1-modi.hsample$`PACF2 for timepoint`) #inla
    ar1q <- quantile(modi.hsample$phi1, c(0.025,0.5,0.975)) #inla
    
    ars <- data.frame(ar1.mean = mean(modi.hsample$phi1),ar2.mean = modi$summary.hyperpar$mean[3],
                      ar1.025 = ar1q[1],ar1.50 = ar1q[2],ar1.975 = ar1q[3],
                      ar2.025 = modi$summary.hyperpar$`0.025quant`[3],
                      ar2.50 = modi$summary.hyperpar$`0.5quant`[3],
                      ar2.975 = modi$summary.hyperpar$`0.975quant`[3], 
                      ar.var=1/modi$summary.hyperpar$mean[1])
  
  return(ars)
  
}

vgam.fitted <- function(inla.dataset, model = "Mt")
{
  #fits vgam model and retrieves fitted values to an inla dataset
  crdata <- filter(inla.dataset, count == 1 )
  
  if(model == "Mt") { print("Mt"); mod <- vglm(cbind(c1,c2) ~ weight, posbernoulli.t, data = crdata)}
  if(model == "Mb"){ print("Mb") ; mod <- vglm(cbind(c1,c2) ~ weight, posbernoulli.b, data = crdata) }
  if(model == "Mtb") { print("Mtb") ; mod <- vglm(cbind(c1,c2) ~ weight, posbernoulli.tb, data = crdata)}
  
  chist2 <- as.data.frame(cbind(timepoint=crdata$timepoint, 
                                totalN = crdata$Ntotal,
                                true.p1=crdata$p1,true.p2=crdata$p2,fitted(mod)))

  colnames(chist2)[(ncol(chist2)-1):(ncol(chist2))] <- c("p1","p2")
  
  return(chist2)
  
}

vgam.ht <- function(vgam.fitted, timepoints = NA)
  # uses HT estimator on VGAM estimated capture probabilities
{
  totaldf <- vgam.fitted
  totaldf$p00.est <- (1-totaldf$p1)*(1-totaldf$p2)
  
  ts <- arrange(partialN(p00.st.tp=data.frame(totaldf$p00.est,1,totaldf$timepoint)),
                timepoint)
  # add zeros
  ts$timepoint <- factor(ts$timepoint, levels=c(1:timepoints))
  adjustedcounts <- as.data.frame.table(xtabs(N.est~timepoint, data=ts))
  colnames(adjustedcounts)[2] <- "count"
  return(adjustedcounts)
  
}

# computing AR period based on coefficients
periodAR <- function(x1,x2) 2*pi/(acos(x1/(2*sqrt(-x2))))

mean.na <- function(x) 
  # mean function to use in apply() settings

{
  m <- mean(x,na.rm = TRUE)
  return(m)
}

# 

