### GENERATE CAPTURE HISTORIES FUNCTIONS ####

which.cat <- function(row) {return(which(row==1))}

# true process variance function
true.var <- function(phi1,phi2,inn.var) return(
                                               (1-phi2)*inn.var/
                                                ((1+phi2)*((1-phi2)^2-phi1^2))
                                               )
# theoretical mean value of the process
theoreticalK <- function(N, phi1,phi2,inn.var) return(log(N)-1/2*true.var(phi1,phi2,inn.var))
# log K = log N - 1/2 * sigma^2

sig.ch <- function(timepoints = 20, params_ar = c(.5,-.5,.04), ar.scale=NA, p.d1=.55,p.d2=.75,
                   w.avg.mn = 15, w.avg.sd =5, model = "Mt", inn.sd = 1.1)
  # returns simulated dataset with parameters defined above where
  # timepoints = length of time series
  # params_ar = phi1, phi2, var
  # ar.scale = multiplication factor for the ar series
  # p.d1 = mean probability of capture on day one
  # p.d2 = mean probability of capture on day two
  # weight.avg.mn = mean of the weight distribution for each timepoint
  # weight.avg.sd = sd of the weight distribution for each timepoint
  # w.avg.mn = 15
  # w.avg.sd = 5
  # inn.sd = sd parameter for individual weight generation
  # model = one of Otis models, including Mth or Mbh
  
{
  ds <- arima.sim(model=list(ar=c(params_ar[1],params_ar[2])),sd=sqrt(params_ar[3]),n=timepoints)
  # find offset
  offset=theoreticalK(N=ar.scale,phi1=params_ar[1],phi2=params_ar[2],inn.var=params_ar[3])
  # round up time series
  ds1 <- round(exp(ds+offset))
  # number of individuals
  no.indiv <- sum(ds1)
  
  # Capture probabilities average for the categories: (1,1), (1,0), (0,1), (0,0)
  if (model == "Mb") prob=c(p.d1*p.d2,p.d1*(1-p.d2),(1-p.d1)*p.d1,(1-p.d1)^2)
  if (model == "Mt") prob=c(p.d1*p.d2,p.d1*(1-p.d2),(1-p.d1)*p.d2,(1-p.d1)*(1-p.d2)) # Mth
  
  # Differences of the alternative specific coefficients are log of prob ratios
  gamma.diff1=log(prob[1]/prob[2])
  if (model == "Mt") gamma.diff2=log(prob[1]/prob[3])
  if (model == "Mb") gamma.diff2=log(prob[3]/prob[1])
  
  #generate realisitic weights
  # sampling weight means for the different timepoints
  weight.avg <- rlnorm(timepoints,log(w.avg.mn),log(w.avg.sd))
  # attribute a weight mean to a given individual belonging to a certain timepoint
  wvec <- rep(weight.avg,times=ds1)
  # generate weight for the individual with previously generated mean and sd = inn.sd
  weight <- rlnorm(no.indiv,log(wvec),log(inn.sd))
  # scale weight
  weight=(weight)/sd(weight)

  # compute individual probabilities based on ratios
  a1=exp(gamma.diff1*weight)
  a2=exp(gamma.diff2*weight)
  p2=a1/(1+a1)
  # par(mfrow=c(1,2))

  # which Otis model to choose from
  if (model == "Mt"){p1=a2/(1+a2)}
  if (model == "Mb"){p1=1-p2*exp(gamma.diff2*weight)}
  
  # create matrix with individual probabilities
  prob.mat.cr=matrix(nrow=4,ncol=no.indiv)
  
  # Mo
  if (model == "Mb"){
  prob.mat.cr[1,]=p1*p2        #11
  prob.mat.cr[2,]=p1*(1-p2)    #10
  prob.mat.cr[3,]=(1-p1)*p1    #01
  prob.mat.cr[4,]=(1-p1)^2     #00
  }    
  
  if(model == "Mt"){
  prob.mat.cr[1,]=p1*p2             #11
  prob.mat.cr[2,]=p1*(1-p2)         #10
  prob.mat.cr[3,]=(1-p1)*p2         #01
  prob.mat.cr[4,]=(1-p1)*(1-p2)     #00
  } 
  
  ### Generate categories for each individual ###
  x=matrix(nrow=no.indiv,ncol=4)
  for (j in 1:no.indiv) x[j,]=t(rmultinom(1,1,prob.mat.cr[,j]))
  #retrieve category
  cats <- apply(x,1,which.cat)
  catlab <- c("1,1","1,0","0,1","0,0")
  # attach labels 
  cat.nr <- sapply(cats,function(x)return(catlab[x]))
  
  # build data frame
  capdata <- data.frame(id=1:no.indiv,weight=weight,p1=p1,p2=p2,chist=cat.nr)
  # add time point for each ind
  capdata$timepoint <- rep(1:timepoints,times=ds1)
  # add individual variable c1 and c2
  splic <- as.numeric(unlist(strsplit(cat.nr,",")))
  capdata$c1 <- splic[seq(1,length(splic)-1,2)]
  capdata$c2 <- splic[seq(2,length(splic),2)]
  
  # Number of individual
  Nframe=data.frame(Ntotal=ds1,logN=ds)
  
  return(list(capdata,Nframe))
}

# crdata <- capdata
generate.cr.data <- function(timepoints = 20, params_ar = c(.5,-.5,.08),
                       ar.scale = 10, p.d1=.55,p.d2=.75,
                       w.avg.mn = 30, w.avg.sd = 5, inn.sd=1.2,
                       plot.d=FALSE, model = "Mt")
  # return capture history dataset without captured individuals 
  # and count dataset associated with it
{ #
  #create N's for different timepoints and locations
  crdatalist <- sig.ch(timepoints=timepoints, params_ar=params_ar, ar.scale=ar.scale,w.avg.mn = w.avg.mn, w.avg.sd = w.avg.sd,
                   inn.sd = inn.sd, model = model, p.d1=p.d1,p.d2=p.d2)
  crdata=crdatalist[[1]]
  truenumbers=crdatalist[[2]]
  # remove no captures
  crdata_nozeros <- crdata[crdata$chist!="0,0",]
  freq00 <- 1-(nrow(crdata_nozeros)/nrow(crdata))
  print(paste("percentage of non-caps =",(1-nrow(crdata_nozeros)/nrow(crdata))))
  print(paste("Real N=",nrow(crdata)))
  
  ### raw count dataset ###
    crdata3 <- crdata_nozeros
    crdata3$timepoint <- factor(crdata3$timepoint, levels=c(1:timepoints))
    crdata3$count <- 1; 
    justcounts <- as.data.frame(xtabs(count ~ timepoint, data = crdata3))
    justcounts$Ntotal <- truenumbers$Ntotal
    justcounts$truelog <- truenumbers$logN
    
    colnames(justcounts)[2] <- "count"

  ## capture history ##
    crdata2 <- crdata_nozeros
    ## obtain capture history variable
    ### TRANSFORM DATASET FOR INLA ###
    ## add count variable
    crdata2$count <- 1
    crdata2$id <- 1:nrow(crdata2)
    ## copy dataset
    crdata3 <- crdata2
    ## add count 0
    crdata3$count <- 0
    ## copy again
    crdata4 <- crdata3
    ## on the new copied set 1, invert all the labels
    crdata3$chist <- as.factor(ifelse(crdata2$chist == "1,1", "1,0", ifelse(crdata2$chist == "1,0", 
                                                                            "0,1", "1,1")))
    ## on the new copied set 2, invert all the labels again
    crdata4$chist <- as.factor(ifelse(crdata3$chist == "1,1", "1,0", ifelse(crdata3$chist == "1,0", 
                                                                            "0,1", "1,1")))
    
    ## final set involves all combinations of observed categories and non-observed
    inladata <- rbind(crdata2,crdata3,crdata4)
    
    inladata$alt.id <- ifelse(inladata$chist=="0,1",1,ifelse(inladata$chist=="1,0",2,3))
    inladata$alt.id2 <- inladata$alt.id
    
    inladata2 <- inladata
    inladata3 <- arrange(inladata2,id,alt.id)
    
    # add id term with NAs
    idre <- rep(NA,nrow(inladata3))
    idre[seq(1,length(idre),3)]<-1:(length(idre)/3)
    inladata3$id_re <- idre
    #add proportion of non-captures just to store it
    inladata3$freq00 <- freq00
    # add timepoint formatted for inla
    inladata3$timepoint.i <- NA
    inladata3$timepoint.i[seq(1,nrow(inladata3),3)] <-inladata3$timepoint[seq(1,nrow(inladata3),3)]
    data4inla <- inladata3
    result=list(justcounts, inladata3)
  return(result)
  
} 