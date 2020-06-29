#### plot phi1 vs phi2 ####
setwd("/Users/pni026/OneDrive - UiT Office 365/")
dir("Results")
finalres <- readRDS("Results/final_results_list_200.rds")
check1 <- summze(finalres) # source from runsimulations
ttt2 <- check1
head(ttt2,50)
ttt2$intar1 <- ttt2$ar1.975-ttt2$ar1.025
ttt2$intar2 <- ttt2$ar2.975-ttt2$ar2.025

hist(ttt2$ar.var[ttt2$method==1], breaks=seq(0,2,.1), col=alpha("red",.5))
hist(ttt2$ar.var[ttt2$method==2], breaks=seq(0,2,.1), col=alpha("orange",.5), add=TRUE)

boxplot(sqrt(ttt2$ar.var)~ttt2$method, breaks=seq(0,2,.1), col=c("green3","orangered2","dodgerblue2","gray90"), names=c("0A","1A","2A",
                                                                                                                        "3A",
                                                                                                                        "0P","1P","2P",
                                                                                                                        "3P"))
par(mfrow=c(1,2))
boxplot(sqrt(ttt2$intar2)~ttt2$method, breaks=seq(0,2,.1), col=cols, 
        names=c("0A","1A","2A","3A","0P","1P","2P","3P"), xlab=expression(phi[2]),ylab="AR 95% Cred Interval Span", ylim = c(.6,1.3))

fit.ar(countdata)
modi <- fit.inla.ar(timeseries)

cols = c("green2","red","deepskyblue2","gray50")
summary(check1)
par(mar=c(5, 4, 4, 2) + 0.1)
mlab<-c("trueN", "ht", "vgam","rcounts","trueN poisson", "ht poisson", "vgam poisson","rcounts poisson")
par(mfrow=c(1,2))
# ratio has to be at least 7 to 12 for poisson next to raw

cx=1.4
tr=.8

for(v in c(.04,.08,.16,.32,.64)){
  ## ABUNDANCE FIGURE ##
  setEPS()
  postscript(paste0("/Users/pni026/OneDrive - UiT Office 365/Results/p1vp2_a_v",round(v*100),".eps"),width = 6, height = 6)
  par(mar=c(5, 5, 4, 2) + 0.1, cex.lab=2)
  plot(1, type="n", ylim=c(-.9, .1), xlim=c(-1,1), 
       xlab = expression(phi["1"]), ylab = expression(phi["2"]), cex = cx, main="")
  
  for(i in 4:1)
  {
    
    filtset <- filter(check1, method==i & inn.var == v)
    points(jitter(filtset$ar1.mean, 5),jitter(filtset$ar2.mean, 5), col=cols[i], pch=19, cex=cx)
    abline(v=filtset$truear1, h=filtset$truear2, lty = 2, col = "gray90")
    
  
  }
  dev.off()
  # 
  ## POISSON FIGURE ##
  setEPS()
  postscript(paste0("/Users/pni026/OneDrive - UiT Office 365/Results/p1vp2_p_v",round(v*100),".eps"),width = 6, height = 6)
  par(mar=c(5, 5, 4, 2) + 0.1, cex.lab=2)
  plot(1, type="n", ylim=c(-.9, .1), xlim=c(-1,1), 
       xlab = expression(phi["1"]), ylab = expression(phi["2"]), cex = cx, main="")
  # add background lines
  abline(v=filtset$truear1, h=filtset$truear2, lty = 2, col = "gray90")
  for(i in 8:5)
  {
    filtset <- filter(check1, method==i & inn.var == v)
    points(jitter(filtset$ar1.mean,5),jitter(filtset$ar2.mean,5), col=cols[i-4], pch=19, cex=cx)
  }
  dev.off()
}

# EPS #
#### image plot ####
require(fields)
colfunc <- colorRampPalette(c("red","yellow"))
brks <- seq(.4,1.03,.02)
colx <-colfunc(length(brks)-1)
par(mfrow=c(2,4))
for(phi in 1:2){
for(v in c(.04,.08,0.16,0.32)){
  filtset1 <- filter(ttt2, inn.var == v)
  cov1 <- aggregate(cbind(ar1cov,ar2cov,ar12cov,RMSE1,RMSE2)~method+truear1, data=filtset1,FUN=mean)
  print(min(cov1$ar1cov,cov1$ar2cov))
  cov2 <- as.matrix(xtabs(ar1cov~truear1+method,data=cov1))
  cov3 <- as.matrix(xtabs(ar2cov~truear1+method,data=cov1))
  
  if (phi==1){
    if(v!=0.32){
      image(seq(-1,1,.5),1:8,cov2, main = "", xlab = expression(phi[1]), ylab = "Method", breaks=brks,col=colx,
        sub = paste0("Var = ",v))
      abline(h=4.5, lty = 1, col = "gray10", cex = 5)}
    else    
      {image.plot(seq(-1,1,.5),1:8,cov2, main =  "", xlab = expression(phi[1]), ylab = "Method", 
                         breaks=brks,col=colx,sub = paste0("Var = ",v))}
  }
  else 
    {image(seq(-1,1,.5),1:8,cov3, main =   "", xlab = expression(phi[1]), ylab = "Method", 
              breaks=brks,col=colx,sub = paste0("Var = ",v))
    abline(h=4.5, lty = 1, col = "gray10", cex = 5)
    }
}
}

for(phi in 1:2){
  for(v in c(.04,.08,0.16,0.32)){
    filtset1 <- filter(ttt2, inn.var == v)
    cov1 <- aggregate(cbind(ar1cov,ar2cov,ar12cov,RMSE1,RMSE2)~method+truear2, data=filtset1,FUN=mean)
    cov2 <- as.matrix(xtabs(ar1cov~truear2+method,data=cov1))
    cov3 <- as.matrix(xtabs(ar2cov~truear2+method,data=cov1))
    if (phi==1){
      if(v!=0.32){
        image(c(-.8,-.5,-.2),1:8,cov2, xlab = expression(phi[2]), ylab = "Method", breaks=brks,col=colx)
        abline(h=4.5, lty = 1, col = "gray10", cex = 5)}
      else{    
        image.plot(c(-.8,-.5,-.2),1:8,cov3, xlab = expression(phi[2]), ylab = "Method", 
                   breaks=brks,col=colx)
      abline(h=4.5, lty = 1, col = "gray10", cex = 5)}
    }
    else {
      image(c(-.8,-.5,-.2),1:8,cov3, xlab = expression(phi[2]), ylab = "Method", 
                  breaks=brks,col=colx)
      abline(h=4.5, lty = 1, col = "gray10", cex = 5)}
    }
  }
### TABLE ###
head(check1)
check2 <- check1
check2$pois <- rep(c(0,0,0,0,1,1,1,1),nrow(check1)/8)
check2$method <- rep(c(1,2,3,4),nrow(check1)/4)
aggdata2 <- aggregate(cbind(ar1cov,RMSE1,ar2cov,RMSE2,ar12cov,RMSE12=RMSE1+RMSE2)~method+inn.var+pois, data=check2, FUN=mean)
aggdata2[,4:9] <- round((aggdata2[,4:9]),2)

(out1 <-as.matrix.data.frame(aggdata2))
write.csv2(out1,file="/Users/pni026/OneDrive - UiT Office 365/Results/aggable_cov_rmse_200sim_round2.csv")

legend("topright")

#### plot var ####

aggvar <- aggregate(cbind(ar1cov,ar2cov,ar12cov,RMSE1,RMSE2,RMSE12=RMSE1+RMSE2)~method+inn.var, data=check1, FUN=mean)
head(test.1)
par(mfrow=c(1,2))
par(mar=c(5, 4, 4, 2) + 0.1)
for(i in 1:4)
{  filtset1 <- filter(aggvar, method == i)
if(i == 1) plot(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), pch=19, ylim=c(0.3,1), xlab = "Variance", ylab = "Coverage")
else points(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), pch=19)
lines(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), lwd=1.5)
}
#legend("topright",c("True N","Multinomial HT", "VGAM","Raw counts", "Poisson"),col=cols, pch=19)

for(i in 1:4)
{  filtset1 <- filter(aggvar, method == i+4)
if(i == 1) plot(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), pch=19, ylim=c(0.3,1), xlab = "Variance", ylab = "Coverage")
else points(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), pch=19)
lines(filtset1$inn.var,filtset1$ar2cov, col = alpha(cols[i],.7), lwd=1.5)
}

# coverage plot
  
# loop through the methods
ylabs = c("Coverage Phi1", "Coverage Phi2","Coverage Phi1,2")
par(mfrow=c(1,3))
  for(c in c(0,4)){
    for(i in 1:4)
      {
      filtset1 <- filter(check1, method == i+c & inn.var == .08)
      filt2 <- aggregate(ar12cov~truear2,FUN=mean, data=filtset1)
      if(i == 1) plot(filt2[,1],filt2[,2], col = alpha(cols[i],.7), pch=19, ylim=c(0.3,1), ylab = "Coverage", xlab=expression(phi[1]))
      else points(filt2[,1],filt2[,2], col = alpha(cols[i],.7), pch=19)
      lines(filt2[,1],filt2[,2], col = alpha(cols[i],.7))
      }
    }
  
#### coolest plots ####
ylabs = c("Coverage Phi1", "Coverage Phi2","Coverage Phi1,2")
par(mfrow=c(1,1))
c=2


# Coverage phi1,2
for(v in c(.08,.16,.32,.64)){
  v=.04
  for(p in 0:1){
  postscript(paste0("MusData/plots/simulations/covphi12_5s_var",v*100,pois[p+1],".eps"),onefile=FALSE, family = "Times", width = 8, height = 6, horizontal = FALSE)
  par(mar=c(8, 5, 4, 2) + 0.1, cex.lab=1.5)
  plot(1, type="n", xlab="",ylab=expression(phi["1,2"] ~ Coverage), ylim=c(0.4, 1), xlim=c(1,15), xaxt="n",
       yaxt="n", frame.plot = FALSE)
  for(i in 1:4)
  {
    int <- 0
    for(chunk in 1:5){
    filtset1 <- filter(check1, method ==i+4 & inn.var == v)
    filt2 <- arrange(filtset1, truear1, truear2)[(1+int):(3+int),]
    
    points((1+int):(3+int),filt2[,11+c], col = alpha(cols[i],1), pch=19, cex=1.2)
    lines((1+int):(3+int),filt2[,11+c], col = alpha(cols[i],1), lwd=1.2)
    int = int+3
    print(int)
    }
  }
  ct=0
  for(j in 0:4){
  Axis(side=1, at=c(1+ct,2+ct,3+ct), labels = c(-.8,-.5,-.2))
  ct=ct+3
  
  }
  Axis(side=2)
  
  abline(v=c(3.5,6.5,9.5,12.5),lty=3)
  
  text(2, 1, labels = expression(phi[1]~"=" ~"-1.0"), cex=1.3)
  text(5, 1, labels = expression(phi[1]~"=" ~"-0.5"), cex=1.3)
  text(8, 1, labels = expression(phi[1]~"=" ~"0.0"), cex=1.3)
  text(11, 1, labels = expression(phi[1]~"=" ~"0.5"), cex=1.3)
  text(14, 1, labels = expression(phi[1]~"=" ~"1.0"), cex=1.3)
  
  mtext(expression(phi[2]), side=1, line=3.5, at=2, cex=1.5)
  mtext(expression(phi[2]), side=1, line=3.5, at=5, cex=1.5)
  mtext(expression(phi[2]), side=1, line=3.5, at=8, cex=1.5)
  mtext(expression(phi[2]), side=1, line=3.5, at=11, cex=1.5)
  mtext(expression(phi[2]), side=1, line=3.5, at=14, cex=1.5)
  
  
  mtext(paste0("Variance = ",v), side=1, line=6, at=8, cex=1.2)
  
  dev.off()
}
}

###
pois<-c("","p")
p=0
for(v in c(.08,.16,.32,.64)){
  for(p in 0:1){
    postscript(paste0("MusData/plots/simulations/covphi12_3s_var",v*100,pois[p+1],".eps"),onefile=FALSE, family = "Times", width = 8, height = 6, horizontal = FALSE)
    par(mar=c(8, 5, 4, 2) + 0.1, cex.lab=1.5)
    plot(1, type="n", xlab="", ylab=expression(phi["1,2"] ~ Coverage), ylim=c(0.4, 1), xlim=c(1,15), xaxt="n",
    yaxt="n", frame.plot = FALSE)
      for(i in 1:4)
      {
        
        int <- 0
        for(chunk in 1:3){
          filtset1 <- filter(check1, method ==(i+ps[p+1]) & inn.var == v)
          filt2 <- arrange(filtset1, truear2, truear1)[(1+int):(5+int),]
          
          points((1+int):(5+int),filt2[,11+c], col = alpha(cols[i],1), pch=19-p*4)
          lines((1+int):(5+int),filt2[,11+c], col = alpha(cols[i],1), lwd=2,lty=1+p)
          int = int+5
          print(int)
      }

        
      }
    ct=0
    for(j in 0:2){
      Axis(side=1, at=c(1+ct,2+ct,3+ct,4+ct,5+ct), labels = c(-1,-.5,0,.5,1))
      ct=ct+5
    }
    
    abline(v=c(5.5,10.5),lty=3)
    text(3, 1, labels = expression(phi[2]~"=" ~"-0.8"), cex=1.3)
    text(8, 1, labels = expression(phi[2]~"=" ~"-0.5"), cex=1.3)
    text(13, 1, labels = expression(phi[2]~"=" ~"-0.2"), cex=1.3)

    Axis(side=2)
    
    mtext(expression(phi[1]), side=1, line=3.5, at=3, cex=1.5)
    mtext(expression(phi[1]), side=1, line=3.5, at=8, cex=1.5)
    mtext(expression(phi[1]), side=1, line=3.5, at=13, cex=1.5)
    mtext(paste0("Variance = ",v), side=1, line=6, at=8, cex=1.2)
    dev.off()
  }
  }




##
for(v in c(.04,.08,.16,.32)){
  # pdf(paste0("/Users/pni026/OneDrive - UiT Office 365/Results/coverage_var",v,"_asphi1.pdf"),width = 10, height = 10)
  par(mfrow=c(2,2))
  for(p in 1:2){
    for(c in 1:2){
      if(c==1){
        plot(1, type="n", xlab=expression(phi[1]), ylab=substitute(Coverage ~ phi[1]), ylim=c(0, 1),  xlim=c(1,15), xaxt="n",
             sub=paste0(Plab[p],"Variance = ",v), frame.plot=F)
        }
      if(c==2){
        plot(1, type="n", xlab=expression(phi[1]), ylab=substitute(Coverage ~ phi[2]), ylim=c(0, 1),  xlim=c(1,15), xaxt="n",
             sub=paste0(Plab[p],"Variance = ",v), frame.plot=F)}
      for(i in 1:4)
      {
        
        int <- 0
        for(chunk in 1:5){
          filtset1 <- filter(check1, method ==(i+ps[p]) & inn.var == v)
          filt2 <- arrange(filtset1, truear1, truear2)[(1+int):(3+int),]
          
          points((1+int):(3+int),filt2[,9+c], col = alpha(cols[i],.7), pch=19)
          lines((1+int):(3+int),filt2[,9+c], col = alpha(cols[i],.7), lwd=2)
          int = int+3
          #print(int)
        }
      }
      Axis(side=1, at=seq(2,14,3), labels = seq(-1,1,.5))
    }}
  # dev.off()
}


ps <- c(0,4)
Plab = c("","(P) ")

for(v in c(.04,.08,.16,.32)){
  # pdf(paste0("/Users/pni026/OneDrive - UiT Office 365/Results/coverage_var",v,".pdf"),width = 10, height = 10)
  par(mfrow=c(2,2))
for(p in 1:2){
for(c in 1:2){
  if(c==1){
  plot(1, type="n", xlab=expression(phi[2]), ylab=substitute(Coverage ~ phi[1]), ylim=c(0, 1), xlim=c(1,15), xaxt="n",
       sub=paste0(Plab[p],"Variance = ",v))
    }
  if(c==2){
    plot(1, type="n", xlab=expression(phi[2]), ylab=substitute(Coverage ~ phi[2]), ylim=c(0, 1), xlim=c(1,15), xaxt="n",
         sub=paste0(Plab[p],"Variance = ",v))}
  
  for(i in 1:4)
  {
    
    int <- 0
    for(chunk in 1:3){
      filtset1 <- filter(check1, method ==(i+ps[p]) & inn.var == v)
      filt2 <- arrange(filtset1, truear2, truear1)[(1+int):(5+int),]
      
      points((1+int):(5+int),filt2[,9+c], col = alpha(cols[i],.7), pch=19)
      lines((1+int):(5+int),filt2[,9+c], col = alpha(cols[i],.7), lwd=2)
      int = int+5
      print(int)
    }
  }
  Axis(side=1, at=c(3,8,13), labels = c(-.8,-.5,-.2))
}}
  # dev.off()
  }
# plot coverage as a function of setting start with -1 phi1 and go through phi2 then -.5 ...
# get the same for RMSE
# heat maps for the RMSE
# plot covphi1 and covphi2 woth two types of lines for each of the 4 methods
# RMSE sqrt(sum(bias2^2+bias2^2))


# caption #
setEPS()
postscript("Results/caption.eps",width = 6, height = 1)
par(mar=c(1, 1, 1, 1))
plot(1, type="n",yaxt="n", xaxt="n", ylab="", xlab="", frame.plot=F)
legend("center", legend=c("Baseline"," CR-INLA"," CR-VGAM ","ObsCount"),col=cols, pch=19, horiz=TRUE, cex=1.1)
dev.off()
par(mar=c(5, 4, 4, 2) + 0.1)
