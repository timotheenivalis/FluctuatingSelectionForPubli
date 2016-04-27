#setwd(dir = "thesis/BirthDates/")
library(R2jags)
source(file = "functionsForBirthDates.R")
AllM<-read.table(file="ForFluctuatingSelectionRaw.txt",header=T,stringsAsFactors=F)
names(AllM)[1] <- "id"
#add julian date
AllM$Origin=as.POSIXct(paste(paste(rep("01/01/",nrow(AllM))),AllM$Calendar_year,sep=""),format="%d/%m/%Y",tz="GMT")
AllM$Julian<-0
for (i in 1:nrow(AllM))
{
  AllM$Julian[i]<-round(julian(as.POSIXct(AllM$Date[i],format="%d/%m/%Y"),origin=as.POSIXct(RawPheno$Origin[i],format="%Y/%m/%d")),digit=0)
}
AllM$RelativeJulian<-AllM$Julian-min(AllM$Julian)

AllMj<-AllM[which(AllM$Age=="J"),]

### looking for oldest juv of each female # ordering by mothers and by juv within mothers
k<-0.039 
A<-31.8

AllMj$A0<-sapply(X = AllMj$Weight,FUN = function(x){max(x,A,na.rm=TRUE)})
AllMj$bd0<-AllMj$RelativeJulian+(1/k)*log(1-AllMj$Weight/(AllMj$A0+1))

orphans<-AllMj$id[which(AllMj$Mother==0)]
for (i in 1:length(orphans))
{
  AllMj$Mother[which(AllMj$id==orphans[i])]<-paste("M",i,sep="")
}

AllMj$MY<-paste(AllMj$Mother,AllMj$Year,sep="")#unique Mother X year


AllMj$ord<-NA
AllMj$Elder<-FALSE
motherlist<-unique(AllMj$MY)
count<-0
for (i in 1:length(motherlist))
{
  offsp<-unique(AllMj$id[which(paste(AllMj$Mother,AllMj$Year,sep="")==motherlist[i])])
  if(length(offsp)!=0)
  {
    #get the mean bd0 and rank them according to it
    meanbd0<-tapply(X = AllMj$bd0[which(AllMj$id %in% offsp)],INDEX = AllMj$id[which(AllMj$id %in% offsp)],mean)
    ord<-order(meanbd0)
    for (j in 1:length(ord))
    {
      AllMj$ord[which(AllMj$id==names(meanbd0[ord[j]]))]<-count+ord[j]
    }
    AllMj$Elder[which(AllMj$id==names(meanbd0[ord[1]]))]<-TRUE
    count<-count+length(ord)
  }
}

AllMj<-AllMj[order(AllMj$ord),]#data frame ready

mass<-AllMj[which(!is.na(AllMj$Weight)),]

maxbd<-max(mass$RelativeJulian)-14

mums<-unique(mass$MY)
MumsPups<-unique(mass$id)
mothernb<-length(mums)

nind<-length(MumsPups)
mother<- sapply(MumsPups,FUN = function(x){which(mums==mass$MY[which(mass$id==x)[1]])})

MaxLitter<-5
pldlist<-F_initPl(MaxLitter = MaxLitter,mums = mums,mother = mother)

pld<-pldlist$pld
pli<-pldlist$pli
initmeanLB<-pldlist$initmeanLB
littermum<-pldlist$littermum

whichind<-vector(length=nrow(mass))
for (i in 1:length(whichind))
{
  whichind[i]<-which(MumsPups==mass$id[i])
}
whichind<-as.integer(as.factor(whichind))

idmum<-(c(1:mothernb)-1)*MaxLitter

IndCohort<-tapply(mass$Cohort,mass$id,function(x){mean(x)-2005})

Aobs<-sapply(X = MumsPups,FUN = function(x){mean(AllM$Weight[which(AllM$id==x & AllM$Age=="A")])})
Aobs[which(is.nan(Aobs))]<-NA
initA<-sapply(Aobs,FUN = function(x){ifelse(is.na(x),rnorm(n = 1,mean = 33,sd = 2),no=NA)})


##### Just estimating A and k per year

sink("models/GrowthAKY")
cat("
    model {
    ######priors and constraints
    #growth

    mA~dunif(25,45)  
    mk~dunif(0.01,0.08)


    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ######likelihood
    for (obs in 1:nobs)
      {
        M[obs]~dnorm(mf[obs],tau)
        mf[obs]<-mA*(1-exp(-mk*(date[obs])))
        date[obs]<-t[obs]-bd[whichind[obs]]
      }
    for (i in 1:nind)
      {
        bd[i]~dunif(-30,maxbd)
      }    
    }
    ",fill = TRUE)
sink()
dataGrowthAKY<-list(M=mass$Weight,
                  t=mass$RelativeJulian,IndCohort=IndCohort,ncohorts=max(IndCohort),
                  nind=length(unique(mass$id)),nobs=length(mass$Weight),
                  whichind=whichind,maxbd=maxbd)
initsGAKY <- function() list(bd=runif(n = length(unique(mass$id)),min = -20,max = maxbd),
                             mA=runif(1,25,45),mk=runif(1,0.01,0.08))

paramsGAKY <- c("mA","mk","bd")
# MCMC settings
ni <- 13000 ; nt <- 10 ; nb <- 3000 ; nc <- 3
GrowthAKY<-jags(dataGrowthAKY,initsGAKY,paramsGAKY,"models/GrowthAKY",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(GrowthAKY)
traceplot(GrowthAKY,varname="mk")

save(GrowthAKY,file = "GrowthAKY")#to be retrieved if need starting values for bd and mean of mA and mk

startingBD<-GrowthAKY$BUGSoutput$mean$bd
plot(startingBD)
mA<-GrowthAKY$BUGSoutput$mean$mA
mk<-GrowthAKY$BUGSoutput$mean$mk


#variable k, A, bd, litter, year by year
sink("models/GrowthAK")
cat("
    model {
    ######priors and constraints
    #growth
    for (i in 1:nind)
      {
        #A[i]~dnorm(mA,tauA)
        A[i]~dunif(mA-5,mA+15)
        k[i]~dnorm(mk,tauk)T(0.01,0.04)
      }
    tauA<-pow(sdA,-2)    
    sdA~dunif(0,10)
    tauk<-pow(sdk,-2)
    sdk~dunif(0,0.01)
  
    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ######likelihood
    for (obs in 1:nobs)
      {
        M[obs]~dnorm(mf[obs],tau)
        mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
        date[obs]<-t[obs]-bd[whichind[obs]]
      }
    for (i in 1:nind)
      {
        bd[i]~dunif(-30,maxbd)
      }    
    }
    ",fill = TRUE)
sink()

dataGrowthAK<-list(M=mass$Weight,
                    t=mass$RelativeJulian,IndCohort=IndCohort,
                    nind=length(unique(mass$id)),nobs=length(mass$Weight),
                    whichind=whichind,maxbd=maxbd,mA=mA,mk=mk)
initsGAK <- function() list(bd=startingBD,A=rnorm(n = length(unique(mass$id)),mean = 35,sd = 1),sdA=sd(initA,na.rm=T),k=rnorm(length(unique(mass$id)),mk,0.001))

paramsGAK <- c("A","k","sdA","bd","sdk")
# MCMC settings
ni <- 13000 ; nt <- 10 ; nb <- 3000 ; nc <- 3
GrowthAK<-jags(dataGrowthAK,initsGAK,paramsGAK,"models/GrowthAK",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(GrowthAK)
save(GrowthAK,file="GrowthAK")

traceplot(GrowthAK,varname="deviance")
traceplot(GrowthAK,varname="sdA")
traceplot(GrowthAK,varname="sdk")

plot(GrowthAK$BUGSoutput$mean$k)
plot(GrowthAK$BUGSoutput$mean$bd)
plot(GrowthAK$BUGSoutput$mean$A)
plot(GrowthAK$BUGSoutput$mean$bd,startingBD)
cor.test(GrowthAK$BUGSoutput$mean$bd,startingBD)

##### Reload models to create initial values########
load(file = "GrowthAK")
load(file = "GrowthAKY")

Ai<-GrowthAK$BUGSoutput$mean$A
ki<-GrowthAK$BUGSoutput$mean$k
ki[434]<-0.0249
bdi<-GrowthAK$BUGSoutput$mean$bd
mA<-GrowthAKY$BUGSoutput$mean$mA
mk<-GrowthAKY$BUGSoutput$mean$mk



plot(GrowthAK$BUGSoutput$mean$A,initA,ylim=c(27,42))
plot(GrowthAK$BUGSoutput$mean$A,Aobs,col="red",pch=16)
cor.test(GrowthAK$BUGSoutput$mean$A,Aobs,use = "complete.obs")

summary(glm(ifelse(is.na(initA),1,0)~GrowthAK$BUGSoutput$mean$A*GrowthAK$BUGSoutput$mean$bd+I(GrowthAK$BUGSoutput$mean$bd^2),family="binomial"))

summary(glm(ifelse(is.na(initA),1,0)~GrowthAK$BUGSoutput$mean$k,family="binomial"))
summary(glm(ifelse(is.na(initA),1,0)~GrowthAK$BUGSoutput$mean$A,family="binomial"))
summary(glm(ifelse(is.na(initA),1,0)~beforeWinter,family="binomial"))
summary(glm(ifelse(is.na(initA),1,0)~IndCohort*GrowthAK$BUGSoutput$mean$A,family="binomial"))
summary(glm(ifelse(is.na(initA),1,0)~GrowthAK$BUGSoutput$mean$k+GrowthAK$BUGSoutput$mean$A + beforeWinter,family="binomial"))
summary(glm(ifelse(is.na(initA),1,0)~GrowthAK$BUGSoutput$mean$A*beforeWinter,family="binomial"))

beforeWinter<-endSeason[IndCohort]-GrowthAK$BUGSoutput$mean$bd
library(ggplot2)
outgg<-data.frame(phi=ifelse(is.na(initA),1,0),bw=beforeWinter,k=GrowthAK$BUGSoutput$mean$k,A=GrowthAK$BUGSoutput$mean$A)

c <- ggplot(outgg, aes(bw,phi))
c + stat_smooth()

c <- ggplot(outgg, aes(k,phi))
c + stat_smooth()

c <- ggplot(outgg, aes(A,phi))
c + stat_smooth()

## now put back the litters!
sink("models/GrowthAKl")
cat("
    model {
    ######priors and constraints
    #growth
    for (i in 1:nind)
      {
        A[i]~dunif(mA-5,mA+15)
        k[i]~dnorm(mk,tauk)T(0.01,0.025)
      }
    tauk<-pow(sdk,-2)
    sdk~dunif(0,0.001)
    
    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ## Birth dates
    for (mum in 1:mothernb)
      {
          #first litter
          meanLB[mum,1]~dunif(-30,maxbd)       
          for (i in 2:MaxLitter)
            {
                meanLB[mum,i]~dunif(gbegin[mum,i],gend[mum,i])
                gbegin[mum,i]<-meanLB[mum,i-1]+20
                gend[mum,i]<-meanLB[mum,i-1]+120
            }
      }
    for (i in 1:nind)
      {
        for(j in 1:MaxLitter)
            {
                pl[i,j]~dgamma(1,1)
                plG[i,j]<-pl[i,j]/sum(pl[i,])
            }
      }
    ######likelihood
    for (obs in 1:nobs)
      {
          M[obs]~dnorm(mf[obs],tau)
          mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
          date[obs]<-t[obs]-meanLB[mother[whichind[obs]],litter[whichind[obs]]]
      }
    for (i in 1:nind)
      {
          litter[i] ~ dcat(plG[i,])
          #correctionlitter[i]<-(mother[i]-1)*MaxLitter
          #littermum[i]<-litter[i]-correctionlitter[i]
      }    
    }
    ",fill = TRUE)
sink()
sparceinitPl<-F_initPlsparse(MaxLitter = MaxLitter,mums = mums,mother = mother)
pld<-sparceinitPl$pld
pli<-sparceinitPl$pli
initmeanLB<-sparceinitPl$initmeanLB
dataGrowthAKl<-list(mothernb=mothernb,M=mass$Weight,
                   t=mass$RelativeJulian,IndCohort=IndCohort,
                   nind=length(unique(mass$id)),nobs=length(mass$Weight),
                   pl=pld,idmum=idmum,MaxLitter=MaxLitter,mother=mother,
                   whichind=whichind,maxbd=maxbd,mA=mA,mk=mk)
initsGAKl <- function() list(bd=bdi,A=Ai,k=ki,pl=pli,meanLB=initmeanLB)

paramsGAKl <- c("A","k","bd","litter","meanLB")
# MCMC settings
ni <- 13000 ; nt <- 10 ; nb <- 3000 ; nc <- 3
GrowthAKl<-jags(dataGrowthAKl,initsGAKl,paramsGAKl,"models/GrowthAKl",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(GrowthAKl)
traceplot(GrowthAKl,varname="litter")
plot(GrowthAKl$BUGSoutput$mean$A)
plot(GrowthAKl$BUGSoutput$sd$A)
plot(GrowthAKl$BUGSoutput$mean$litter)


sd(runif(n = 10000,min = 0,max = 20))

hist(Aobs)

### put back phi and Aobs
sink("models/GrowthAKlphi")
cat("
    model {
    ######priors and constraints
    #survival
    meanmu~dnorm(0,0.001)
    mean.phi<-1/(1+exp(-meanmu))
    BetaA~dnorm(0,0.001)
    BetaD~dnorm(0,0.001)
    BetaAD~dnorm(0,0.001)
    BetaS~dnorm(0,0.001)

    #growth
    for (i in 1:nind)
      {
          A[i]~dunif(mA-5,mA+16)
          k[i]~dnorm(mk,tauk)T(0.01,0.025)
      }
    tauk<-pow(sdk,-2)
    sdk~dunif(0,0.001)
    
    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ## Birth dates
    for (mum in 1:mothernb)
      {
        #first litter
        meanLB[mum,1]~dunif(-30,maxbd)       
        for (i in 2:MaxLitter)
          {
            meanLB[mum,i]~dunif(gbegin[mum,i],gend[mum,i])
            gbegin[mum,i]<-meanLB[mum,i-1]+20
            gend[mum,i]<-meanLB[mum,i-1]+120
          }
      }
    for (i in 1:nind)
      {
        for(j in 1:MaxLitter)
          {
            pl[i,j]~dgamma(1,1)
            plG[i,j]<-pl[i,j]/sum(pl[i,])
          }
      }
    ######likelihood
    for (obs in 1:nobs)
      {
        M[obs]~dnorm(mf[obs],tau)
        mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
        date[obs]<-t[obs]-meanLB[mother[whichind[obs]],litter[whichind[obs]]]
      }
    for (i in 1:nind)
      {
        DeltaWinter[i]<-endSeason[IndCohort[i]]-meanLB[mother[i],litter[i]]
        litter[i] ~ dcat(plG[i,])
        phi[i]~dbin(p[i],1)
        logit(p[i])<-meanmu+BetaA*A[i]+BetaD*DeltaWinter[i]+BetaAD*A[i]*DeltaWinter[i]+BetaS*sex[i]
        pa[i]<-phi2[i]*A[i]
      }    
    covAPhi<-mean(pa)-mean(A)*mean(phi2)    
    }
    ",fill = TRUE)
sink()

Aobs<-sapply(Aobs,FUN=function(x){ifelse(x<31,31,x)})
initA<-Ai
initA[which(!is.na(Aobs))]<-NA

sparceinitPl<-F_initPlsparse(MaxLitter = MaxLitter,mums = mums,mother = mother)
pld<-sparceinitPl$pld
pli<-sparceinitPl$pli
initmeanLB<-sparceinitPl$initmeanLB
dataGrowthAKlphi<-list(mothernb=mothernb,M=mass$Weight,endSeason=endSeason,
                    t=mass$RelativeJulian,IndCohort=IndCohort,
                    nind=length(unique(mass$id)),nobs=length(mass$Weight),
                    pl=pld,idmum=idmum,MaxLitter=MaxLitter,mother=mother,
                    whichind=whichind,maxbd=maxbd,mA=mA,mk=mk,phi=phis,phi2=phi2,A=Aobs)
initsGAKlphi <- function() list(A=initA,k=ki,pl=pli,meanLB=initmeanLB,BetaA=0,BetaD=0,BetaAD=0,BetaS=0,meanmu=0)

paramsGAKlphi <- c("A","k","litter","meanLB","BetaA","BetaD","BetaAD","meanmu","covAPhi","BetaS")
# MCMC settings
ni <- 13000 ; nt <- 20 ; nb <- 3000 ; nc <- 3
GrowthAKlphi<-jags(dataGrowthAKlphi,initsGAKlphi,paramsGAKlphi,"models/GrowthAKlphi",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(GrowthAKlphi)
traceplot(GrowthAKlphi,varname="meanmu")
GrowthAKlphi.mcmc<-as.mcmc(GrowthAKlphi)
SGrowthAKlphi.mcmc<-summary(GrowthAKlphi.mcmc)
SGrowthAKlphi.mcmc$statistics["BetaD",];SGrowthAKlphi.mcmc$quantiles["BetaD",]
SGrowthAKlphi.mcmc$statistics["BetaA",];SGrowthAKlphi.mcmc$quantiles["BetaA",]
SGrowthAKlphi.mcmc$statistics["BetaAD",];SGrowthAKlphi.mcmc$quantiles["BetaAD",]
SGrowthAKlphi.mcmc$statistics["meanmu",];SGrowthAKlphi.mcmc$quantiles["meanmu",]


xA<-seq(from = 30,to = 51,by = 1)
yPhi0<-GrowthAKlphi$BUGSoutput$mean$meanmu+GrowthAKlphi$BUGSoutput$mean$BetaA*xA+GrowthAKlphi$BUGSoutput$mean$BetaD*0
yPhi150<-GrowthAKlphi$BUGSoutput$mean$meanmu+GrowthAKlphi$BUGSoutput$mean$BetaA*xA+GrowthAKlphi$BUGSoutput$mean$BetaD*150+GrowthAKlphi$BUGSoutput$mean$BetaAD*150*xA


plot(x=xA,y=1/(1+exp(-yPhi0)))
points(x=xA,y=1/(1+exp(-yPhi150)))


#### Again with interaction A:D ####
#### Also put sex in survival! 
#### And the direct covariance between A and Phi!
sink("models/GrowthAKlphi2")
cat("
    model {
    ######priors and constraints
    #survival
    meanmu~dnorm(0,0.001)
    mean.phi<-1/(1+exp(-meanmu))
    BetaA~dnorm(0,0.001)
    BetaD~dnorm(0,0.001)
    BetaAD~dnorm(0,0.001)
    BetaS~dnorm(0,0.001)
    
    tauphi2~dunif(0,10)
    dBetaA~dnorm(0,0.001)
    mu~dunif(0,1)
    #growth
    for (i in 1:nind)
{
    A[i]~dunif(mA-5,mA+16)
    k[i]~dnorm(mk,tauk)T(0.01,0.025)
}
    tauk<-pow(sdk,-2)
    sdk~dunif(0,0.001)
    
    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ## Birth dates
    for (mum in 1:mothernb)
{
    #first litter
    meanLB[mum,1]~dunif(-30,maxbd)       
    for (i in 2:MaxLitter)
{
    meanLB[mum,i]~dunif(gbegin[mum,i],gend[mum,i])
    gbegin[mum,i]<-meanLB[mum,i-1]+20
    gend[mum,i]<-meanLB[mum,i-1]+120
}
}
    for (i in 1:nind)
{
    for(j in 1:MaxLitter)
{
    pl[i,j]~dgamma(1,1)
    plG[i,j]<-pl[i,j]/sum(pl[i,])
}
}
    ######likelihood
    for (obs in 1:nobs)
{
    M[obs]~dnorm(mf[obs],tau)
    mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
    date[obs]<-t[obs]-meanLB[mother[whichind[obs]],litter[whichind[obs]]]
}
    for (i in 1:nind)
{
    DeltaWinter[i]<-endSeason[IndCohort[i]]-meanLB[mother[i],litter[i]]
    litter[i] ~ dcat(plG[i,])
    phi[i]~dbin(p[i],1)
    logit(p[i])<-meanmu+BetaA*A[i]+BetaD*DeltaWinter[i]+BetaAD*A[i]*DeltaWinter[i]+BetaS*sex[i]
    pa[i]<-phi[i]*A[i]
}    
    covAPhi<-mean(pa)-mean(A)*mean(phi)
    }
    ",fill = TRUE)
sink()

Aobs<-sapply(Aobs,FUN=function(x){ifelse(x<31,31,x)})
initA<-Ai
initA[which(!is.na(Aobs))]<-NA

sex<-tapply(X = mass$Sex,INDEX = mass$id,FUN = function(x){mean(x=="Male")})
indPhi<-tapply(X = mass$Phi,INDEX = mass$id,FUN = mean)
CohortPhi<-tapply(indPhi,IndCohort,mean)
#CohortPhiM<-tapply(indPhi[which(sex==1)],IndCohort[which(sex==1)],mean)
#CohortPhiF<-tapply(indPhi[which(sex==0)],IndCohort[which(sex==0)],mean)
phi2<-phis
#phi2[which(sex==0)]<-indPhi[which(sex==0)]/CohortPhiM[IndCohort[which(sex==0)]]
#phi2[which(sex==1)]<-indPhi[which(sex==1)]/CohortPhiF[IndCohort[which(sex==1)]]
phi2<-indPhi/CohortPhi[IndCohort]
plot(phi2)

sparceinitPl<-F_initPlsparse(MaxLitter = MaxLitter,mums = mums,mother = mother)
pld<-sparceinitPl$pld
pli<-sparceinitPl$pli
initmeanLB<-sparceinitPl$initmeanLB
dataGrowthAKlphi2<-list(mothernb=mothernb,M=mass$Weight,endSeason=endSeason,
                       t=mass$RelativeJulian,IndCohort=IndCohort,sex=sex,
                       nind=length(unique(mass$id)),nobs=length(mass$Weight),
                       pl=pld,MaxLitter=MaxLitter,mother=mother,
                       whichind=whichind,maxbd=maxbd,mA=mA,mk=mk,phi=phis,A=Aobs)
initsGAKlphi2 <- function() list(A=initA,k=ki,pl=pli,meanLB=initmeanLB,BetaA=0,BetaD=0,dBetaA=0,BetaS=0,meanmu=0,mu=0.5,tauphi2=0.1)

paramsGAKlphi2 <- c("A","k","litter","meanLB","BetaA","BetaD","BetaAD","meanmu","covAPhi","BetaS")
# MCMC settings
ni <- 230000 ; nt <- 200 ; nb <- 30000 ; nc <- 3
GrowthAKlphi2<-jags(dataGrowthAKlphi2,initsGAKlphi2,paramsGAKlphi2,"models/GrowthAKlphi2",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


print(GrowthAKlphi2)
traceplot(GrowthAKlphi2,varname="deviance")
traceplot(GrowthAKlphi2,varname="covAPhi")
traceplot(GrowthAKlphi2,varname="BetaD")
traceplot(GrowthAKlphi2,varname="BetaA")
traceplot(GrowthAKlphi2,varname="BetaAD")
traceplot(GrowthAKlphi2,varname="BetaS")

GrowthAKlphi2.mcmc<-as.mcmc(GrowthAKlphi2)
SGrowthAKlphi2.mcmc<-summary(GrowthAKlphi2.mcmc)
SGrowthAKlphi2.mcmc$statistics["covAPhi",]

SGrowthAKlphi.mcmc$statistics["meanmu",]
SGrowthAKlphi.mcmc$statistics["BetaA",]
GrowthAKlphi2$BUGSoutput$mean$BetaA
mean(as.mcmc((GrowthAKlphi2$BUGSoutput$sims.array[,1,"BetaA"]*var(GrowthAKlphi2$BUGSoutput$mean$A))))
HPDinterval(as.mcmc((GrowthAKlphi2$BUGSoutput$sims.array[,1,"BetaA"]*var(GrowthAKlphi2$BUGSoutput$mean$A))))
var(GrowthAKlphi2$BUGSoutput$mean$A)
cov(GrowthAKlphi2$BUGSoutput$mean$A,phi2)
cov(GrowthAKlphi2$BUGSoutput$mean$A[which(sex==1)],phi2[which(sex==1)])
cov(GrowthAKlphi2$BUGSoutput$mean$A[which(sex==0)],phi2[which(sex==0)])
cov(GrowthAKlphi2$BUGSoutput$mean$A[which(sex==1)],phis[which(sex==1)])
cov(GrowthAKlphi2$BUGSoutput$mean$A[which(sex==0)],phis[which(sex==0)])


xm<-25:50
ym<-GrowthAKlphi2$BUGSoutput$mean$meanmu+mean(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1,"BetaA"])*xm
ymH<-GrowthAKlphi2$BUGSoutput$mean$meanmu+(mean(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1,"BetaA"])+sd(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1,"BetaA"]))*xm
ymL<-GrowthAKlphi2$BUGSoutput$mean$meanmu+(mean(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1,"BetaA"])-sd(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1,"BetaA"]))*xm

predphi<-1/(1+exp(-ym))
predphiH<-1/(1+exp(-ymH))
predphiL<-1/(1+exp(-ymL))
plot(x=xm,y=predphi,type="l",lwd=3,xlab="Asymptotic mass",ylab="Predicted survival",ylim=c(0,1))
lines(x = xm,predphiH)
lines(x = xm,predphiL)

summary(glm(phis~GrowthAKlphi2$BUGSoutput$mean$A*sex))

plot(y=tapply(indPhi,IndCohort,mean),tapply(GrowthAKlphi2$BUGSoutput$mean$A,IndCohort,mean))
summary(lm(tapply(indPhi,IndCohort,mean)~tapply(GrowthAKlphi2$BUGSoutput$mean$A,IndCohort,mean)))

230000
autocorr(as.mcmc(GrowthAKlphi2$BUGSoutput$sims.array[1:1000,1:3,"BetaA"]),lags = c(1,10,20,30,40,50))

#@ run for longer @####
GrowthAKlphi2b<-jags(dataGrowthAKlphi2,initsGAKlphi2,paramsGAKlphi2,"models/GrowthAKlphi2",
                    n.chains = 3, n.thin = 6000, n.iter = 6300000, n.burnin = 300000, working.directory = getwd())
traceplot(GrowthAKlphi2b,varname="deviance")
traceplot(GrowthAKlphi2b,varname="covAPhi")
traceplot(GrowthAKlphi2b,varname="BetaD")
traceplot(GrowthAKlphi2b,varname="BetaA")
traceplot(GrowthAKlphi2b,varname="BetaAD")
traceplot(GrowthAKlphi2b,varname="BetaS")

mean(as.mcmc(as.vector(GrowthAKlphi2b$BUGSoutput$sims.array[,1:3,"BetaA"]))>0)
mean(as.mcmc(as.vector(GrowthAKlphi2b$BUGSoutput$sims.array[,1:3,"BetaAD"]))<0)

mD<-50#mean(endSeason)-mean(GrowthAKlphi2b$BUGSoutput$mean$meanLB[,1])
xm<-25:50
ym<-GrowthAKlphi2b$BUGSoutput$mean$meanmu+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaA"])*xm+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaD"])*mD+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaAD"])*mD*xm
ymH<-GrowthAKlphi2b$BUGSoutput$mean$meanmu+
  (mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaA"])+sd(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaA"]))*xm+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaD"])*mD+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaAD"])*mD*xm
ymL<-GrowthAKlphi2b$BUGSoutput$mean$meanmu+
  (mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaA"])-sd(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaA"]))*xm+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaD"])*mD+
  mean(GrowthAKlphi2b$BUGSoutput$sims.array[1:1000,1,"BetaAD"])*mD*xm

predphi<-1/(1+exp(-ym))
predphiH<-1/(1+exp(-ymH))
predphiL<-1/(1+exp(-ymL))
plot(x=xm,y=predphi,type="l",lwd=3,xlab="Asymptotic mass",ylab="Predicted survival",ylim=c(0,1))
lines(x = xm,predphiH)
lines(x = xm,predphiL)
#@ clean predictions @ ####
fullArray<-rbind(GrowthAKlphi2b$BUGSoutput$sims.array[,1,],GrowthAKlphi2b$BUGSoutput$sims.array[,2,],GrowthAKlphi2b$BUGSoutput$sims.array[,3,])
vcovbeta<-var(fullArray[,c("meanmu","BetaA","BetaD","BetaS","BetaAD")])
fixmeans<-matrix(data=as.numeric(GrowthAKlphi2b$BUGSoutput$mean[c("meanmu","BetaA","BetaD","BetaS","BetaAD")]))
newdat <- expand.grid(A=30:52,D=0:183,S=0:1)

mm <- as.matrix(data.frame(intercept=1,newdat,AD=newdat$A*newdat$D))#design matrix
newdat$logitphi <-  mm %*% fixmeans
pvar1 <- diag(mm %*% tcrossprod(vcovbeta,mm))
newdat <- data.frame(
  newdat
  , plo = newdat$logitphi-1.96*sqrt(pvar1)
  , phi = newdat$logitphi+1.96*sqrt(pvar1)
)
newdatRef<-newdat[which(newdat$S==0 & newdat$D==120),]
plot(1/(1+exp(-newdatRef$logitphi)),x=newdatRef$A,ylim=c(0,1),type="l",lwd=3,xlab="Asymptotic mass",ylab="Survival probability",las=1,col="blue")
lines(1/(1+exp(-newdatRef$plo)),x=newdatRef$A,col="blue")
lines(1/(1+exp(-newdatRef$phi)),x=newdatRef$A,col="blue")

newdatRef<-newdat[which(newdat$S==0 & newdat$D==68),]
lines(1/(1+exp(-newdatRef$logitphi)),x=newdatRef$A,ylim=c(0,1),lwd=3,col="red")
lines(1/(1+exp(-newdatRef$plo)),x=newdatRef$A,col="red")
lines(1/(1+exp(-newdatRef$phi)),x=newdatRef$A,col="red")

mean(1/(1+exp(-newdat[which(newdat$S==0 & newdat$D==80 & newdat$A<52),"logitphi"]))-1/(1+exp(-newdat[which(newdat$S==0 & newdat$D==80 & newdat$A>30),"logitphi"])))
lm(1/(1+exp(-newdat[which(newdat$S==0 & newdat$D==80 & newdat$A),"logitphi"]))~I(1:23))
lm(1/(1+exp(-newdat[which(newdat$S==0 & newdat$D==80 & newdat$A),"logitphi"]))[6:17]~I(6:17))

newdatRef<-newdat[which(newdat$S==0 & newdat$D==171),]
lines(1/(1+exp(-newdatRef$logitphi)),x=newdatRef$A,ylim=c(0,1),lwd=3,col="gray")
lines(1/(1+exp(-newdatRef$plo)),x=newdatRef$A,col="gray")
lines(1/(1+exp(-newdatRef$phi)),x=newdatRef$A,col="gray")

newdatRef<-newdat[which(newdat$S==0 & newdat$D==6),]
lines(1/(1+exp(-newdatRef$logitphi)),x=newdatRef$A,ylim=c(0,1),lwd=3,col="gray")
lines(1/(1+exp(-newdatRef$plo)),x=newdatRef$A,col="gray")
lines(1/(1+exp(-newdatRef$phi)),x=newdatRef$A,col="gray")


points(x = GrowthAKlphi2b$BUGSoutput$mean$A,y=phis)
       #y=rep(-0.01,length(GrowthAKlphi2b$BUGSoutput$mean$A)),pch=1)

c(325,312,276,286,292,261,288,282)-(tapply(X = IndBD$bd,IndBD$cohort,mean)+min(AllMj$Julian))
c(325,312,276,286,292,261,288,282)-(tapply(X = IndBD$bd,IndBD$cohort,min)+min(AllMj$Julian))
c(325,312,276,286,292,261,288,282)-(tapply(X = IndBD$bd,IndBD$cohort,max)+min(AllMj$Julian))

meanDyears<-c(325,312,276,286,292,261,288,282)-(tapply(X = IndBD$bd,IndBD$cohort,mean)+min(AllMj$Julian))
mean(meanDyears[-c(1,2)])

lastfree<-c(325,312,276,286,292,261,288,282)
firstfree<-c(151,148,131,124,171,146,163,153)

litterBD<-GrowthAKlphi2b$BUGSoutput$mean$meanLB
litternb<-GrowthAKlphi2b$BUGSoutput$median$litter
pbd<-1:length(IndCohort)
for (i in 1:length(IndCohort))
  {
    pbd[i]<-litterBD[mother[i],litternb[i]]
  }
IndBD<-data.frame(id=names(IndCohort),cohort=IndCohort+2005,bd=pbd)
head(IndBD)
points(x=IndBD$cohort,y=IndBD$bd+min(AllMj$Julian))

rawsnow<-read.table("rawsnow.txt",header=T)
rawsnow$hto000j0[which(is.na(rawsnow$hto000j0))]<-0
rawsnow<-rawsnow[which(!duplicated(paste(rawsnow$year,rawsnow$julian))),]

plot(c(lastfree,firstfree),x=rep(2006:2013,2),ylim=c(100,350),pch=16,col="white")

plot(x = 2006:2013,y=rep(0,8),ylim=c(90,350),pch=16,col="white")
for (y in 2006:2013)
{
  thisyear<-rawsnow[which(rawsnow$year==y ),]#& rawsnow$julian>firstfree[y-2005] & rawsnow$julian<lastfree[y-2005]),]
  points(thisyear$julian,x=rep(y,nrow(thisyear)),
         col=rgb(red = 0,green = 0,blue = 1,alpha =(1-(max(rawsnow$hto000j0,na.rm=T)-thisyear$hto000j0)/max(rawsnow$hto000j0,na.rm=T))^(1/5)),pch=16)
  #lines(x =c(y,y),y = c(0,firstfree[y-2005]),col=rgb(0,0,1,0.8),lwd=3)
  #lines(x =c(y,y),y = c(lastfree[y-2005],350),col=rgb(0,0,1,0.8),lwd=3)
}

af<-c(114,128,126,128,114,109,130,113,114,99,130,123,115,92,121,114,110)
al<-c(318,276,308,325,308,277,324,324,325,312,276,286,292,261,288,282,294)
plot(y=c(af,al),rep(c(1998:2014),2))
plot(al-af)
t.test((al-af)[3:10],y = (al-af)[10:17],var.equal = F)


####@ No interaction @####
sink("models/GrowthAKlphi3")
cat("
    model {
    ######priors and constraints
    #survival
    meanmu~dnorm(0,0.001)
    mean.phi<-1/(1+exp(-meanmu))
    BetaA~dnorm(0,0.001)
    BetaD~dnorm(0,0.001)
    BetaS~dnorm(0,0.001)
    
    tauphi2~dunif(0,10)
    dBetaA~dnorm(0,0.001)
    mu~dunif(0,1)
    #growth
    for (i in 1:nind)
{
    A[i]~dunif(mA-5,mA+16)
    k[i]~dnorm(mk,tauk)T(0.01,0.025)
}
    tauk<-pow(sdk,-2)
    sdk~dunif(0,0.001)
    
    tau<-pow(sdMass,-2)
    sdMass<-2.05
    
    ## Birth dates
    for (mum in 1:mothernb)
{
    #first litter
    meanLB[mum,1]~dunif(-30,maxbd)       
    for (i in 2:MaxLitter)
{
    meanLB[mum,i]~dunif(gbegin[mum,i],gend[mum,i])
    gbegin[mum,i]<-meanLB[mum,i-1]+20
    gend[mum,i]<-meanLB[mum,i-1]+120
}
}
    for (i in 1:nind)
{
    for(j in 1:MaxLitter)
{
    pl[i,j]~dgamma(1,1)
    plG[i,j]<-pl[i,j]/sum(pl[i,])
}
}
    ######likelihood
    for (obs in 1:nobs)
{
    M[obs]~dnorm(mf[obs],tau)
    mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
    date[obs]<-t[obs]-meanLB[mother[whichind[obs]],litter[whichind[obs]]]
}
    for (i in 1:nind)
{
    DeltaWinter[i]<-endSeason[IndCohort[i]]-meanLB[mother[i],litter[i]]
    litter[i] ~ dcat(plG[i,])
    phi[i]~dbin(p[i],1)
    logit(p[i])<-meanmu+BetaA*A[i]+BetaD*DeltaWinter[i]+BetaS*sex[i]
}    
    }
    ",fill = TRUE)
sink()

GrowthAKlphi3<-jags(dataGrowthAKlphi2,initsGAKlphi2,paramsGAKlphi2,"models/GrowthAKlphi3",
                     n.chains = 3, n.thin = 60, n.iter = 63000, n.burnin = 3000, working.directory = getwd())
traceplot(GrowthAKlphi3,varname="deviance")
traceplot(GrowthAKlphi3,varname="BetaD")
traceplot(GrowthAKlphi3,varname="BetaA")
traceplot(GrowthAKlphi3,varname="BetaS")

GrowthAKlphi3$BUGSoutput$mean$BetaA

plot(as.mcmc(as.vector(GrowthAKlphi3$BUGSoutput$sims.array[,1:3,"BetaA"])))
mean(as.mcmc(as.vector(GrowthAKlphi3$BUGSoutput$sims.array[,1:3,"BetaA"])))
HPDinterval(as.mcmc(as.vector(GrowthAKlphi3$BUGSoutput$sims.array[,1:3,"BetaA"])))
mean(as.mcmc(as.vector(GrowthAKlphi3$BUGSoutput$sims.array[,1:3,"BetaA"]))>0)


####@ Heritability of asymptotic size ####

litterBD<-GrowthAKlphi2b$BUGSoutput$mean$meanLB
litternb<-GrowthAKlphi2b$BUGSoutput$median$litter
pbd<-1:length(IndCohort)
for (i in 1:length(IndCohort))
{
  pbd[i]<-litterBD[mother[i],litternb[i]]
}

IndBD<-data.frame(id=names(IndCohort),cohort=IndCohort+2005,bd=pbd,A=GrowthAKlphi2b$BUGSoutput$mean$A)
IndBD2<-merge(x = IndBD,y=AllMj[!duplicated(AllMj$id),c("id","Sex","Mother","Father")],by = "id",all.x = T,all.y=F)
IndBD2$MeanAoff<-NA
IndBD2$SexRatio<-NA
for (i in 1:nrow(IndBD2))
{
  IndBD2$MeanAoff[i]<-mean(IndBD2$A[which(IndBD2$Mother==IndBD2$id[i] | IndBD2$Father==IndBD2$id[i])])
  IndBD2$SexRatio[i]<-mean(IndBD2$Sex[which(IndBD2$Mother==IndBD2$id[i] | IndBD2$Father==IndBD2$id[i])]=="Male")
}
plot(IndBD2$A,IndBD2$MeanAoff)
plot(IndBD2$A~IndBD2$bd)
2*lm(IndBD2$MeanAoff~IndBD2$A+as.factor(IndBD2$cohort))$coefficient[2]

summary(lm(MeanAoff~A+Sex+SexRatio+as.factor(cohort),data=IndBD2))

summary(lm(IndBD2$A~IndBD2$Sex+as.factor(IndBD2$cohort)))



arrayA3chains<-GrowthAKlphi2b$BUGSoutput$sims.array[,,which(substr(names(GrowthAKlphi2b$BUGSoutput$sims.array[1,1,]),1,1)=="A")]
arrayA<-rbind(arrayA3chains[,1,],arrayA3chains[,2,],arrayA3chains[,3,])
h2<-1:nrow(arrayA)
for (j in 1:nrow(arrayA))
{
  IndBD2$A<-arrayA[j,]
  for (i in 1:nrow(IndBD2))
  {
    IndBD2$MeanAoff[i]<-mean(IndBD2$A[which(IndBD2$Mother==IndBD2$id[i] | IndBD2$Father==IndBD2$id[i])])
  }
  h2[j]<-2*lm(IndBD2$MeanAoff~IndBD2$A+as.factor(IndBD2$cohort))$coefficient[2]
}
plot(as.mcmc(h2))
mean(h2)

trendA<-1:nrow(arrayA)
for (j in 1:nrow(arrayA))
{
  IndBD2$A<-arrayA[j,]
  trendA[j]<-lm(IndBD2$A~IndBD2$cohort)$coefficient[2]
}
plot(trendA)
