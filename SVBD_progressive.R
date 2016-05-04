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

##### Reload models to create initial values########
load(file = "GrowthAK")
load(file = "GrowthAKY")

Ai<-GrowthAK$BUGSoutput$mean$A
ki<-GrowthAK$BUGSoutput$mean$k
ki[444]<-0.0249
bdi<-GrowthAK$BUGSoutput$mean$bd
mA<-GrowthAKY$BUGSoutput$mean$mA
mk<-GrowthAKY$BUGSoutput$mean$mk


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

save(GrowthAKl,file="GrowthAKl")

plot(GrowthAK$BUGSoutput$mean$A,GrowthAKl$BUGSoutput$mean$A)
cor.test(GrowthAK$BUGSoutput$mean$A,GrowthAKl$BUGSoutput$mean$A)
#Hopefully GrowthAKl is closer to reality, but in any case, correlation of 78% with and without litters!!!

### put back Aobs
sink("models/GrowthAKlAobs")
cat("
    model {
    ######priors and constraints
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
        litter[i] ~ dcat(plG[i,])
      }    
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
dataGrowthAKlphi<-list(mothernb=mothernb,M=mass$Weight,
                    t=mass$RelativeJulian,IndCohort=IndCohort,
                    nind=length(unique(mass$id)),nobs=length(mass$Weight),
                    pl=pld,idmum=idmum,MaxLitter=MaxLitter,mother=mother,
                    whichind=whichind,maxbd=maxbd,mA=mA,mk=mk,A=Aobs)
initsGAKlphi <- function() list(A=initA,k=ki,pl=pli,meanLB=initmeanLB)

paramsGAKlphi <- c("A","k","litter","meanLB")
# MCMC settings
ni <- 13000 ; nt <- 20 ; nb <- 3000 ; nc <- 3
GrowthAKlAobs<-jags(dataGrowthAKlphi,initsGAKlphi,paramsGAKlphi,"models/GrowthAKlAobs",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(GrowthAKlAobs)

plot(GrowthAK$BUGSoutput$mean$A,GrowthAKlAobs$BUGSoutput$mean$A)
plot(GrowthAKlAobs$BUGSoutput$mean$A,GrowthAKl$BUGSoutput$mean$A,pch=16)

cor.test(GrowthAK$BUGSoutput$mean$A,GrowthAKl$BUGSoutput$mean$A)
cor.test(GrowthAKlAobs$BUGSoutput$mean$A,GrowthAKl$BUGSoutput$mean$A)
cor.test(GrowthAKlAobs$BUGSoutput$mean$A,GrowthAK$BUGSoutput$mean$A)
points(y=GrowthAKl$BUGSoutput$mean$A,x=Aobs,col="red")

cor.test(GrowthAK$BUGSoutput$mean$A,Aobs)

cov(GrowthAK$BUGSoutput$mean$A,Aobs, use = "complete.obs")/sqrt(var(GrowthAK$BUGSoutput$mean$A)*(var(Aobs, na.rm=T)-2.05^2))

cor.test(GrowthAKl$BUGSoutput$median$A,Aobs)


Aestimates <- data.frame(id= names(Aobs), A=GrowthAK$BUGSoutput$mean$A)
head(Aestimates)

testmerge <- merge(x=YearPheno, y=Aestimates,by.x = "ID", by.y="id",all.x = TRUE)
dim(testmerge)
dim(YearPheno)

testmerge$A[testmerge$Age=="A"] <- testmerge$Mass[testmerge$Age=="A"]

plot(testmerge$A,testmerge$Mass)

summary(glm(data=testmerge, Phi ~1+Sex+Age*A, family=binomial))
summary(glm(data=testmerge, Fitness ~1+Sex+Age+Mass))


write.table(x = Aestimates , file = "AtoMerge.txt",quote = FALSE, row.names = FALSE)
