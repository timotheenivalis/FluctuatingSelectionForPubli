#setwd("C:/Users/Timoth?e/Dropbox")
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
library(nlme)
library(R2jags)
library(reshape2)
mass<-read.table(file="ForFluctuatingSelectionRaw.txt",header=T,stringsAsFactors=F)
mass<-mass[which(mass$Age=="J"),]


###add julian date
  mass$Origin=as.POSIXct(paste(paste(rep("01/01/",dim(mass)[1])),mass$Calendar_year,sep=""),format="%d/%m/%Y",tz="GMT")
mass$Julian<-0
for (i in 1:dim(mass)[1])
{
  mass$Julian[i]<-round(julian(as.POSIXct(mass$Date[i],format="%d/%m/%Y"),origin=as.POSIXct(mass$Origin[i],format="%Y/%m/%d")),digit=0)
}
mass$RelativeJulian<-mass$Julian-min(mass$Julian)
plot(100*mass$Weight/mass$Body_Length,x=mass$RelativeJulian)

ind1<-mass[which(mass$ID_Individual==names(which(table(mass$ID_Individual)==7)[1])),]

mass06<-mass[which(mass$Calendar_year==2006),]


maxt<-max(table(mass06$ID_Individual))
ind06<-unique(mass06$ID_Individual)
Wind<-matrix(NA,nrow=length(ind06),ncol=maxt)
Tind<-matrix(NA,nrow=length(ind06),ncol=maxt)
for(i in 1:length(ind06))
{
  subind<-mass06[which(mass06$ID_Individual==ind06[i]),]
  for (d in 1:dim(subind)[1])
  {
    Wind[i,d]<-100*subind$Weight[d]/subind$Tail_Length[d]
    Tind[i,d]<-subind$RelativeJulian[d]
  }
}
Tind[which(is.na(Tind))]<-160

sink("Growth1")
cat("
    model {
    #priors and constraints
    A~dunif(10,40)
    k~dunif(0.01,0.99)
    tau<-pow(sd,-2)
    sd~dunif(0,10)
    for (ind in 1:nind)
    {
    b[ind]~dunif(0,60)
    }
    #likelihood
    for (ind in 1:nind)
    {
    for (i in 1:Temps)    
    {
    M[ind,i]~dnorm(mf[ind,i],tau)
    mf[ind,i]<-A*(1-exp(-k*(t[ind,i]-b[ind])))
    }
    }
    
    }
    ",fill = TRUE)
sink()

dataGrowth<-list(M=Wind,t=Tind,Temps=maxt,nind=length(ind06))
initsG <- function() list(A = runif(n=1,10,40), k=runif(n=1,0.01,0.99), b=runif(n=ind06,0,60))
paramsG <- c("A","k","sd")
ni <- 5000 ; nt <- 10 ; nb <- 1000 ; nc <- 3
Growth1<-jags(dataGrowth,initsG,paramsG,"Growth1", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(Growth1)

sink("Growth1-0")
cat("
    model {
    #priors and constraints
    I~dunif(10,100)
    B~dunif(0.01,0.99)
    B2~dunif(-1,1)
    tau<-pow(sd,-2)
    sd~dunif(0,10)
    for (ind in 1:nind)
    {
    b[ind]~dunif(0,60)
    }
    #likelihood
    for (ind in 1:nind)
    {
    for (i in 1:Temps)    
    {
    M[ind,i]~dnorm(mf[ind,i],tau)
    mf[ind,i]<-I+B*(t[ind,i]-b[ind])+B2*pow(t[ind,i]-b[ind], 2)
    }
    }
    
    }
    ",fill = TRUE)
sink()
dataGrowth<-list(M=Wind,t=Tind,Temps=maxt,nind=length(ind06))
initsG <- function() list(I = runif(n=1,10,100), B=runif(n=1,0.01,0.99),B2=0, b=runif(n=ind06,0,60))
paramsG <- c("I","B","B2","sd")
ni <- 5000 ; nt <- 10 ; nb <- 1000 ; nc <- 3
Growth1_0<-jags(dataGrowth,initsG,paramsG,"Growth1-0", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(Growth1_0)

curve(expr = 45.5 + 0.318*x - 0.003*(x^2), from = 0, to = 100)
curve(expr = 31.969868  + 0.479558*x - 0.002658*(x^2), from = 0, to = 100)

plot(Wind[,1], x=Tind[,1])
summary(lm(Wind[,1]~ Tind[,1]+I(Tind[,1]^2)))

min(Tind)
max(Tind)
