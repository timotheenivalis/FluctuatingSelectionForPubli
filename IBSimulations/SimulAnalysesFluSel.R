setwd("/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/IBSimulations/")

#### PRICE EQUATION #### 
evolsel <- read.table("CheckPriceEq/EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta) # weird pattern due to shrinkage of G variation? 
plot(evolsel$DeltaBV, x=evolsel$Var)

#### BREEDER'S EQUATION WITH H²=1####
evolsel <- read.table("CheckBreederEd_H2_1//EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(corv, ylim=c(0,1))
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta)
plot(evolsel$DeltaBV, x=evolsel$Var)

summary(lm(DeltaBV ~ Sel, data=evolsel))
#### BREEDER'S EQUATION WITH H²=0.17####
evolsel <- read.table("CheckBreederEd_H2_017/EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(corv, ylim=c(-1,1))
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta)
plot(evolsel$DeltaBV, x=evolsel$Var)

library(lme4)
slm <- summary(lm(DeltaBV ~ Sel, data=evolsel))
slm$coefficients
coefficients(lm(DeltaBV ~ Sel, data=evolsel))
summary(lm(DeltaBVRecruits ~ SelRecruits, data=evolsel))

head(evolsel)

## With various h2 ####
h2values <- seq(0.1,1,0.1)
foldernames <- paste0("tryh2eq",h2values)
slopeSE<- data.frame(h2 =h2values, slope=NA,SE=NA, ratio=NA)
for (i in 1:length(foldernames))
{
  evolsel <- read.table(file = paste0(foldernames[i],"/EvolSel.txt"), header=TRUE)
  slm <- summary(lm(DeltaBV ~ Sel, data=evolsel))
  slopeSE[i,2:3] <- slm$coefficients['Sel',c(1,2)]
  slopeSE[i,4] <- mean(evolsel$DeltaBV)/mean(evolsel$Sel)
}
mean(evolsel$DeltaBV/evolsel$Sel)
plot(evolsel$Sel,evolsel$DeltaBV)
abline(0,1)

plot(slopeSE[,1], slopeSE[,2])
arrows(x0=slopeSE[,1], y0 = slopeSE[,2]+1.96*slopeSE[,3], y1=slopeSE[,2]-1.96*slopeSE[,3],code = 3, angle=90)
abline(a = 0, b=1)

plot(slopeSE[,1], slopeSE[,4])
slopeSE

captures <- read.table("tryh2eq1/Captures.txt",header=TRUE)
tapply(X = captures$Age, captures$RUN, function(x) mean(x==1))
mean(captures$Age==1)


es <- read.table(file = "EvolSel.txt", header=TRUE)
mean(es$DeltaBV)/mean(es$Sel)
mean(es$Beta)


#### Analysing power to detect flu evol
library(pedantics)
  library(MCMCglmm)
ALLcaptures <- read.table(file = "ForBetaPower/Captures.txt", header=TRUE)
ALLevolsel <- read.table(file = "ForBetaPower/EvolSel.txt", header=TRUE)
ALLPedigree <- read.table(file = "ForBetaPower/Pedigree.txt", header=TRUE)
Pedigree <- ALLPedigree[ALLPedigree$RUN==1, -1]
names(Pedigree) <- c("animal", "dam", "sire")
Pedigree[Pedigree==0] <- NA
ped<-orderPed(Pedigree)


captures <- ALLcaptures[ALLcaptures$RUN==1,]
evolsel <- ALLevolsel[ALLevolsel$RUN==1,]

plot(evolsel$Sel)
evolsel$Sel>0
YearSelOverMedian <- which(evolsel$Sel>median(evolsel$Sel))

captures$fitness <- captures$Rho+ 2*captures$Phi 
captures$pheno <- captures$BVZ + captures$EVZ
captures$animal <- captures$ID
captures$Z1 <- NA
captures$Z2 <- NA
for (i in 1:nrow(captures))
{
  if (captures$Year[i] %in% YearSelOverMedian)
  {
    captures$Z1[i] <- captures$pheno[i]
  }else{
    captures$Z2[i] <- captures$pheno[i]
  }
}

cov(captures$Z1, captures$fitness, use = "complete.obs")
cov(captures$Z2, captures$fitness, use = "complete.obs")


cov(captures$BVZ[captures$Year %in% YearSelOverMedian], captures$fitness[captures$Year %in% YearSelOverMedian], use = "complete.obs")
cov(captures$BVZ[! captures$Year %in% YearSelOverMedian], captures$fitness[! captures$Year %in% YearSelOverMedian], use = "complete.obs")

priorTwoPeriodsZEXP <- list(G=list(G1=list(V=diag(3), nu=1, alpha.mu=rep(0,3), alpha.V=diag(3)),
                                G2=list(V=diag(3), nu=4, alpha.mu=rep(0,3), alpha.V=diag(3)),
                                G3=list(V=diag(3), nu=4, alpha.mu=rep(0,3), alpha.V=diag(3))),
                         R=list(V=diag(3), nu=4))
priorTwoPeriodsZ <- list(G=list(G1=list(V=diag(3), nu=4),
                                G2=list(V=diag(3), nu=4),
                                G3=list(V=diag(3), nu=4)),
                         R=list(V=diag(3), nu=4))
mcmcTwoPZ<- MCMCglmm(cbind(fitness,Z1,Z2) ~ trait-1+at.level(trait,c(1)):(Sex*Age),
                                 random=~us(trait):animal+us(trait):ID+idh(trait):Year,
                                 rcov=~us(trait):units, family=c("gaussian",rep("gaussian",2)),
                                 prior=priorTwoPeriodsZ,
                                 pedigree=ped,data=captures,verbose=TRUE,nitt=55000,burnin=5000,thin=100)
summary(mcmcTwoPZ)
plot(mcmcTwoPZ)

autocorr(mcmcTwoPZ$VCV)


#### Univariate models####

priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))


mcmcBLUPSBMI0 <- MCMCglmm(pheno ~1,
                          random=~animal+ID+Year,
                          rcov=~units,
                          prior=priorBLUPS0,
                          pedigree=ped,data=captures,verbose=TRUE,nitt=12000,burnin=2000,thin=10,pr=TRUE)
summary(mcmcBLUPSBMI0)
autocorr(mcmcBLUPSBMI0$VCV)

BV<-mcmcBLUPSBMI0$Sol[,grep(pattern = "animal*",x = colnames(mcmcBLUPSBMI0$Sol))]
animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
colnames(BV) <- animalID

pmBV<-data.frame(animalID,posterior.mode(BV))
names(pmBV)<-c("ID","pBV")
mpmBV<-merge(x = pmBV,y = captures,by="ID",all.y=TRUE, all.x = FALSE)
plot(mpmBV$Year,mpmBV$pBV)

BVextend <- BV[,mpmBV$ID]# duplicates posterior distribution for ind present in multiple years

lmBV<-as.mcmc(apply(BVextend,MARGIN = 1,function(x){coef(lm(x~1+mpmBV$Year))[2]}))


library(mgcv)

bvplotlist <- list()
for (i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  gm0 <- gam(bv~1+s(t),data=damdat)
  plotgm0 <- plot.gam(gm0,n = 20)
  bvplotlist[[i]] <- cbind(plotgm0[[1]]$x,plotgm0[[1]]$fit)
}
plot(x=0,ylim=c(-1,1),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2]-x[1,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
abline(h=0)

bvpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVextend), ncol = 10))
names(bvpairwise) <- 1:10
for (i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  tmeanbv <- tapply(damdat$bv,damdat$t,mean)
  bvpairwise[i,] <- tmeanbv[-1]-tmeanbv[-11]
}

boxplot(bvpairwise)
abline(h=0)

evolOver <- as.mcmc(rowMeans(bvpairwise[,which(1:10 %in% YearSelOverMedian)]))
evolUnder <- as.mcmc(rowMeans(bvpairwise[,which(! 1:10 %in% YearSelOverMedian)]))

plot(evolOver-evolUnder)
HPDinterval(evolOver-evolUnder)
