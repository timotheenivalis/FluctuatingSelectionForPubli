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
ALLcaptures <- read.table(file = "ForBetaPower/Captures.txt", header=TRUE)
ALLevolsel <- read.table(file = "ForBetaPower/EvolSel.txt", header=TRUE)
ALLPedigree <- read.table(file = "ForBetaPower/Pedigree.txt", header=TRUE)
Pedigree <- ALLPedigree[ALLPedigree$RUN==1, -1]
names(Pedigree) <- c("animal", "dam", "sire")
Pedigree[Pedigree==0] <- NA
ped<-orderPed(Pedigree)
pedigreeStats(ped)


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

priorTwoPeriodsZ <- list(G=list(G1=list(V=diag(3), nu=100),
                                G2=list(V=diag(3), nu=100),
                                G3=list(V=diag(3), nu=100)),
                         R=list(V=diag(3), nu=100))

mcmcTwoPZ<- MCMCglmm(cbind(fitness,Z1,Z2) ~ trait-1+at.level(trait,c(1)):(Sex*Age),
                                 random=~us(trait):animal+us(trait):ID+idh(trait):Year,
                                 rcov=~us(trait):units, family=c("poisson",rep("gaussian",2)),
                                 prior=priorTwoPeriodsZ,
                                 pedigree=ped,data=captures,verbose=TRUE,nitt=13500,burnin=3500,thin=10)
summary(mcmcTwoPZ)
plot(mcmcTwoPZ)

autocorr(mcmcTwoPZ$VCV)
