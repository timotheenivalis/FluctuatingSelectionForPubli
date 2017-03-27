
setwd("C:/Users/Thimothee Admin/Documents/thesis/Mass/FluctuatingSelectionForPubli/")
library(MCMCglmm)
library(MASS)
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)
ped <- read.table(file = "ped.txt", header=TRUE, stringsAsFactors = FALSE)



summary(lmer(formula =BMI ~ 1 + Sex *Age+Age*RJst+(1|ID) + (1|Year), data=YearPheno, REML=FALSE))

priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))


mcmcBLUPSBMI0 <- MCMCglmm(BMI ~Sex*Age+Age*RJst,
                        random=~animal+ID+Mother+Year,
                        rcov=~units,
                        prior=priorBLUPS0,
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100,pr=TRUE)
summary(mcmcBLUPSBMI0)
autocorr(mcmcBLUPSBMI0$VCV)



h2A <- mcmcBLUPSBMI0$VCV[,"animal"]/(mcmcBLUPSBMI0$VCV[,"animal"]+mcmcBLUPSBMI0$VCV[,"ID"]+mcmcBLUPSBMI0$VCV[,"Mother"]+mcmcBLUPSBMI0$VCV[,"units"])
plot(h2A)
posterior.mode(h2A)
HPDinterval(h2A)
posterior.mode(mcmcBLUPSBMI0$VCV[,"animal"])
HPDinterval(mcmcBLUPSBMI0$VCV[,"animal"])

BV<-mcmcBLUPSBMI0$Sol[,grep(pattern = "animal*",x = colnames(mcmcBLUPSBMI0$Sol))]
animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
colnames(BV) <- animalID

pmBV<-data.frame(animalID,posterior.mode(BV))
names(pmBV)<-c("ID","pBV")
mpmBV<-merge(x = pmBV,y = YearPheno,by="ID",all.y=TRUE, all.x = FALSE)
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
plot(x=0,xlim=c(2006,2015),ylim=c(-5,5),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2]-x[1,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
abline(h=0)

plot(x=0,xlim=c(2006,2015),ylim=c(-5,5),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
abline(h=0)


bvpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVextend), ncol = 2015-2006))
names(bvpairwise) <- 2007:2015
for (i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  tmeanbv <- tapply(damdat$bv,damdat$t,mean)
  bvpairwise[i,] <- tmeanbv[-1]-tmeanbv[-10]
}

boxplot(bvpairwise)
abline(h=0)

bvchange <- vector()
for(i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  bvchange[i] <- mean(damdat$bv[damdat$t==2010])-mean(damdat$bv[damdat$t==2006])
}
plot(as.mcmc(bvchange))
HPDinterval(as.mcmc(bvchange))
mean(as.mcmc(bvchange)>0)

bvchange <- vector()
for(i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  bvchange[i] <- mean(damdat$bv[damdat$t==2015])-mean(damdat$bv[damdat$t==2010])
}
plot(as.mcmc(bvchange))
HPDinterval(as.mcmc(bvchange))
mean(as.mcmc(bvchange)<0)


#### Drift simulations ####
pmBV<-data.frame(animalID,posterior.mode(BV))
names(pmBV)<-c("ID","pBV")
mpmBV<-merge(x = pmBV,y = YearPheno,by="ID",all.y=TRUE, all.x = FALSE)
plot(mpmBV$Year,mpmBV$pBV)

BVdtemp <- t(sapply(X = as.numeric(mcmcBLUPSBMI0$VCV[,"animal"]),FUN = function(x) {rbv(pedigree = ped,G = x)} ))
colnames(BVdtemp) <- ped[,1]

BVdextend <- BVdtemp[,mpmBV$ID]# duplicates posterior distribution for ind present in multiple years
bvdpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVdextend), ncol = 2015-2006))
names(bvdpairwise) <- 2007:2015
for (i in 1:nrow(BVdextend))
{
  damdatd<- data.frame(bv=BVdextend[i,],t=mpmBV$Year)
  tmeanbvd <- tapply(damdatd$bv,damdatd$t,mean)
  bvdpairwise[i,] <- tmeanbvd[-1]-tmeanbvd[-10]
}
str(bvdpairwise)
HighDrit <- apply(X = bvdpairwise, MARGIN = 2, function(x){quantile(x, probs = 0.025)})
LowDrift <- apply(X = bvdpairwise, MARGIN = 2, function(x){quantile(x, probs = 0.975)})

polygon(x = c(2006,2008:2014,2016,2016,2014:2008,2006) , y = c(LowDrift,rev(HighDrit)),fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.3), lty=2)


##### Selection gradients
SelGRAByYear <- vector(length = 2015-2006)
SeSelGRAByYear <- vector(length = 2015-2006)
CISelGRAByYear <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  m0 <- glm(FitnessYear ~ 1 + BMIst + Sex *Age + Age*RJst, data=YearPheno[YearPheno$Year==t,], family=quasipoisson, na.action = "na.omit")
  SelGRAByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelGRAByYear[t-2005] <- sm0$coefficients[2,2]
  if(t<2015) {CISelGRAByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "BMIst"))}else{
    CISelGRAByYear[,t-2005] <- c(SelGRAByYear[t-2005]-1.96*SeSelGRAByYear[t-2005],SelGRAByYear[t-2005]+1.96*SeSelGRAByYear[t-2005])
  }
}


#### Correlation between selection and evolution#### 

SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {cor(x,ZSelAByYear)})
hist(SelToG)
plot(as.mcmc(SelToG))
posterior.mode(as.mcmc(SelToG))
HPDinterval(as.mcmc(SelToG))
mean(SelToG<0)

SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {coefficients(lm(x~1+SelGRAByYear[-10]))[2]})


SelPhiToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {cor(x[-9],SelGRAByYear[-c(1,10)])})
hist(SelPhiToG)
plot(as.mcmc(SelPhiToG))
posterior.mode(as.mcmc(SelPhiToG))
HPDinterval(as.mcmc(SelToG))
mean(SelToG<0)
plot(x=2006:2014,SelAByYearPhi[-10])
plot(x=2007:2015,y=bvpairwise[1,])
plot(x=SelAByYearPhi[-10],y=bvpairwise[1,])

rbind(SelAByYearPhi,bvpairwise[1,])
plot(colMeans(bvpairwise),SelAByYearPhi[-c(1)])


vAF <- 0.1
vAM <- posterior.mode(mcmcBLUPSBMI0$VCV[,"animal"])
Gcovvar <- diag(3)*(vAM)
Gcovvar[1,1] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSBMI0$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSBMI0$VCV[,"ID"]),mean(mcmcBLUPSBMI0$VCV[,"ID"]))
Icovvar <- diag(3)*(vIM)
Icovvar[1,1] <- vIF

vMF <- 0.1
vMM <- ifelse(posterior.mode(mcmcBLUPSBMI0$VCV[,"Mother"])>=0, posterior.mode(mcmcBLUPSBMI0$VCV[,"Mother"]),mean(mcmcBLUPSBMI0$VCV[,"Mother"]))
Mcovvar <- diag(3)*(vMM)
Mcovvar[1,1] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSBMI0$VCV[,"units"])
Rcovvar <- diag(3)*(vRM)
Rcovvar[1,1] <- vRF

priorTwoPeriodsA <- list(G=list(G1=list(V=Gcovvar, nu=1),
                                G2=list(V=Icovvar, nu=1),
                                G3=list(V=Mcovvar, nu=1),
                                G4=list(V=diag(1), nu=0.001)),
                         R=list(V=Rcovvar, nu=1))

mcmcBivTwoPeriods <- MCMCglmm(cbind(Fitness,BMI1,BMI2) ~ trait-1+at.level(trait,c(1)):(Sex*Age+RJst)+at.level(trait,c(2:3)):(Sex*Age+Age*RJst),
                              random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(1))):Year,
                              rcov=~us(trait):units, family=c("poisson",rep("gaussian",2)),
                              prior=priorTwoPeriodsA,
                              pedigree=ped,data=YearPheno,verbose=TRUE,nitt=650000,burnin=150000,thin=500)

summary(mcmcBivTwoPeriods)
plot(mcmcBivTwoPeriods)

save.image("~/thesis/Mass/FluctuatingSelectionForPubli/EnvQGBMI.RData")

BEBMI1 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI1.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI1.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI1.units"]
plot(BEBMI1)

BEBMI2 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI2.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI2.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitBMI2.units"]
plot(BEBMI2)

plot(BEBMI1 -BEBMI2)

###SURVIVAL


vAF <- 0.1
vAM <- posterior.mode(mcmcBLUPSBMI0$VCV[,"animal"])
Gcovvar <- diag(3)*(vAM)
Gcovvar[3,3] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSBMI0$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSBMI0$VCV[,"ID"]),mean(mcmcBLUPSBMI0$VCV[,"ID"]))
Icovvar <- diag(3)*(vIM)
Icovvar[3,3] <- vIF

vMF <- 0.1
vMM <- ifelse(posterior.mode(mcmcBLUPSBMI0$VCV[,"Mother"])>=0, posterior.mode(mcmcBLUPSBMI0$VCV[,"Mother"]),mean(mcmcBLUPSBMI0$VCV[,"Mother"]))
Mcovvar <- diag(3)*(vMM)
Mcovvar[3,3] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSBMI0$VCV[,"units"])
Rcovvar <- diag(3)*(vRM)
Rcovvar[3,3] <- 1

priorTwoPeriodsPHI <- list(G=list(G1=list(V=Gcovvar, nu=1),
                                G2=list(V=Icovvar, nu=1),
                                G3=list(V=Mcovvar, nu=1),
                                G4=list(V=diag(3), nu=1)),
                         R=list(V=Rcovvar, nu=1, fix=3))

mcmcBivTwoPeriodsPHI <- MCMCglmm(cbind(BMIPhi1,BMIPhi2,Phi) ~ trait-1+at.level(trait,c(3)):(Sex*Age+RJst)+at.level(trait,c(1:2)):(Sex*Age+Age*RJst),
                              random=~us(trait):animal+us(trait):ID+us(trait):Mother+idh(trait):Year,
                              rcov=~us(trait):units, family=c(rep("gaussian",2), "categorical"),
                              prior=priorTwoPeriodsPHI,
                              pedigree=ped,data=YearPheno,verbose=TRUE,nitt=650000,burnin=150000,thin=500)

summary(mcmcBivTwoPeriodsPHI)
save.image("~/thesis/Mass/FluctuatingSelectionForPubli/EnvQGBMI.RData")


###Repro
mcmcBivTwoPeriodsRHO <- MCMCglmm(cbind(Rho,BMIRho1,BMIRho2) ~ trait-1+at.level(trait,c(1)):(Sex*Age+RJst)+at.level(trait,c(2:3)):(Sex*Age+Age*RJst),
                              random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(1))):Year,
                              rcov=~us(trait):units, family=c("poisson",rep("gaussian",2)),
                              prior=priorTwoPeriodsA,
                              pedigree=ped,data=YearPheno,verbose=TRUE,nitt=650000,burnin=150000,thin=500)

summary(mcmcBivTwoPeriods)
save.image("~/thesis/Mass/FluctuatingSelectionForPubli/EnvQGBMI.RData")


#### MODEL WITH GENETIC GROUPS ####
library(MCMCglmm)
library(nadiv)


priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))


mcmcBLUPSBMI0_GG <- MCMCglmm(BMI ~Sex*Age+Age*RJst + GGImm,
                          random=~animal+ID+Mother+Year,
                          rcov=~units,
                          prior=priorBLUPS0,
                          pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100,pr=TRUE)
summary(mcmcBLUPSBMI0_GG)

BV<-mcmcBLUPSBMI0_GG$Sol[,grep(pattern = "animal*",x = colnames(mcmcBLUPSBMI0_GG$Sol))]
animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
colnames(BV) <- animalID

pmBV<-data.frame(animalID,posterior.mode(BV))
names(pmBV)<-c("ID","pBV")
mpmBV<-merge(x = pmBV,y = YearPheno,by="ID",all.y=TRUE, all.x = FALSE)
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
plot(x=0,xlim=c(2006,2016),ylim=c(-5,5),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2]-x[1,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
abline(h=0)

plot(x=0,xlim=c(2006,2016),ylim=c(-5,5),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
abline(h=0)

bvpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVextend), ncol = 2016-2006))
names(bvpairwise) <- 2007:2016
for (i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  tmeanbv <- tapply(damdat$bv,damdat$t,mean)
  bvpairwise[i,] <- tmeanbv[-1]-tmeanbv[-11]
}

boxplot(bvpairwise)
abline(h=0)

bvchange <- vector()
for(i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  bvchange[i] <- mean(damdat$bv[damdat$t==2014])-mean(damdat$bv[damdat$t==2006])
}
plot(as.mcmc(bvchange))
HPDinterval(as.mcmc(bvchange))
mean(as.mcmc(bvchange)>0)

bvrebound <- vector()
bvrebounddiff <- vector()
for(i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  lmrebound <- lm(bv ~ 1 + t + I(t^2), data=damdat)
  bvrebound[i] <- coefficients(lmrebound)[3]
  bvrebounddiff[i] <- mean(damdat$bv[damdat$t==2016])-mean(damdat$bv[damdat$t==2014])
}
plot(bvrebounddiff)
mean(bvrebounddiff<0)


#### Correlation selection evolution ####
SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {cor(x,ZSelAByYear)})
hist(SelToG)
plot(as.mcmc(SelToG))
posterior.mode(as.mcmc(SelToG))
HPDinterval(as.mcmc(SelToG))
mean(SelToG<0)

SelDiffYear <- vector(length = 2015-2006)
for (t in 2006:2015)
{
  m0 <- glm(FitnessZ ~ 1 + BMIst + Sex *Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson, na.action = "na.omit")
  SelDiffYear[t-2005] <- coefficients(m0)[2]*var(YearPheno[YearPheno$Year==t,"BMIst"], na.rm = TRUE)
}

SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {cor(x,SelDiffYear)})
hist(SelToG)
plot(as.mcmc(SelToG))
posterior.mode(as.mcmc(SelToG))
HPDinterval(as.mcmc(SelToG))
mean(SelToG<0)


vAF <- 0.1
vAM <- posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"animal"])
Gcovvar <- diag(3)*(vAM)
Gcovvar[1,1] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"ID"]),mean(mcmcBLUPSBMI0$VCV[,"ID"]))
Icovvar <- diag(3)*(vIM)
Icovvar[1,1] <- vIF

vMF <- 0.1
vMM <- ifelse(posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"Mother"])>=0, posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"Mother"]),mean(mcmcBLUPSBMI0_GG$VCV[,"Mother"]))
Mcovvar <- diag(3)*(vMM)
Mcovvar[1,1] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSBMI0_GG$VCV[,"units"])
Rcovvar <- diag(3)*(vRM)
Rcovvar[1,1] <- vRF

priorTwoPeriodsA <- list(G=list(G1=list(V=Gcovvar, nu=1),
                                G2=list(V=Icovvar, nu=1),
                                G3=list(V=Mcovvar, nu=1),
                                G4=list(V=diag(1), nu=0.001)),
                         R=list(V=Rcovvar, nu=1))

mcmcBivTwoPeriods_GG <- MCMCglmm(cbind(FitnessZ,BMI1,BMI2) ~ trait-1+at.level(trait,c(1)):(Sex*Age+RJst + GGImm)+at.level(trait,c(2:3)):(Sex*Age+Age*RJst + GGImm),
                              random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(1))):Year,
                              rcov=~us(trait):units, family=c("poisson",rep("gaussian",2)),
                              prior=priorTwoPeriodsA,
                              pedigree=ped,data=YearPheno,verbose=TRUE,nitt=650000,burnin=150000,thin=500)

summary(mcmcBivTwoPeriods_GG)
plot(mcmcBivTwoPeriods_GG)

save(mcmcBivTwoPeriods_GG, file = "mcmcBivTwoPeriods_GG")


load("mcmcBivTwoPeriods_GG")

corGM <- mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI1.animal"]/sqrt(mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.animal"]*mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.animal"])
HPDinterval(as.mcmc(corGM))
posterior.mode(corGM)

BetaG1 <-  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.animal"]/mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.animal"]
HPDinterval(BetaG1)
posterior.mode(BetaG1)

BetaG2 <-  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.animal"]/mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.animal"]
HPDinterval(BetaG2)
posterior.mode(BetaG2)

covE1 <- mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.ID"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.Mother"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.units"]
plot(covE1)
BetaE1 <- covE1 /(mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.ID"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.Mother"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.units"])
plot(BetaE1)

covE2 <- mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.ID"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.Mother"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.units"]
plot(covE2)
BetaE2 <- covE2 /(mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.ID"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.Mother"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.units"])
plot(BetaE2)

diffBeta1 <- BetaG1 - BetaE1
plot(diffBeta1)
HPDinterval(diffBeta1)
mean(diffBeta1>0)

diffBeta2 <- BetaG2 - BetaE2
plot(diffBeta2)

covP1 <- mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.animal"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.ID"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.Mother"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI1.units"]
plot(covP1)
BetaP1 <- covP1 /(mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.animal"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.ID"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.Mother"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI1:traitBMI1.units"])
plot(BetaP1)

covP2 <- mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.animal"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.ID"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.Mother"]+
  mcmcBivTwoPeriods_GG$VCV[,"traitFitnessZ:traitBMI2.units"]
plot(covP2)
BetaP2 <- covP2 /(mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.animal"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.ID"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.Mother"]+
                    mcmcBivTwoPeriods_GG$VCV[,"traitBMI2:traitBMI2.units"])
plot(BetaP2)
