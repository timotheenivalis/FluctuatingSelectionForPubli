library(MCMCglmm)
library(MASS)
YearPheno <- read.table(file = "YearPheno.txt", header=T)
ped <- read.table(file = "ped.txt", header=TRUE)

priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))

mcmcBLUPSA0 <- MCMCglmm(A ~Sex+Age,#*(RJst+RJ2st)
                        random=~animal+ID+Mother+Year,
                        rcov=~units,
                        prior=priorBLUPS0,
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100,pr=TRUE)
summary(mcmcBLUPSA0)
autocorr(mcmcBLUPSA0$VCV)

save(mcmcBLUPSA0, file = "mcmcBLUPSA0")

h2A <- mcmcBLUPSA0$VCV[,"animal"]/(mcmcBLUPSA0$VCV[,"animal"]+mcmcBLUPSA0$VCV[,"ID"]+mcmcBLUPSA0$VCV[,"Mother"]+mcmcBLUPSA0$VCV[,"units"])
plot(h2A)
posterior.mode(h2A)
HPDinterval(h2A)
posterior.mode(mcmcBLUPSA0$VCV[,"animal"])
HPDinterval(mcmcBLUPSA0$VCV[,"animal"])

BV<-mcmcBLUPSA0$Sol[,grep(pattern = "animal*",x = colnames(mcmcBLUPSA0$Sol))]
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
plot(x=0,xlim=c(2006,2015),ylim=c(-2,2),type="n")
trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})

bvpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVextend), ncol = 2015-2006))
names(bvpairwise) <- 2007:2015
for (i in 1:nrow(BVextend))
{
  damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
  lm(bv~1+t,data=damdat[damdat$t==2006 | damdat$t==2007,])
  mean(damdat$bv[damdat$t==2007])-mean(damdat$bv[damdat$t==2006])
  tmeanbv <- tapply(damdat$bv,damdat$t,mean)
  bvpairwise[i,] <- tmeanbv[-1]-tmeanbv[-10]
}

boxplot(bvpairwise)
abline(h=0)

##### Selection gradients
SelGRAByYear <- vector(length = 2015-2006)
SeSelGRAByYear <- vector(length = 2015-2006)
CISelGRAByYear <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  m0 <- glm(FitnessYear ~ 1 + Ast + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson)
  SelGRAByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelGRAByYear[t-2005] <- sm0$coefficients[2,2]
  if(t<2015) {CISelGRAByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "Ast"))}else{
    CISelGRAByYear[,t-2005] <- c(SelGRAByYear[t-2005]-1.96*SeSelGRAByYear[t-2005],SelGRAByYear[t-2005]+1.96*SeSelGRAByYear[t-2005])
  }
}

SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {cor(x,SelGRAByYear[-10])})
hist(SelToG)
plot(as.mcmc(SelToG))
posterior.mode(as.mcmc(SelToG))
HPDinterval(as.mcmc(SelToG))
mean(SelToG<0)

SelToG <- apply(bvpairwise,MARGIN = 1, FUN = function(x) {coefficients(lm(x~1+SelGRAByYear[-10]))[2]})

#### BIV ANIMAL MODEL ALL YEARS #### (SCARY) # For A ####

nuVM <- 0.001
nuVF <- 0.001
nucovMF <- 0.001
spePriornu <- diag(11)*(nuVM-nuVF)+nuVF
spePriornu[1,] <- nucovMF
spePriornu[,1] <- nucovMF

vAF <- 0.1
vAM <- posterior.mode(mcmcBLUPSA0$VCV[,"animal"])
Gcovvar <- diag(11)*(vAM)
Gcovvar[1,1] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSA0$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSA0$VCV[,"ID"]),mean(mcmcBLUPSA0$VCV[,"ID"]))
Icovvar <- diag(11)*(vIM)
Icovvar[1,1] <- vIF

vMF <- 0.1
vMM <- posterior.mode(mcmcBLUPSA0$VCV[,"Mother"])
Mcovvar <- diag(11)*(vMM)
Mcovvar[1,1] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSA0$VCV[,"units"])
Rcovvar <- diag(11)*(vRM)
Rcovvar[1,1] <- vRF

priorAllYearsA <- list(G=list(G1=list(V=Gcovvar, nu=spePriornu, fix =2),
                              G2=list(V=Icovvar, nu=spePriornu, fix =2),
                              G3=list(V=Mcovvar, nu=spePriornu, fix =2),
                              G4=list(V=diag(1), nu=0.001)),
                       R=list(V=Rcovvar, nu=spePriornu, fix =2))

mcmcBivAllYears <- MCMCglmm(cbind(Fitness,A2006,A2007,A2008,A2009,A2010,A2011,A2012,A2013,A2014,A2015) ~ trait-1+Sex+Age*(RJst+RJ2st)+at.level(trait,c(1)):(Sex+Age*(RJst+RJ2st)),
                            random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(1))):Year,
                            rcov=~us(trait):units, family=rep("gaussian",11),
                            prior=priorAllYearsA,
                            pedigree=ped,data=YearPheno,verbose=TRUE,nitt=13,burnin=3,thin=1)
summary(mcmcBivAllYears)
save(mcmcBivAllYears,file = "mcmcBivAllYearsA0")

nuVM <- 0.001
nuVF <- 0.001
nucovMF <- 0.001
spePriornu <- diag(11)*(nuVM-nuVF)+nuVF
spePriornu[1,] <- nucovMF
spePriornu[,1] <- nucovMF

vAF <- 0.1
vAM <- posterior.mode(mcmcBLUPSA0$VCV[,"animal"])
Gcovvar <- diag(11)*(vAM)
Gcovvar[1,1] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSA0$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSA0$VCV[,"ID"]),mean(mcmcBLUPSA0$VCV[,"ID"]))
Icovvar <- diag(11)*(vIM)
Icovvar[1,1] <- vIF

vMF <- 0.1
vMM <- posterior.mode(mcmcBLUPSA0$VCV[,"Mother"])
Mcovvar <- diag(11)*(vMM)
Mcovvar[1,1] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSA0$VCV[,"units"])
Rcovvar <- diag(11)*(vRM)
Rcovvar[1,1] <- vRF

priorAllYearsA <- list(G=list(G1=list(V=Gcovvar, nu=spePriornu, fix =2 ),
                              G2=list(V=Mcovvar, nu=11),
                              G3=list(V=diag(1), nu=0.001)),
                       R=list(V=Rcovvar, nu=spePriornu, fix =2))

mcmcBivAllYears <- MCMCglmm(cbind(Fitness,A2006,A2007,A2008,A2009,A2010,A2011,A2012,A2013,A2014,A2015) ~ trait-1+Sex+Age*(RJst+RJ2st)+at.level(trait,c(1)):(Sex+Age*(RJst+RJ2st)),
                            random=~us(trait):animal+idh(trait):Mother+us(at.level(trait,c(1))):Year,
                            rcov=~us(trait):units, family=c("poisson",rep("gaussian",10)),
                            prior=priorAllYearsA,
                            pedigree=ped,data=YearPheno,verbose=TRUE,nitt=130000,burnin=30000,thin=100)
summary(mcmcBivAllYears)
save(mcmcBivAllYears,file = "mcmcBivAllYearsA0")
plot(mcmcBivAllYears)

## unstructured
nuVM <- 10^8
nuVF <- 1
nucovMF <- 1
spePriornu <- diag(11)*(nuVM-nuVF)+nuVF
spePriornu[1,-1] <- nucovMF
spePriornu[-1,1] <- nucovMF

priorAllYearsAexp <- list(G=list(G1=list(V=Gcovvar, nu=spePriornu, fix =2, alpha.mu =diag(x = Gcovvar), alpha.V = Gcovvar*1000),
                                 G2=list(V=Mcovvar, nu=1, alpha.mu = diag(x = Mcovvar), alpha.V = Mcovvar*1000),
                                 G3=list(V=diag(1), nu=0.001)),
                          R=list(V=Rcovvar, nu=spePriornu, fix =2))

mcmcBivAllYears2 <- MCMCglmm(cbind(Fitness,A2006,A2007,A2008,A2009,A2010,A2011,A2012,A2013,A2014,A2015) ~ trait-1+Sex+Age*(RJst+RJ2st)+at.level(trait,c(1)):(Sex+Age*(RJst+RJ2st)),
                             random=~us(trait):animal+idh(trait):Mother+us(at.level(trait,c(1))):Year,
                             rcov=~us(trait):units, family=c("poisson",rep("gaussian",10)),
                             prior=priorAllYearsAexp,
                             pedigree=ped,data=YearPheno,verbose=TRUE,nitt=1300000,burnin=300000,thin=1000)

save(mcmcBivAllYears2,file = "mcmcBivAllYearsA2")
summary(mcmcBivAllYears2)
plot(mcmcBivAllYears2$VCV)

load("mcmcBivAllYearsA2")
colnames(mcmcBivAllYears2$VCV)
sg <- (mcmcBivAllYears2$VCV[,2:11])
sr <- (mcmcBivAllYears2$VCV[,135:144])

vg <- (mcmcBivAllYears2$VCV[,c(13,25,37,49,61,73,85,97,109,121)])
vr <- (mcmcBivAllYears2$VCV[,c(146,158,170,182,194,206,218,230,242,254)])
 
bg <- sg/vg                           
be <- (sr+sg)/(vr+vg)

diffb <- bg - be
plot(diffb)


priorAllYearsAexpnoG <- list(G=list(G1=list(V=diag(11), nu=1, fix=2),
                                    G2=list(V=diag(11), nu=1),
                                 G2=list(V=diag(1), nu=0.001)),
                          R=list(V=diag(11), nu=1, fix =2))

mcmcBivAllYearsNoG <- MCMCglmm(cbind(Fitness,A2006,A2007,A2008,A2009,A2010,A2011,A2012,A2013,A2014,A2015) ~ trait-1+Sex+Age*(RJst+RJ2st)+at.level(trait,c(1)):(Sex+Age*(RJst+RJ2st)),
                             random=~us(trait):ID +idh(trait):Mother+us(at.level(trait,c(1))):Year,
                             rcov=~us(trait):units, family=c("poisson",rep("gaussian",10)),
                             prior=priorAllYearsAexpnoG,
                             data=YearPheno,verbose=TRUE,nitt=13000,burnin=3000,thin=10)
summary(mcmcBivAllYearsNoG)

#### Based on two groups of years #####

vAF <- 0.1
vAM <- posterior.mode(m0mcBLUPSA0$VCV[,"animal"])
Gcovvar <- diag(3)*(vAM)
Gcovvar[1,1] <- vAF

vIF <- 0.1
vIM <- ifelse(posterior.mode(mcmcBLUPSA0$VCV[,"ID"])>=0, posterior.mode(mcmcBLUPSA0$VCV[,"ID"]),mean(mcmcBLUPSA0$VCV[,"ID"]))
Icovvar <- diag(3)*(vIM)
Icovvar[1,1] <- vIF

vMF <- 0.1
vMM <- posterior.mode(mcmcBLUPSA0$VCV[,"Mother"])
Mcovvar <- diag(3)*(vMM)
Mcovvar[1,1] <- vMF

vRF <- 0.1
vRM <- posterior.mode(mcmcBLUPSA0$VCV[,"units"])
Rcovvar <- diag(3)*(vRM)
Rcovvar[1,1] <- vRF

priorTwoPeriodsA <- list(G=list(G1=list(V=Gcovvar, nu=1),
                                G2=list(V=Icovvar, nu=1),
                                G3=list(V=Mcovvar, nu=1),
                                G4=list(V=diag(1), nu=0.001)),
                         R=list(V=Rcovvar, nu=1))

mcmcBivTwoPeriods <- MCMCglmm(cbind(Fitness,A1,A2) ~ trait-1+Sex+Age*(RJst+RJ2st)+at.level(trait,c(1)):(Sex+Age*(RJst+RJ2st)),
                              random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(1))):Year,
                              rcov=~us(trait):units, family=c("poisson",rep("gaussian",2)),
                              prior=priorTwoPeriodsA,
                              pedigree=ped,data=YearPheno,verbose=TRUE,nitt=650000,burnin=150000,thin=500)

summary(mcmcBivTwoPeriods)
save(mcmcBivTwoPeriods,file="mcmcBivTwoPeriods")

load("mcmcBivTwoPeriods")

corGM <- mcmcBivTwoPeriods$VCV[,"traitA2:traitA1.animal"]/sqrt(mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.animal"]*mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.animal"])
HPDinterval(corGM)
posterior.mode(corGM)

BetaG1 <-  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.animal"]/mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.animal"]
HPDinterval(BetaG1)
posterior.mode(BetaG1)

BetaG2 <-  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.animal"]/mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.animal"]
HPDinterval(BetaG2)
posterior.mode(BetaG2)

covE1 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.units"]
plot(covE1)
BetaE1 <- covE1 /(mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.ID"]+
                    mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.Mother"]+
                    mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.units"])
plot(BetaE1)

covE2 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.units"]
plot(covE2)
BetaE2 <- covE2 /(mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.ID"]+
                    mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.Mother"]+
                    mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.units"])
plot(BetaE2)

diffBeta1 <- BetaG1 - BetaE1
plot(diffBeta1)
HPDinterval(diffBeta1)
mean(diffBeta1>0)

diffBeta2 <- BetaG2 - BetaE2
plot(diffBeta2)

covP1 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.animal"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA1.units"]
plot(covP1)
BetaP1 <- covP1 /(mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.animal"]+
                    mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.ID"]+
                    mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.Mother"]+
                    mcmcBivTwoPeriods$VCV[,"traitA1:traitA1.units"])
plot(BetaP1)

covP2 <- mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.animal"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.ID"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.Mother"]+
  mcmcBivTwoPeriods$VCV[,"traitFitness:traitA2.units"]
plot(covP2)
BetaP2 <- covP2 /(mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.animal"]+
                    mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.ID"]+
                    mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.Mother"]+
                    mcmcBivTwoPeriods$VCV[,"traitA2:traitA2.units"])
plot(BetaP2)


  