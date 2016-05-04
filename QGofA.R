
library(MCMCglmm)
library(MASS)
YearPheno <- read.table(file = "YearPheno.txt", header=T)
ped <- read.table(file = "ped.txt", header=TRUE)

priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))

mcmcBLUPSA0 <- MCMCglmm(A ~Sex+Age*(RJst+RJ2st),
                        random=~animal+ID+Mother+Year,
                        rcov=~units,
                        prior=priorBLUPS0,
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100,pr=TRUE)
summary(mcmcBLUPSA0)
autocorr(mcmcBLUPSA0$VCV)

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
