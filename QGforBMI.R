library(MCMCglmm)
library(MASS)
YearPheno <- read.table(file = "YearPheno.txt", header=T)
ped <- read.table(file = "ped.txt", header=TRUE)

summary(lmer(formula =BMI ~ 1 + Sex + Age +(1|ID) + (1|Year), data=YearPheno, REML=FALSE))

priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))


mcmcBLUPSBMI0 <- MCMCglmm(BMI ~Sex+Age,#*(RJst+RJ2st)
                        random=~animal+ID+Mother+Year,
                        rcov=~units,
                        prior=priorBLUPS0,
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100,pr=TRUE)
summary(mcmcBLUPSBMI0)
autocorr(mcmcBLUPSBMI0$VCV)
