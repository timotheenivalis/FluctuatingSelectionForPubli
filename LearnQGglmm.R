priorBLUPSFitness0<-list(G=list(G1=list(V=1, nu=1.1),G2=list(V=1, nu=1.1),G3=list(V=1, nu=1.1)),
                  R=list(V=1, nu=1.1))

mcmcBLUPSFitness0 <- MCMCglmm(FitnessZ ~Sex*Age+RJst*Age,
                        random=~animal+ID+Year,
                        rcov=~units,
                        prior=priorBLUPSFitness0, family= "poisson",
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100)
summary(mcmcBLUPSFitness0)
toosmall <- 1:1000
posterior.mode(as.mcmc(mcmcBLUPSFitness0$VCV[,"animal"][toosmall]))
HPDinterval(as.mcmc(mcmcBLUPSFitness0$VCV[,"animal"][toosmall]))

save(mcmcBLUPSFitness0, file="mcmcBLUPSFitness0")
mean(YearPheno$FitnessZ)

pred <- mcmcBLUPSFitness0$X%*% as.numeric(mcmcBLUPSFitness0$Sol[1,c("(Intercept)","SexMale","AgeJ","SexMale:AgeJ", "RJst", "AgeJ:RJst")]) 

library(QGglmm)

QGmean(mu = mcmcBLUPSFitness0$Sol[1,"(Intercept)"], var = sum(mcmcBLUPSFitness0$VCV[1,-3]),link.inv = exp, predict=pred)

FitnessParam <- list()
for (i in 1:length(mcmcBLUPSFitness0$Sol[,"(Intercept)"]))
  {
  pred <- mcmcBLUPSFitness0$X%*% as.numeric(mcmcBLUPSFitness0$Sol[i,c("(Intercept)","SexMale","AgeJ","SexMale:AgeJ", "RJst", "AgeJ:RJst")]) 
FitnessParam[[i]] <- QGparams(mu = mcmcBLUPSFitness0$Sol[i,"(Intercept)"],
         var.a = mcmcBLUPSFitness0$VCV[i,"animal"],
         var.p = sum(mcmcBLUPSFitness0$VCV[i,c(1,2)]),
         model = "Poisson.log", 
         predict = as.numeric(pred))
}
HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
mean(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
plot(type="l",as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))

HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))
posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))
mean(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))

plot(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))

toosmall <- which(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))>10^-3 & unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))<100)

posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))
mean(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))
HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))
plot(type="l",as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))


pred <- mcmcBLUPSFitness0$X%*% as.numeric(mcmcBLUPSFitness0$Sol[i,c("(Intercept)","SexMale","AgeJ","SexMale:AgeJ", "RJst", "AgeJ:RJst")]) 
QGparams(mu = mcmcBLUPSFitness0$Sol[i,c("(Intercept)")],
         var.a = mcmcBLUPSFitness0$VCV[i,"animal"],
         var.p = sum(mcmcBLUPSFitness0$VCV[i,c(1,2)]),
         model = "Poisson.log",predict = as.numeric(pred))

var(YearPheno$FitnessZ)
mean(YearPheno$FitnessZ)
