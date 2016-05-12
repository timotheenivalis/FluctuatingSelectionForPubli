priorBLUPSFitness0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001),G4=list(V=0.1, nu=0.0001)),
                  R=list(V=0.1, nu=0.0001))

mcmcBLUPSFitness0 <- MCMCglmm(Fitness ~Sex*Age,
                        random=~animal+ID+Mother+Year,
                        rcov=~units,
                        prior=priorBLUPSFitness0, family= "poisson",
                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100)
summary(mcmcBLUPSFitness0)

posterior.mode(as.mcmc(mcmcBLUPSFitness0$VCV[,"animal"][toosmall]))
HPDinterval(as.mcmc(mcmcBLUPSFitness0$VCV[,"animal"][toosmall]))

save(mcmcBLUPSFitness0, file="mcmcBLUPSFitness0")
mean(YearPheno$Fitness)

pred <- mcmcBLUPSFitness0$X%*% as.numeric(mcmcBLUPSFitness0$Sol[1,c("(Intercept)","SexMale","AgeJ","SexMale:AgeJ")]) 

QGmean(mu = mcmcBLUPSFitness0$Sol[1,"(Intercept)"], var = sum(mcmcBLUPSFitness0$VCV[1,]),link.inv = exp, predict=pred)

FitnessParam <- list()
for (i in 1:length(mcmcBLUPSFitness0$Sol[,"(Intercept)"]))
FitnessParam[[i]] <- QGparams(mu = mcmcBLUPSFitness0$Sol[i,"(Intercept)"],
         var.a = mcmcBLUPSFitness0$VCV[i,"animal"],
         var.p = sum(mcmcBLUPSFitness0$VCV[i,-4]),
         model = "Poisson.log", 
         predict = as.numeric(pred))
HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
mean(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))
plot(type="l",as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))))

HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))
posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))
mean(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))

plot(as.mcmc(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))[toosmall]))

toosmall <- which(unlist(lapply(FitnessParam,function(x){x["var.a.obs"]}))>10^-3)

posterior.mode(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))
HPDinterval(as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))
plot(type="l",as.mcmc(unlist(lapply(FitnessParam,function(x){x["h2.obs"]}))[toosmall]))

