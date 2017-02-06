# offspring size VS. number
plot(YearPheno$OS, YearPheno$Rho)
summary(glm(Rho ~ 1 + Sex+OS, data=YearPheno[YearPheno$Age=="A",], family="poisson"))
summary(glm(Fitness ~ 1 + Sex+OS+I(OS^2), data=YearPheno[YearPheno$Age=="A",], family="poisson"))

priorSizeNb <- list(G=list(G1=list(V=diag(2), nu=1),
                           G2=list(V=diag(2), nu=1),
                           G3=list(V=diag(2), nu=1),
                           G4=list(V=diag(2), nu=0.001)),
                    R=list(V=diag(2), nu=1))

mcmcSizeNb<- MCMCglmm(cbind(Rho,OS) ~ trait-1+Sex+RJst+RJ2st+at.level(trait,c(1)):(Sex+RJst+RJ2st),
                      random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(trait):Year,
                      rcov=~us(trait):units, family=c("poisson","gaussian"),
                      prior=priorSizeNb,
                      pedigree=ped,data=YearPheno[YearPheno$Age=="A",],verbose=TRUE,nitt=26000,burnin=6000,thin=20)

summary(mcmcSizeNb)


# the better model
YearPheno$RhoAd <- YearPheno$Rho* ifelse(test = YearPheno$Age=="A",1,NA)
YearPheno$AJuv <- YearPheno$A* ifelse(test = YearPheno$Age=="J",1,NA)

priorSizeNb2 <- list(G=list(G1=list(V=diag(2), nu=0.001),
                            G2=list(V=diag(2), nu=0.001),
                            G3=list(V=diag(2), nu=0.001)),
                     R=list(V=diag(2), nu=0.001))
mcmcSizeNb2<- MCMCglmm(cbind(RhoAd,AJuv) ~ trait-1+trait:Sex,
                       random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(trait):Year,
                       rcov=~idh(trait):units, family=c("poisson","gaussian"),
                       prior=priorSizeNb2,
                       pedigree=ped,data=YearPheno,verbose=TRUE,nitt=130000,burnin=30000,thin=100)

summary(mcmcSizeNb2)
plot(mcmcSizeNb2$VCV)