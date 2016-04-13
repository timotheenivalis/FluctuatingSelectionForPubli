YearPheno <- read.table(file = "YearPheno.txt", header=T)

summary(lmer(Mass ~ 1+Sex*Age + (1|Year), data=YearPheno))

summary(glm(Phi ~ 1 + Mass + Sex , data=YearPheno[YearPheno$Age=="A" & YearPheno$Year<max(YearPheno$Year),], family=binomial))
summary(glmer(Phi ~ 1 + as.factor(Year) + Mass + Sex +(0+Mass|Year), data=YearPheno[YearPheno$Age=="A" YearPheno$Year<max(YearPheno$Year),], family=binomial))

summary(glm(Rho ~ 1 + StMass + Sex , data=YearPheno, family=poisson))
summary(glmer(Rho ~ 1 + as.factor(Year) + StMass + Sex +(0+Mass|Year), data=YearPheno, family=poisson))


summary(glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=poisson))
summary(glmer(Fitness ~ 1 + as.factor(Year) + StMass + Sex +(0+Mass|Year), data=YearPheno, family=poisson))
summary(glmer(Fitness ~ 1 + StMass + Sex + Age +(1+Mass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson))

prior0 <- list(G=list(G1= list(V=0.00001*diag(2), nu=2)), R=list(V=diag(1),nu=1))
m0 <- MCMCglmm(fixed = Fitness ~ 1 + Sex + Age + StMass,
               random = ~ us(1+StMass):Year, prior= prior0, family="poisson",
               data=YearPheno[!is.na(YearPheno$StMass),], nitt = 110000, burnin = 10000, thin = 100)
summary(m0)
autocorr(m0$VCV)
plot(m0)
