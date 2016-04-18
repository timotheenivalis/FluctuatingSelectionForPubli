library(MASS)

m0<-glm(Fitness ~ 1 +Age , data=YearPheno[YearPheno$Year==2015,], family = quasipoisson)
summary(m0)
tapply(YearPheno[YearPheno$Year==2015,"Fitness"], YearPheno[YearPheno$Year==2015,"Age"], mean)
confint(object = m0, method="boot")
YearPheno <- read.table(file = "YearPheno.txt", header=T)

SelByYear <- vector(length = 2015-2006)
SeSelByYear <- vector(length = 2015-2006)
CISelByYear <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
  {
  m0 <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson)
  SelByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelByYear[t-2005] <- sm0$coefficients[2,2]
  if(t<2015) {CISelByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "StMass"))}else{
    CISelByYear[,t-2005] <- c(SelByYear[t-2005]-1.96*SeSelByYear[t-2005],SelByYear[t-2005]+1.96*SeSelByYear[t-2005])
  }
}
plot(SelByYear, x=2006:2015, ylim=c(-1.5,0.8), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYear[1,],
       y1 = CISelByYear[2,], angle = 90,length = 0.1)
m0all <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=quasipoisson)
abline(h=coefficients(m0all)[2], lty=2)
sm0all <- summary(m0all)
lowm0all <- coefficients(m0all)[2]+1.96*sm0all$coefficients[2,2]
highm0all <- coefficients(m0all)[2]-1.96*sm0all$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0all,lowm0all, highm0all, highm0all),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)




summary(glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=poisson))

mmRRfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1+Mass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)
mmRIfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)
mmRnoCorfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1|Year) + (0+StMass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)

summary(mmRRfitness)
summary(mmRnoCorfitness)
logLik(mmRRfitness)
logLik(mmRIfitness)
anova(mmRIfitness,mmRRfitness)
anova(mmRIfitness,mmRnoCorfitness)
anova(mmRRfitness,mmRnoCorfitness)
  CImmRnoCorfitness <- confint(mmRnoCorfitness)

points(x=2006:2015,y=unlist(coefficients(mmRnoCorfitness)$Year["StMass"]))



prior0 <- list(G=list(G1= list(V=0.00001*diag(2), nu=2)), R=list(V=diag(1),nu=1))
m0 <- MCMCglmm(fixed = Fitness ~ 1 + Sex + Age + StMass,
               random = ~ us(1+StMass):Year, prior= prior0, family="poisson",
               data=YearPheno[!is.na(YearPheno$StMass),], nitt = 110000, burnin = 10000, thin = 100)
summary(m0)
autocorr(m0$VCV)
plot(m0)


##### Fitness components ####

SelByYearPhi <- vector(length = 2015-2005)
SeSelByYearPhi <- vector(length = 2015-2005)
CISelByYearPhi <- matrix(NA,nrow=2,ncol=2015-2005)

for (t in 2006:2014)
{
  m0 <- glm(Phi ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=binomial)
  SelByYearPhi[t-2005] <- coefficients(m0)[2]
  CISelByYearPhi[,t-2005]<-as.numeric(confint(m0,parm="StMass"))
  sm0<-summary(m0)
  SeSelByYearPhi[t-2005] <- sm0$coefficients[2,2]
}
SelByYearPhi[10] <- NA
SeSelByYearPhi[10] <- NA
plot(SelByYearPhi, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYearPhi[1,],
       y1 = CISelByYearPhi[2,], angle = 90,length = 0.1)
m0allphi <- glm(Phi ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year<2015,], family=binomial)
abline(h=coefficients(m0allphi)[2], lty=2)
sm0allphi <- summary(m0allphi)
lowm0allphi <- coefficients(m0allphi)[2]+1.96*sm0allphi$coefficients[2,2]
highm0allphi <- coefficients(m0allphi)[2]-1.96*sm0allphi$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allphi,lowm0allphi, highm0allphi, highm0allphi),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

summary(glm(Phi ~ 1 + Mass + Sex , data=YearPheno[YearPheno$Age=="A" & YearPheno$Year<max(YearPheno$Year),], family=binomial))
summary(glmer(Phi ~ 1 + as.factor(Year) + Mass + Sex +(0+Mass|Year), data=YearPheno[YearPheno$Age=="A"& YearPheno$Year<max(YearPheno$Year),], family=binomial))

summary(glmer(Phi ~ 1 + StMass + Sex + Age +(1|Year)+(0+StMass|Year), data=YearPheno, family=binomial))
summary(glmer(Phi ~ 1 + StMass + Sex + Age + (1|Year) + (0+StMass|Year), data=YearPheno, family=binomial))


#####☻ RHO ######

SelByYearRho <- vector(length = 2015-2005)
SeSelByYearRho <- vector(length = 2015-2005)
CISelByYearRho <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  m0 <- glm(Rho ~ 1 + StMass + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Age=="A",], family=quasipoisson)
  CISelByYearRho[,t-2005]<-as.numeric(confint(m0,parm="StMass"))
  SelByYearRho[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelByYearRho[t-2005] <- sm0$coefficients[2,2]
}
plot(SelByYearRho, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
sd(SelByYearRho)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYearRho[1,],
       y1 = CISelByYearRho[2,], angle = 90,length = 0.1)
m0allRho <- glm(Rho ~ 1 + StMass + Sex , data=YearPheno[YearPheno$Age=="A",], family=quasipoisson)
abline(h=coefficients(m0allRho)[2], lty=2)
sm0allRho <- summary(m0allRho)
lowm0allRho <- coefficients(m0allRho)[2]+1.96*sm0allRho$coefficients[2,2]
highm0allRho <- coefficients(m0allRho)[2]-1.96*sm0allRho$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allRho,lowm0allRho, highm0allRho, highm0allRho),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)

summary(glm(Rho ~ 1 + StMass + Sex  , data=YearPheno[YearPheno$Age=="A",], family=poisson))
summary(glmer(Rho ~ 1  + StMass + Sex  +(1|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson))
summary(glmer(Rho ~ 1  + StMass + Sex  +(1|Year)+(0+StMass|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson))

var(SelByYearRho)+var(SelByYearPhi)+2*cov(SelByYearRho,SelByYearPhi)
var(SelByYear)
var(SelByYearRho+SelByYearPhi)
cor(SelByYearRho,SelByYearPhi)
