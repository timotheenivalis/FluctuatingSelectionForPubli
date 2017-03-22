library(lme4)
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)
head(YearPheno)

repetBMI_A <- lmer(formula =BMI ~ 1 + RelativeJulian + I(RelativeJulian^2) + (1|ID) + (1|Year), data=YearPheno[YearPheno$Age=="A",])
summary(repetBMI_A)
576.11/(576.11+508.32)#R2 in adults

repetBMI <- lmer(formula =BMI ~ 1 + Age * RelativeJulian+ Age*I(RelativeJulian^2) + (1|ID) + (1|Year), data=YearPheno)
summary(repetBMI)
434.78/(434.78+789.91)#R2 at all ages

repetBMIst <- lmer(formula =BMIst ~ 1 + Age + (1|ID) + (1|Year), data=YearPheno)
summary(repetBMIst)
0.11986/(0.11986+0.23844)

######### Selection #####

SelAByYear <- vector(length = 2016-2006)
SeSelAByYear <- vector(length = 2016-2006)
CISelAByYear <- matrix(NA,nrow=2,ncol=2016-2005)
for (t in 2006:2016)
{
  m0 <- glm(Fitness ~ 1 + BMIst + Sex * Age + Age*RJst, data=YearPheno[YearPheno$Year==t,], family=quasipoisson, na.action = "na.omit")
  SelAByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelAByYear[t-2005] <- sm0$coefficients[2,2]
  #if(t!=2015) {
    CISelAByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "BMIst"))
  #  }else{
  #  CISelAByYear[,t-2005] <- c(SelAByYear[t-2005]-1.96*SeSelAByYear[t-2005],SelAByYear[t-2005]+1.96*SeSelAByYear[t-2005])
  #}
}
plot(SelAByYear, x=2006:2016, ylim=c(-1.5,0.8), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2016,x1 = 2006:2016,code = 3, y0 = CISelAByYear[1,],
       y1 = CISelAByYear[2,], angle = 90,length = 0.1)
m0all <- glm(Fitness ~ 1 + BMIst + Sex +Age , data=YearPheno, family=quasipoisson)
abline(h=coefficients(m0all)[2], lty=2)
sm0all <- summary(m0all)
lowm0all <- coefficients(m0all)[2]+1.96*sm0all$coefficients[2,2]
highm0all <- coefficients(m0all)[2]-1.96*sm0all$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(lowm0all,lowm0all, highm0all, highm0all),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)


mmARRfitness <- glmer(Fitness ~ 1 + BMIst + Sex *  Age  + Age*RJst +(1+BMIst|Year),
                      data=YearPheno, family=poisson, na.action = "na.omit")
mmARIfitness <- glmer(Fitness ~ 1 + BMIst + Sex * Age  + Age*RJst +(1|Year), data=YearPheno, family=poisson, na.action = "na.omit")
mmARnoCorfitness <- glmer(Fitness ~ 1 + BMIst + Sex *  Age  + Age*RJst +(1|Year) + (0+BMIst|Year),
                          data=YearPheno, family=poisson, na.action = "na.omit")

summary(mmARRfitness)
anova(mmARRfitness,mmARnoCorfitness)
smmARnoCorfitness <- summary(mmARnoCorfitness)
fitnessAanova <- anova(mmARIfitness,mmARnoCorfitness)

CImmARnoCorfitness <- confint(mmARnoCorfitness)


##### Fitness components ####
### AD + Juv
SelAByYearPhi <- vector(length = 2016-2005)
SeSelAByYearPhi <- vector(length = 2016-2005)
CISelAByYearPhi <- matrix(NA,nrow=2,ncol=2016-2005)

for (t in 2006:2015)
{
  m0 <- glm(Phi ~ 1 + BMIst + Sex *Age  + Age*RJst , data=YearPheno[YearPheno$Year==t,], family=binomial, na.action = "na.omit")
  SelAByYearPhi[t-2005] <- coefficients(m0)[2]
  CISelAByYearPhi[,t-2005]<-as.numeric(confint(m0,parm="BMIst"))
  sm0<-summary(m0)
  SeSelAByYearPhi[t-2005] <- sm0$coefficients[2,2]
}
SelAByYearPhi[11] <- NA
SeSelAByYearPhi[11] <- NA
plot(SelAByYearPhi, x=2006:2016, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearPhi[1,],
       y1 = CISelAByYearPhi[2,], angle = 90,length = 0.1)
m0allphi <- glm(Phi ~ 1 + BMIst + Sex *Age + Age*RJst, data=YearPheno[YearPheno$Year<2015,], family=binomial, na.action = "na.omit")
abline(h=coefficients(m0allphi)[2], lty=2)
sm0allphi <- summary(m0allphi)
lowm0allphi <- coefficients(m0allphi)[2]+1.96*sm0allphi$coefficients[2,2]
highm0allphi <- coefficients(m0allphi)[2]-1.96*sm0allphi$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(lowm0allphi,lowm0allphi, highm0allphi, highm0allphi),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

mmRIphi <- glmer(Phi ~ 1 + BMIst + Sex * Age + Age*RJst +(1|Year), data=YearPheno, family=binomial, na.action = "na.omit")
mmRnoCorphi <- glmer(Phi ~ 1 + BMIst + Sex * Age + Age*RJst+  (1|Year) + (0+BMIst|Year), data=YearPheno, family=binomial, na.action = "na.omit")
smmRnoCorphi <- summary(mmRnoCorphi)
CImmRnoCorphi <- confint(mmRnoCorphi)

PhiAanova <- anova(mmRIphi,mmRnoCorphi)

#####  RHO ######

SelAByYearRho <- vector(length = 2016-2005)
SeSelAByYearRho <- vector(length = 2016-2005)
CISelAByYearRho <- matrix(NA,nrow=2,ncol=2016-2005)
for (t in 2006:2016)
{
  m0 <- glm(Rho ~ 1 + BMIst + Sex + RJst, data=YearPheno[YearPheno$Year==t & YearPheno$Age=="A",], family=quasipoisson, na.action = "na.omit")
  CISelAByYearRho[,t-2005]<-as.numeric(confint(m0,parm="BMIst"))
  SelAByYearRho[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelAByYearRho[t-2005] <- sm0$coefficients[2,2]
}
plot(SelAByYearRho, x=2006:2016, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2016,x1 = 2006:2016,code = 3, y0 = CISelAByYearRho[1,],
       y1 = CISelAByYearRho[2,], angle = 90,length = 0.1)
m0allRho <- glm(Rho ~ 1 + BMIst + Sex +RJst, data=YearPheno[YearPheno$Age=="A",], family=quasipoisson, na.action = "na.omit")
abline(h=coefficients(m0allRho)[2], lty=2)
sm0allRho <- summary(m0allRho)
lowm0allRho <- coefficients(m0allRho)[2]+1.96*sm0allRho$coefficients[2,2]
highm0allRho <- coefficients(m0allRho)[2]-1.96*sm0allRho$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(lowm0allRho,lowm0allRho, highm0allRho, highm0allRho),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)

summary(glm(Rho ~ 1 + BMIst + Sex  +RJst, data=YearPheno[YearPheno$Age=="A",], family=poisson, na.action = "na.omit"))
mmRIrho <- glmer(Rho ~ 1  + BMIst + Sex  +(1|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson, na.action = "na.omit")
mmRnoCorrho <- glmer(Rho ~ 1  + BMIst + Sex  +(1|Year)+(0+BMIst|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson, na.action = "na.omit")

var(SelAByYear)
var(SelAByYearRho+SelAByYearPhi, na.rm = T)
cor.test(SelAByYearRho,SelAByYearPhi)
smmRnoCorrho <- summary(mmRnoCorrho)
RhoAanova <- anova(mmRIrho,mmRnoCorrho)

CImmRnoCorrho <- confint(mmRnoCorrho)



#### Dynamics of phenotype! ####
plot(tapply(X = YearPheno$BMI, INDEX = YearPheno$Year, function(x){mean(x,na.rm=TRUE)}))
plot(tapply(X = YearPheno$BMI, INDEX = YearPheno$Year, function(x){sd(x,na.rm=TRUE)}))


# save.image("~/thesis/Mass/FluctuatingSelectionForPubli/EnvSelBMI.RData")


################################################################################################%
################################### ZYGOTE SELECTION ############################################
################################################################################################%
library(lme4)
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)
head(YearPheno)
######### Selection #####

ZSelAByYear <- vector(length = 2015-2006)
ZSeSelAByYear <- vector(length = 2015-2006)
ZCISelAByYear <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  Zm0 <- glm(FitnessZ ~ 1 + BMIst + Sex * Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson)
  ZSelAByYear[t-2005] <- coefficients(Zm0)[2]
  Zsm0<-summary(Zm0)
  ZSeSelAByYear[t-2005] <- Zsm0$coefficients[2,2]
  if(t<2015) {ZCISelAByYear[,t-2005] <- as.numeric(confint(object = Zm0,parm = "BMIst"))}else{
    ZCISelAByYear[,t-2005] <- c(ZSelAByYear[t-2005]-1.96*ZSeSelAByYear[t-2005],ZSelAByYear[t-2005]+1.96*ZSeSelAByYear[t-2005])
  }
}
plot(ZSelAByYear, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = ZCISelAByYear[1,],
       y1 = ZCISelAByYear[2,], angle = 90,length = 0.1)
Zm0all <- glm(FitnessZ ~ 1 + BMIst + Sex +Age , data=YearPheno, family=quasipoisson)
abline(h=coefficients(Zm0all)[2], lty=2)
Zsm0all <- summary(Zm0all)
Zlowm0all <- coefficients(Zm0all)[2]+1.96*Zsm0all$coefficients[2,2]
Zhighm0all <- coefficients(Zm0all)[2]-1.96*Zsm0all$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(Zlowm0all,Zlowm0all, Zhighm0all, Zhighm0all),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)


ZmmARRfitness <- glmer(FitnessZ ~ 1 + BMIst + Sex *  Age +(1+BMIst|Year),
                      data=YearPheno, family=poisson, na.action = "na.omit")
ZmmARIfitness <- glmer(FitnessZ ~ 1 + BMIst + Sex * Age  +(1|Year), data=YearPheno, family=poisson, na.action = "na.omit")
ZmmARnoCorfitness <- glmer(FitnessZ ~ 1 + BMIst + Sex *  Age  +(1|Year) + (0+BMIst|Year),
                          data=YearPheno, family=poisson, na.action = "na.omit")

summary(ZmmARRfitness)
anova(ZmmARRfitness,ZmmARnoCorfitness)
anova(ZmmARRfitness,ZmmARIfitness)

ZsmmARnoCorfitness <- summary(ZmmARnoCorfitness)
ZfitnessAanova <- anova(ZmmARIfitness,ZmmARRfitness)
ZfitnessAanova$`Pr(>Chisq)`*0.75

ZCImmARnoCorfitness <- confint(ZmmARnoCorfitness)
ZCImmARRfitness <- confint(ZmmARRfitness)

##### Fitness components ####
#### AD + Juv survival #### SAME AS BEFORE
SelAByYearPhi <- vector(length = 2015-2005)
SeSelAByYearPhi <- vector(length = 2015-2005)
CISelAByYearPhi <- matrix(NA,nrow=2,ncol=2015-2005)

for (t in 2006:2015)
{
  m0 <- glm(Phi ~ 1 + BMIst + Sex *Age , data=YearPheno[YearPheno$Year==t,], family=binomial)
  SelAByYearPhi[t-2005] <- coefficients(m0)[2]
  CISelAByYearPhi[,t-2005]<-as.numeric(confint(m0,parm="BMIst"))
  sm0<-summary(m0)
  SeSelAByYearPhi[t-2005] <- sm0$coefficients[2,2]
}

var(SelAByYearPhi, na.rm = TRUE)

plot(SelAByYearPhi, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearPhi[1,],
       y1 = CISelAByYearPhi[2,], angle = 90,length = 0.1)
m0allphi <- glm(Phi ~ 1 + BMIst + Sex *Age , data=YearPheno[YearPheno$Year<2015,], family=binomial)
abline(h=coefficients(m0allphi)[2], lty=2)
sm0allphi <- summary(m0allphi)
lowm0allphi <- coefficients(m0allphi)[2]+1.96*sm0allphi$coefficients[2,2]
highm0allphi <- coefficients(m0allphi)[2]-1.96*sm0allphi$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(lowm0allphi,lowm0allphi, highm0allphi, highm0allphi),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

mmRIphi <- glmer(Phi ~ 1 + BMIst + Sex * Age +(1|Year), data=YearPheno, family=binomial)
mmRnoCorphi <- glmer(Phi ~ 1 + BMIst + Sex * Age + (1|Year) + (0+BMIst|Year), data=YearPheno, family=binomial)
smmRnoCorphi <- summary(mmRnoCorphi)
CImmRnoCorphi <- confint(mmRnoCorphi)

PhiAanova <- anova(mmRIphi,mmRnoCorphi)

#####  RHO ######

ZSelAByYearRho <- vector(length = 2015-2005)
ZSeSelAByYearRho <- vector(length = 2015-2005)
ZCISelAByYearRho <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{

  m0 <- glm(RhoZ ~ 1 + BMIst + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Phi==1,], family=quasipoisson)
  ZCISelAByYearRho[,t-2005]<-as.numeric(confint(m0,parm="BMIst"))
  ZSelAByYearRho[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  ZSeSelAByYearRho[t-2005] <- sm0$coefficients[2,2]
}
plot(ZSelAByYearRho, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = ZCISelAByYearRho[1,],
       y1 = ZCISelAByYearRho[2,], angle = 90,length = 0.1)
Zm0allRho <- glm(RhoZ ~ 1 + BMIst + Sex , data=YearPheno[YearPheno$Phi==1,], family=quasipoisson)
abline(h=coefficients(Zm0allRho)[2], lty=2)
Zsm0allRho <- summary(Zm0allRho)
Zlowm0allRho <- coefficients(Zm0allRho)[2]+1.96*Zsm0allRho$coefficients[2,2]
Zhighm0allRho <- coefficients(Zm0allRho)[2]-1.96*Zsm0allRho$coefficients[2,2]
polygon(x=c(2005,2017,2017,2005),y=c(Zlowm0allRho,Zlowm0allRho, Zhighm0allRho, Zhighm0allRho),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)

ZmmRIrho <- glmer(RhoZ ~ 1  + BMIst + Sex  +(1|Year), data=YearPheno[ YearPheno$Phi==1,], family=poisson)
ZmmRnoCorrho <- glmer(RhoZ ~ 1  + BMIst + Sex  +(1|Year)+(0+BMIst|Year), data=YearPheno[ YearPheno$Phi==1,], family=poisson)
ZmmRrho <- glmer(RhoZ ~ 1  + BMIst + Sex  +(1+BMIst|Year), data=YearPheno[YearPheno$Phi==1,], family=poisson)

var(ZSelAByYear)
var(ZSelAByYearRho+SelAByYearPhi, na.rm = T)
cor.test(ZSelAByYearRho,SelAByYearPhi)
ZsmmRnoCorrho <- summary(ZmmRnoCorrho)
ZRhoAanovaNoCore <- anova(ZmmRIrho,ZmmRnoCorrho)
ZRhoAanova <- anova(ZmmRIrho,ZmmRrho)

ZCImmRnoCorrho <- confint(ZmmRnoCorrho)


