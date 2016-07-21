library(MASS)
# YearPheno <- read.table(file = "YearPheno.txt", header=T)
# 
# SelByYear <- vector(length = 2015-2006)
# SeSelByYear <- vector(length = 2015-2006)
# CISelByYear <- matrix(NA,nrow=2,ncol=2015-2005)
# for (t in 2006:2015)
# {
#   m0 <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson)
#   SelByYear[t-2005] <- coefficients(m0)[2]
#   sm0<-summary(m0)
#   SeSelByYear[t-2005] <- sm0$coefficients[2,2]
#   if(t<2015) {CISelByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "StMass"))}else{
#     CISelByYear[,t-2005] <- c(SelByYear[t-2005]-1.96*SeSelByYear[t-2005],SelByYear[t-2005]+1.96*SeSelByYear[t-2005])
#   }
# }
# plot(SelByYear, x=2006:2015, ylim=c(-1.5,0.8), xlab="Year", ylab = "Selection gradient")
# abline(h=0)
# arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYear[1,],
#        y1 = CISelByYear[2,], angle = 90,length = 0.1)
# m0all <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=quasipoisson)
# abline(h=coefficients(m0all)[2], lty=2)
# sm0all <- summary(m0all)
# lowm0all <- coefficients(m0all)[2]+1.96*sm0all$coefficients[2,2]
# highm0all <- coefficients(m0all)[2]-1.96*sm0all$coefficients[2,2]
# polygon(x=c(2005,2016,2016,2005),y=c(lowm0all,lowm0all, highm0all, highm0all),
#         fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)
# 
# 
# 
# 
# summary(glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=poisson))
# 
# mmRRfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1+Mass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)
# mmRIfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)
# mmRnoCorfitness <- glmer(Fitness ~ 1 + StMass + Sex + Age +(1|Year) + (0+StMass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson)
# 
# summary(mmRRfitness)
# summary(mmRnoCorfitness)
# logLik(mmRRfitness)
# logLik(mmRIfitness)
# anova(mmRIfitness,mmRRfitness)
# anova(mmRIfitness,mmRnoCorfitness)
# anova(mmRRfitness,mmRnoCorfitness)
# CImmRnoCorfitness <- confint(mmRnoCorfitness)
# 
# points(x=2006:2015,y=unlist(coefficients(mmRnoCorfitness)$Year["StMass"]))
# 
# 
# 
# prior0 <- list(G=list(G1= list(V=0.00001*diag(2), nu=2)), R=list(V=diag(1),nu=1))
# m0 <- MCMCglmm(fixed = Fitness ~ 1 + Sex + Age + StMass,
#                random = ~ us(1+StMass):Year, prior= prior0, family="poisson",
#                data=YearPheno[!is.na(YearPheno$StMass),], nitt = 110000, burnin = 10000, thin = 100)
# summary(m0)
# autocorr(m0$VCV)
# plot(m0)
# 
# 
# ##### Fitness components ####
# 
# SelByYearPhi <- vector(length = 2015-2005)
# SeSelByYearPhi <- vector(length = 2015-2005)
# CISelByYearPhi <- matrix(NA,nrow=2,ncol=2015-2005)
# 
# for (t in 2006:2014)
# {
#   m0 <- glm(Phi ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=binomial)
#   SelByYearPhi[t-2005] <- coefficients(m0)[2]
#   CISelByYearPhi[,t-2005]<-as.numeric(confint(m0,parm="StMass"))
#   sm0<-summary(m0)
#   SeSelByYearPhi[t-2005] <- sm0$coefficients[2,2]
# }
# SelByYearPhi[10] <- NA
# SeSelByYearPhi[10] <- NA
# plot(SelByYearPhi, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
# abline(h=0)
# arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYearPhi[1,],
#        y1 = CISelByYearPhi[2,], angle = 90,length = 0.1)
# m0allphi <- glm(Phi ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year<2015,], family=binomial)
# abline(h=coefficients(m0allphi)[2], lty=2)
# sm0allphi <- summary(m0allphi)
# lowm0allphi <- coefficients(m0allphi)[2]+1.96*sm0allphi$coefficients[2,2]
# highm0allphi <- coefficients(m0allphi)[2]-1.96*sm0allphi$coefficients[2,2]
# polygon(x=c(2005,2016,2016,2005),y=c(lowm0allphi,lowm0allphi, highm0allphi, highm0allphi),
#         fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )
# 
# summary(glm(Phi ~ 1 + Mass + Sex , data=YearPheno[YearPheno$Age=="A" & YearPheno$Year<max(YearPheno$Year),], family=binomial))
# summary(glmer(Phi ~ 1 + as.factor(Year) + Mass + Sex +(0+Mass|Year), data=YearPheno[YearPheno$Age=="A"& YearPheno$Year<max(YearPheno$Year),], family=binomial))
# 
# mmRIphi <- glmer(Phi ~ 1 + StMass + Sex + Age +(1|Year)+(0+StMass|Year), data=YearPheno, family=binomial)
# mmRnoCorphi <- glmer(Phi ~ 1 + StMass + Sex + Age + (1|Year) + (0+StMass|Year), data=YearPheno, family=binomial)
# summary(mmRnoCorphi)
# CImmRnoCorphi <- confint(mmRnoCorphi)
# 
# ##### RHO ######
# 
# SelByYearRho <- vector(length = 2015-2005)
# SeSelByYearRho <- vector(length = 2015-2005)
# CISelByYearRho <- matrix(NA,nrow=2,ncol=2015-2005)
# for (t in 2006:2015)
# {
#   m0 <- glm(Rho ~ 1 + StMass + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Age=="A",], family=quasipoisson)
#   CISelByYearRho[,t-2005]<-as.numeric(confint(m0,parm="StMass"))
#   SelByYearRho[t-2005] <- coefficients(m0)[2]
#   sm0<-summary(m0)
#   SeSelByYearRho[t-2005] <- sm0$coefficients[2,2]
# }
# plot(SelByYearRho, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
# abline(h=0)
# sd(SelByYearRho)
# arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelByYearRho[1,],
#        y1 = CISelByYearRho[2,], angle = 90,length = 0.1)
# m0allRho <- glm(Rho ~ 1 + StMass + Sex , data=YearPheno[YearPheno$Age=="A",], family=quasipoisson)
# abline(h=coefficients(m0allRho)[2], lty=2)
# sm0allRho <- summary(m0allRho)
# lowm0allRho <- coefficients(m0allRho)[2]+1.96*sm0allRho$coefficients[2,2]
# highm0allRho <- coefficients(m0allRho)[2]-1.96*sm0allRho$coefficients[2,2]
# polygon(x=c(2005,2016,2016,2005),y=c(lowm0allRho,lowm0allRho, highm0allRho, highm0allRho),
#         fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)
# 
# summary(glm(Rho ~ 1 + StMass + Sex  , data=YearPheno[YearPheno$Age=="A",], family=poisson))
# mmRIrho <- glmer(Rho ~ 1  + StMass + Sex  +(1|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson)
# mmRnoCorrho <- glmer(Rho ~ 1  + StMass + Sex  +(1|Year)+(0+StMass|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson)
# 
# var(SelByYearRho)+var(SelByYearPhi)+2*cov(SelByYearRho,SelByYearPhi)
# var(SelByYear)
# var(SelByYearRho+SelByYearPhi)
# cor(SelByYearRho,SelByYearPhi)
# 
# CImmRnoCorrho <- confint(mmRnoCorrho)
# 
# ###########################################
# ##################  QG ####################
# 
# library(MCMCglmm)
# library(pedantics)
# ped <- read.table(file = "ped.txt", header=TRUE)
# 
# YearPheno$animal <- YearPheno$ID
# priorBLUPS0<-list(G=list(G1=list(V=0.1, nu=0.0001),G2=list(V=0.1, nu=0.0001),G3=list(V=0.1, nu=0.0001)),
#                   R=list(V=0.1, nu=0.0001))
# mcmcBLUPS0 <- MCMCglmm(Mass ~Sex+Age,
#                        random=~animal+Mother+Year,
#                        rcov=~units,
#                        prior=priorBLUPS0,
#                        pedigree=ped,data=YearPheno,verbose=TRUE,nitt=12000,burnin=2000,thin=10,pr=TRUE)
# summary(mcmcBLUPS0)
# 
# 
# BV<-mcmcBLUPS0$Sol[,grep(pattern = "animal*",x = colnames(mcmcBLUPS0$Sol))]
# animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
# colnames(BV) <- animalID
# 
# pmBV<-data.frame(animalID,posterior.mode(BV))
# names(pmBV)<-c("ID","pBV")
# mpmBV<-merge(x = pmBV,y = YearPheno,by="ID",all.y=TRUE, all.x = FALSE)
# plot(mpmBV$Year,mpmBV$pBV)
# 
# BVextend <- BV[,mpmBV$ID]# duplicates posterior distribution for ind present in multiple years
# 
# lmBV<-as.mcmc(apply(BVextend,MARGIN = 1,function(x){coef(lm(x~1+mpmBV$Year))[2]}))
# 
# library(mgcv)
# 
# bvplotlist <- list()
# for (i in 1:nrow(BVextend))
# {
#   damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
#   gm0 <- gam(bv~1+s(t),data=damdat)
#   plotgm0 <- plot.gam(gm0,n = 20)
#   bvplotlist[[i]] <- cbind(plotgm0[[1]]$x,plotgm0[[1]]$fit)
# }
# plot(x=0,xlim=c(2006,2015),ylim=c(-2,2),type="n")
# trashidontwantyou<-lapply(bvplotlist, function(x){lines(x[,1],x[,2], col=rgb(0.1,0.1,0.1,alpha = 0.1))})
# 
# bvpairwise <- as.data.frame(matrix(NA, nrow= nrow(BVextend), ncol = 2015-2006))
# names(bvpairwise) <- 2007:2015
# for (i in 1:nrow(BVextend))
# {
#   damdat<- data.frame(bv=BVextend[i,],t=mpmBV$Year)
#   lm(bv~1+t,data=damdat[damdat$t==2006 | damdat$t==2007,])
#   mean(damdat$bv[damdat$t==2007])-mean(damdat$bv[damdat$t==2006])
#   tmeanbv <- tapply(damdat$bv,damdat$t,mean)
#   bvpairwise[i,] <- tmeanbv[-1]-tmeanbv[-10]
# }
# 
# boxplot(bvpairwise)
# abline(h=0)
# 
# #### BIV ANIMAL MODEL ALL YEARS #### (SCARY) # For Mass ####
# 
# nucovM <- 3
# nucovF <- 3
# nucovMF <- 10^8
# spePriornu <- matrix(c(nucovF,nucovMF,nucovMF,nucovMF,nucovMF,nucovM,nucovMF,nucovM,nucovMF),nrow = 3,byrow = T)
# 
# covAM <- 0
# covAMF <- 0
# vAF <- 1
# vAM <- posterior.mode(mcmcBLUPS0$VCV[,"animal"])
# Gcovvar <- matrix(c(vAF,covAMF,covAMF,covAMF,vAM,covAM,covAMF,covAM,vAM),nrow = 3,byrow = T)
# 
# covAM <- 0
# covAMF <- 0
# vAF <- 1
# vIM <- posterior.mode(mcmcBLUPS0$VCV[,"ID"])
# Gcovvar <- matrix(c(vAF,covAMF,covAMF,covAMF,vAM,covAM,covAMF,covAM,vAM),nrow = 3,byrow = T)
# 
# 
# priorAllYears <- list(G=list(G1=list(V=Gcovvar, nu=11, alpha.mu=rep(0,11), alpha.V=diag(11)*1000, fix =2),
#                              G2=list(V=diag(11), nu=11, alpha.mu=rep(0,11), alpha.V=diag(11)*1000, fix =2),
#                              G3=list(V=diag(11), nu=11, alpha.mu=rep(0,11), alpha.V=diag(11)*1000, fix =2),
#                              G4=list(V=diag(1), nu=0.001)),
#                       R=list(V=diag(11), nu=11, fix =2))
# 
# mcmcBivAllYears <- MCMCglmm(cbind(Fitness,M2006,M2007,M2008,M2009,M2010,M2011,M2012,M2013,M2014,M2015) ~ trait-1+Sex+Age+at.level(trait,c(1)):(Sex+Age),
#                             random=~us(trait):animal+us(trait):ID+us(trait):Mother+us(at.level(trait,c(11))):Year,
#                             rcov=~us(trait):units, family=rep("gaussian",11),
#                             prior=priorAllYears,
#                             pedigree=ped,data=YearPheno,verbose=TRUE,nitt=120000,burnin=20000,thin=100)
# summary(mcmcBivAllYears)
# save(mcmcBivAllYears,file = "mcmcBivAllYears1")
# 

############################################################################################################
######################################### Asymptotic mass analysis #########################################
############################################################################################################
library(MASS)
library(lme4)
YearPheno <- read.table(file = "YearPheno.txt", header=T)

SelAByYear <- vector(length = 2015-2006)
SeSelAByYear <- vector(length = 2015-2006)
CISelAByYear <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  m0 <- glm(Fitness ~ 1 + Ast + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=quasipoisson)
  SelAByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelAByYear[t-2005] <- sm0$coefficients[2,2]
  if(t<2015) {CISelAByYear[,t-2005] <- as.numeric(confint(object = m0,parm = "Ast"))}else{
    CISelAByYear[,t-2005] <- c(SelAByYear[t-2005]-1.96*SeSelAByYear[t-2005],SelAByYear[t-2005]+1.96*SeSelAByYear[t-2005])
  }
}
plot(SelAByYear, x=2006:2015, ylim=c(-1.5,0.8), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYear[1,],
       y1 = CISelAByYear[2,], angle = 90,length = 0.1)
m0all <- glm(Fitness ~ 1 + Ast + Sex +Age , data=YearPheno, family=quasipoisson)
abline(h=coefficients(m0all)[2], lty=2)
sm0all <- summary(m0all)
lowm0all <- coefficients(m0all)[2]+1.96*sm0all$coefficients[2,2]
highm0all <- coefficients(m0all)[2]-1.96*sm0all$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0all,lowm0all, highm0all, highm0all),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)


mmARRfitness <- glmer(Fitness ~ 1 + Ast + Sex +  Age +(1+Ast|Year),
                      data=YearPheno, family=poisson, na.action = "na.omit")
mmARIfitness <- glmer(Fitness ~ 1 + Ast + Sex +  Age  +(1|Year), data=YearPheno, family=poisson, na.action = "na.omit")
mmARnoCorfitness <- glmer(Fitness ~ 1 + Ast + Sex +  Age  +(1|Year) + (0+Ast|Year),
                          data=YearPheno, family=poisson, na.action = "na.omit")

summary(mmARRfitness)
anova(mmARRfitness,mmARnoCorfitness)
smmARnoCorfitness <- summary(mmARnoCorfitness)
fitnessAanova <- anova(mmARIfitness,mmARnoCorfitness)

CImmARnoCorfitness <- confint(mmARnoCorfitness)
                         
##### Fitness components ####
### AD + Juv
SelAByYearPhi <- vector(length = 2015-2005)
SeSelAByYearPhi <- vector(length = 2015-2005)
CISelAByYearPhi <- matrix(NA,nrow=2,ncol=2015-2005)

for (t in 2006:2014)
{
  m0 <- glm(Phi ~ 1 + Ast + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=binomial)
  SelAByYearPhi[t-2005] <- coefficients(m0)[2]
  CISelAByYearPhi[,t-2005]<-as.numeric(confint(m0,parm="Ast"))
  sm0<-summary(m0)
  SeSelAByYearPhi[t-2005] <- sm0$coefficients[2,2]
}
SelAByYearPhi[10] <- NA
SeSelAByYearPhi[10] <- NA
plot(SelAByYearPhi, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearPhi[1,],
       y1 = CISelAByYearPhi[2,], angle = 90,length = 0.1)
m0allphi <- glm(Phi ~ 1 + Ast + Sex +Age , data=YearPheno[YearPheno$Year<2015,], family=binomial)
abline(h=coefficients(m0allphi)[2], lty=2)
sm0allphi <- summary(m0allphi)
lowm0allphi <- coefficients(m0allphi)[2]+1.96*sm0allphi$coefficients[2,2]
highm0allphi <- coefficients(m0allphi)[2]-1.96*sm0allphi$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allphi,lowm0allphi, highm0allphi, highm0allphi),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

summary(glm(Phi ~ 1 + A + Sex , data=YearPheno[YearPheno$Age=="A" & YearPheno$Year<max(YearPheno$Year),], family=binomial))
summary(glmer(Phi ~ 1 + as.factor(Year) + A + Sex +(0+A|Year), data=YearPheno[YearPheno$Age=="A"& YearPheno$Year<max(YearPheno$Year),], family=binomial))

mmRIphi <- glmer(Phi ~ 1 + Ast + Sex + Age +(1|Year), data=YearPheno, family=binomial)
mmRnoCorphi <- glmer(Phi ~ 1 + Ast + Sex + Age + (1|Year) + (0+Ast|Year), data=YearPheno, family=binomial)
  smmRnoCorphi <- summary(mmRnoCorphi)
CImmRnoCorphi <- confint(mmRnoCorphi)

PhiAanova <- anova(mmRIphi,mmRnoCorphi)

#### Ad only
SelAByYearPhiAd <- vector(length = 2015-2005)
SeSelAByYearPhiAd <- vector(length = 2015-2005)
CISelAByYearPhiAd <- matrix(NA,nrow=2,ncol=2015-2005)

for (t in 2006:2014)
{
  m0 <- glm(Phi ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Age=="A",], family=binomial)
  SelAByYearPhiAd[t-2005] <- coefficients(m0)[2]
  CISelAByYearPhiAd[,t-2005]<-as.numeric(confint(m0,parm="Ast"))
  sm0<-summary(m0)
  SeSelAByYearPhiAd[t-2005] <- sm0$coefficients[2,2]
}
SelAByYearPhiAd[10] <- NA
SeSelAByYearPhiAd[10] <- NA
plot(SelAByYearPhiAd, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearPhiAd[1,],
       y1 = CISelAByYearPhiAd[2,], angle = 90,length = 0.1)
m0allphiAd <- glm(Phi ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Year<2015 & YearPheno$Age=="A",], family=binomial)
abline(h=coefficients(m0allphiAd)[2], lty=2)
sm0allphiAd <- summary(m0allphiAd)
lowm0allphiAd <- coefficients(m0allphiAd)[2]+1.96*sm0allphiAd$coefficients[2,2]
highm0allphiAd <- coefficients(m0allphiAd)[2]-1.96*sm0allphiAd$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allphiAd,lowm0allphiAd, highm0allphiAd, highm0allphiAd),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

mmRIphiAd <- glmer(Phi ~ 1 + Ast + Sex +(1|Year), data=YearPheno[YearPheno$Age=="A",], family=binomial)
mmRnoCorphiAd <- glmer(Phi ~ 1 + Ast + Sex + (1|Year) + (0+Ast|Year), data=YearPheno[YearPheno$Age=="A",], family=binomial)
smmRnoCorphiAd <- summary(mmRnoCorphiAd)
CImmRnoCorphiAd <- confint(mmRnoCorphiAd)

anova(mmRIphiAd,mmRnoCorphiAd)

### Juvenile only
SelAByYearPhiJuv <- vector(length = 2015-2005)
SeSelAByYearPhiJuv <- vector(length = 2015-2005)
CISelAByYearPhiJuv <- matrix(NA,nrow=2,ncol=2015-2005)

for (t in 2006:2014)
{
  m0 <- glm(Phi ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Age=="J",], family=binomial)
  SelAByYearPhiJuv[t-2005] <- coefficients(m0)[2]
  CISelAByYearPhiJuv[,t-2005]<-as.numeric(confint(m0,parm="Ast"))
  sm0<-summary(m0)
  SeSelAByYearPhiJuv[t-2005] <- sm0$coefficients[2,2]
}
SelAByYearPhiJuv[10] <- NA
SeSelAByYearPhiJuv[10] <- NA
plot(SelAByYearPhiJuv, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearPhiJuv[1,],
       y1 = CISelAByYearPhiJuv[2,], angle = 90,length = 0.1)
m0allphiJuv <- glm(Phi ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Year<2015 & YearPheno$Age=="J",], family=binomial)
abline(h=coefficients(m0allphiJuv)[2], lty=2)
sm0allphiJuv <- summary(m0allphiJuv)
lowm0allphiJuv <- coefficients(m0allphiJuv)[2]+1.96*sm0allphiJuv$coefficients[2,2]
highm0allphiJuv <- coefficients(m0allphiJuv)[2]-1.96*sm0allphiJuv$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allphiJuv,lowm0allphiJuv, highm0allphiJuv, highm0allphiJuv),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2 )

mmRIphiJuv <- glmer(Phi ~ 1 + Ast + Sex +(1|Year), data=YearPheno[YearPheno$Age=="J",], family=binomial)
mmRnoCorphiJuv <- glmer(Phi ~ 1 + Ast + Sex + (1|Year) + (0+Ast|Year), data=YearPheno[YearPheno$Age=="J",], family=binomial)
smmRnoCorphiJuv <- summary(mmRnoCorphiJuv)
CImmRnoCorphiJuv <- confint(mmRnoCorphiJuv)

anova(mmRIphiJuv,mmRnoCorphiJuv)


#####  RHO ######

SelAByYearRho <- vector(length = 2015-2005)
SeSelAByYearRho <- vector(length = 2015-2005)
CISelAByYearRho <- matrix(NA,nrow=2,ncol=2015-2005)
for (t in 2006:2015)
{
  m0 <- glm(Rho ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Year==t & YearPheno$Age=="A",], family=quasipoisson)
  CISelAByYearRho[,t-2005]<-as.numeric(confint(m0,parm="Ast"))
  SelAByYearRho[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelAByYearRho[t-2005] <- sm0$coefficients[2,2]
}
plot(SelAByYearRho, x=2006:2015, ylim=c(-2,2), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = CISelAByYearRho[1,],
       y1 = CISelAByYearRho[2,], angle = 90,length = 0.1)
m0allRho <- glm(Rho ~ 1 + Ast + Sex , data=YearPheno[YearPheno$Age=="A",], family=quasipoisson)
abline(h=coefficients(m0allRho)[2], lty=2)
sm0allRho <- summary(m0allRho)
lowm0allRho <- coefficients(m0allRho)[2]+1.96*sm0allRho$coefficients[2,2]
highm0allRho <- coefficients(m0allRho)[2]-1.96*sm0allRho$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0allRho,lowm0allRho, highm0allRho, highm0allRho),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2)

summary(glm(Rho ~ 1 + Ast + Sex  , data=YearPheno[YearPheno$Age=="A",], family=poisson))
mmRIrho <- glmer(Rho ~ 1  + Ast + Sex  +(1|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson)
mmRnoCorrho <- glmer(Rho ~ 1  + Ast + Sex  +(1|Year)+(0+Ast|Year), data=YearPheno[YearPheno$Age=="A",], family=poisson)

var(SelAByYear)
var(SelAByYearRho+SelAByYearPhi, na.rm = T)
cor.test(SelAByYearRho,SelAByYearPhi)
smmRnoCorrho <- summary(mmRnoCorrho)
RhoAanova <- anova(mmRIrho,mmRnoCorrho)

CImmRnoCorrho <- confint(mmRnoCorrho)


#### Dynamics of phenotype! ####
plot(tapply(X = YearPheno$A, INDEX = YearPheno$Year, function(x){mean(x,na.rm=TRUE)}))
plot(tapply(X = YearPheno$A, INDEX = YearPheno$Year, function(x){sd(x,na.rm=TRUE)}))
