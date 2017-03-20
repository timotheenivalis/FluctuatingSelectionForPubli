library(lme4)
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)


#### Global options ####
nbsimuls <- 5 # number of replicates for every scenario
SD_Beta_Levels <- c(0, 0.1, 0.5, 1, 2, 4) #levels of standard deviation over median value for selection
realindnumb <- tapply(YearPheno$ID, YearPheno$Year, function(x){length(unique(x))}) # the sample size per-year
realindnumbRho <- tapply(YearPheno$ID[YearPheno$Phi==1], YearPheno$Year[YearPheno$Phi==1], function(x){length(unique(x))}) # the sample size per-year for reproduction


#### Annualized fitness ####
ZmmARRfitness <- glmer(FitnessZ ~ 1  + BMIst + Sex  +(1+BMIst|Year), data=YearPheno, family=poisson)
SelCoeff_fitness <- getME(ZmmARRfitness, "fixef")["BMIst"]
SDint_fitness <- getME(ZmmARRfitness, "theta")["Year.(Intercept)"]

SDSel_fitness <- SD_Beta_Levels*abs(SelCoeff_fitness)

fitnesspv <- list()
for (vsr in 1:length(SDSel_fitness))
{
  pvalues <- vector(length = nbsimuls)
  
  for (j in 1:nbsimuls)
  {
    nyears <- length(realindnumb)
    nindyear <- realindnumb
    inddata <- c(NA, NA, NA)
    for (i in 1:nyears)
    {
      phen <- rnorm(nindyear[i],mean = 0,sd = 1)
      sel <- rnorm(n = 1, 0, SDSel_fitness[vsr]) # fluctuating selection
      shift <- rnorm(n = 1, 0, SDint_fitness) # uncorrelated intercept variation
      rho <- rpois(n = nindyear[i],lambda = exp(log(2)+shift+sel*phen))
      year <- rep(x = i,times=nindyear[i])
      inddata<- rbind(inddata, cbind(phen, rho, year ))
    }
    inddata <- inddata[!is.na(inddata[,1]),]
    inddata <- as.data.frame(inddata)
    library(lme4)
    
    mfull <- glmer(formula = rho ~ phen + (1+phen|year), data=inddata, family=poisson)
    mnull <- glmer(formula = rho ~ phen + (1|year), data=inddata, family=poisson)
    amn<-anova(mfull, mnull)
    pvalues[j] <- amn$`Pr(>Chisq)`[2]
  }
  fitnesspv[[vsr]] <- pvalues
}

plot(x=SD_Beta_Levels, y=unlist(lapply(fitnesspv, function(x){mean(x<0.075)}))) #1.5 df

#### Annual Reproductive Success ####
ZmmRrho <- glmer(RhoZ ~ 1  + BMIst + Sex  +(1|Year) + (0+BMIst|Year), data=YearPheno[YearPheno$Phi==1,], family=poisson)
SelCoeff_rho <- getME(ZmmRrho, "fixef")["BMIst"]
SDint_rho <- getME(ZmmRrho, "theta")["Year.(Intercept)"]

SDSel_rho <- SD_Beta_Levels*abs(SelCoeff_rho)

rhopv <- list()
for (vsr in 1:length(SDSel_rho))
{
  pvalues <- vector(length = nbsimuls)
  
  for (j in 1:nbsimuls)
  {
    nyears <- length(realindnumbRho)
    nindyear <- realindnumb
    inddata <- c(NA, NA, NA)
    for (i in 1:nyears)
    {
      phen <- rnorm(nindyear[i],mean = 0,sd = 1)
      sel <- rnorm(n = 1, 0, SDSel_rho[vsr]) # fluctuating selection
      shift <- rnorm(n = 1, 0, SDint_rho) # uncorrelated intercept variation
      rho <- rpois(n = nindyear[i],lambda = exp(log(2)+shift+sel*phen))
      year <- rep(x = i,times=nindyear[i])
      inddata<- rbind(inddata, cbind(phen, rho, year ))
    }
    inddata <- inddata[!is.na(inddata[,1]),]
    inddata <- as.data.frame(inddata)
    library(lme4)
    
    mfull <- glmer(formula = rho ~ phen + (1+phen|year), data=inddata, family=poisson)
    mnull <- glmer(formula = rho ~ phen + (1|year), data=inddata, family=poisson)
    amn<-anova(mfull, mnull)
    pvalues[j] <- amn$`Pr(>Chisq)`[2]
  }
  rhopv[[vsr]] <- pvalues
}

plot(x=SD_Beta_Levels, y=unlist(lapply(rhopv, function(x){mean(x<0.075)}))) #1.5 df

#### Survival ####
mmRnoCorphi <- glmer(Phi ~ 1 + BMIst + Sex * Age + (1|Year) + (0+BMIst|Year), data=YearPheno, family=binomial)
SelCoeff_phi <- getME(mmRnoCorphi, "fixef")["BMIst"]
SDint_phi <- getME(mmRnoCorphi, "theta")["Year.(Intercept)"]

SDSel_phi<- SD_Beta_Levels*abs(SelCoeff_phi)

phipv <- list()
for (vsr in 1:length(SDSel_phi))
{
  pvalues <- vector(length = nbsimuls)
  for (j in 1:nbsimuls)
  {
    nyears <- 10
    nindyear <- realindnumb
    inddata <- c(NA, NA, NA)
    yearrealizedbeta <- vector(length = nyears)
    for (i in 1:nyears)
    {
      phen <- rnorm(nindyear[i],mean = 0,sd = 1)
      sel <- rnorm(1,0,SDSel_phi[vsr])
      shift <- rnorm(1,0,SDint_phi)
      phi <- rbinom(n = nindyear[i],size = 1,prob = 1/(1+exp(-(-1+shift+sel*phen))))
      year <- rep(x = i,times=nindyear[i])
      inddata<- rbind(inddata, cbind(phen, phi, year ))
      yearrealizedbeta[i] <- cov(phen, phi)/var(phen) 
    }
    inddata <- inddata[!is.na(inddata[,1]),]
    inddata <- as.data.frame(inddata)
    
    mfull <- glmer(formula = phi ~ phen + (1+phen|year), data=inddata, family=poisson)
    mnull <- glmer(formula = phi ~ phen + (1|year), data=inddata, family=poisson)
    amn<-anova(mfull, mnull)
    pvalues[j] <- amn$`Pr(>Chisq)`[2]
  }
  phipv[[vsr]] <- pvalues
}

plot(x=SD_Beta_Levels, y=unlist(lapply(fitnesspv, function(x){mean(x<0.075)}))) #1.5 df