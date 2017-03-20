library(lme4)
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)

realindnumb <- tapply(YearPheno$ID, YearPheno$Year, function(x){length(unique(x))})

VA <- 1
nind <- 1000

indbv <- rnorm(n = nind, mean = 0, sd = sqrt(VA))

indLatentfitness <- ifelse(indbv>0,1,0)+abs(rnorm(n = nind, mean = 0, sd=0.1))
indLatentfitness <- indLatentfitness/mean(indLatentfitness)

plot(indLatentfitness, x=indbv)
0.99*cov(indLatentfitness, indbv)
drift <- indbv+rnorm(n = nind, 0, VA/2)

mean(drift*indLatentfitness)-mean(indbv)
0.99*cov(indLatentfitness, indbv)

cov(indLatentfitness, indbv)

mean(indbv*indLatentfitness)-mean(indbv)


nbsimuls <- 50
varselrho <- c(0,0.01,0.025,0.05,0.1)
varintrho <- 0.3

rhopv <- list()
for (vsr in 1:length(varselrho))
{
  pvalues <- vector(length = nbsimuls)
  
  for (j in 1:nbsimuls)
  {
    nyears <- 10
    nindyear <- realindnumb
    inddata <- c(NA, NA, NA)
    for (i in 1:nyears)
    {
      phen <- rnorm(nindyear[i],mean = 0,sd = 1)
      sel <- rnorm(1,0,sqrt(varselrho[vsr]))
      shift <- rnorm(1,0,sqrt(varintrho))
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

plot(x=varselrho, y=unlist(lapply(rhopv, function(x){mean(x<0.05)})))



nbsimuls <- 100
varselphi <- c(0,0.05,0.1,0.25,0.5,1)
varintphi <- 0.7
phivar <-list()
phipv <- list()
for (vsr in 1:length(varselphi))
{
  pvalues <- vector(length = nbsimuls)
  realizedvar <- vector(length = nbsimuls)
  for (j in 1:nbsimuls)
  {
    nyears <- 10
    nindyear <- realindnumb
    inddata <- c(NA, NA, NA)
    yearrealizedbeta <- vector(length = nyears)
    for (i in 1:nyears)
    {
      phen <- rnorm(nindyear[i],mean = 0,sd = 1)
      sel <- rnorm(1,0,sqrt(varselphi[vsr]))
      shift <- rnorm(1,0,sqrt(varintphi))
      phi <- rbinom(n = nindyear[i],size = 1,prob = 1/(1+exp(-(-1+shift+sel*phen))))
      year <- rep(x = i,times=nindyear[i])
      inddata<- rbind(inddata, cbind(phen, phi, year ))
      yearrealizedbeta[i] <- cov(phen, phi)/var(phen) 
    }
    realizedvar[j] <- var(yearrealizedbeta)
    inddata <- inddata[!is.na(inddata[,1]),]
    inddata <- as.data.frame(inddata)
    library(lme4)
    
    mfull <- glmer(formula = phi ~ phen + (1+phen|year), data=inddata, family=poisson)
    mnull <- glmer(formula = phi ~ phen + (1|year), data=inddata, family=poisson)
    amn<-anova(mfull, mnull)
    pvalues[j] <- amn$`Pr(>Chisq)`[2]
  }
  phipv[[vsr]] <- pvalues
  phivar[[vsr]] <- realizedvar
}

plot(x=sqrt(varselphi)/0.4, y=unlist(lapply(phipv, function(x){mean(x<0.075)})))#mixture of 2 and 1 degree of freedom
abline(v=1, lty=2)
abline(v=2, lty=2)
plot(x=varselphi, y=unlist(lapply(phivar, function(x){mean(x)})))
plot(x=varselphi, y=unlist(lapply(phivar, function(x){max(x)})))
