library(MCMCglmm)

setwd("/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/IBSimulations/")

FullDataSet <- read.table(file = "tryh2eq0.5/Captures.txt", header = TRUE)
FullEvolSel <- read.table(file = "tryh2eq0.5/EvolSel.txt", header = TRUE)


#### Load data ####

evolsel <- FullEvolSel#[FullEvolSel$RUN==1,]
captures <- FullDataSet#[FullDataSet$RUN==1,]

captures$fitness <- captures$Rho+ 2*captures$Phi 
captures$pheno <- captures$BVZ + captures$EVZ


#### Check if we measure the right selection/evolution ####
selyear<-vector(length = max(captures$Year)-1)
postsyear <- vector(length = max(captures$Year)-1)
evolyear<-vector(length = max(captures$Year)-1)

for (i in 1:max(captures$Year))
{
  ny<-length(captures$pheno[captures$Year==i-1])
  selyear[i] <- ((ny-1)/ny)*cov(captures$pheno[captures$Year==i-1], captures$fitness[captures$Year==i-1]/mean(captures$fitness[captures$Year==i-1]))
  evolyear[i] <- mean(captures$BVZ[captures$Year==i])-mean(captures$BVZ[captures$Year==i-1])
}

mean(evolyear)/mean(selyear)
lm(evolyear ~ selyear)

cor(selyear, evolsel$Sel)
plot(selyear, evolyear)

plot(evolsel$DeltaBV, evolyear)
cor(evolsel$DeltaBV, evolyear)
plot(evolsel$DeltaBV, evolsel$Sel)
plot(evolsel$Sel, selyear)

#### Is selection positive of negative ? ####

evolsel$Sel>0

