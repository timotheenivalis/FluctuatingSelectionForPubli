setwd("/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/IBSimulations/")

#### PRICE EQUATION #### 
evolsel <- read.table("CheckPriceEq/EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta) # weird pattern due to shrinkage of G variation? 
plot(evolsel$DeltaBV, x=evolsel$Var)

#### BREEDER'S EQUATION WITH H²=1####
evolsel <- read.table("CheckBreederEd_H2_1//EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(corv, ylim=c(0,1))
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta)
plot(evolsel$DeltaBV, x=evolsel$Var)

#### BREEDER'S EQUATION WITH H²=0.17####
evolsel <- read.table("CheckBreederEd_H2_017/EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(corv, ylim=c(0,1))
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta)
plot(evolsel$DeltaBV, x=evolsel$Var)

