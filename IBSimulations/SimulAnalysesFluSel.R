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

summary(lm(DeltaBV ~ Sel, data=evolsel))
#### BREEDER'S EQUATION WITH H²=0.17####
evolsel <- read.table("CheckBreederEd_H2_017/EvolSel.txt", header=TRUE)
corv <- vector(length = max(evolsel$RUN))
for (RUN in 1:max(evolsel$RUN))
{
  esR <- evolsel[evolsel$RUN==RUN,]
  corv[RUN] <- cor(esR$DeltaBV, esR$Sel)
}
plot(corv, ylim=c(-1,1))
plot(evolsel$DeltaBV, x=evolsel$Sel)
plot(evolsel$DeltaBV, x=evolsel$Beta)
plot(evolsel$DeltaBV, x=evolsel$Var)

library(lme4)
slm <- summary(lm(DeltaBV ~ Sel, data=evolsel))
slm$coefficients
coefficients(lm(DeltaBV ~ Sel, data=evolsel))
summary(lm(DeltaBVRecruits ~ SelRecruits, data=evolsel))

head(evolsel)

## With various h2 ####
h2values <- seq(0.1,1,0.1)
foldernames <- paste0("tryh2eq",h2values)
slopeSE<- data.frame(h2 =h2values, slope=NA,SE=NA)
for (i in 1:length(foldernames))
{
  evolsel <- read.table(file = paste0(foldernames[i],"/EvolSel.txt"), header=TRUE)
  slm <- summary(lm(DeltaBV ~ Sel, data=evolsel))
  slopeSE[i,2:3] <- slm$coefficients['Sel',c(1,2)]
}

plot(slopeSE[,1], slopeSE[,2])
arrows(x0=slopeSE[,1], y0 = slopeSE[,2]+1.96*slopeSE[,3], y1=slopeSE[,2]-1.96*slopeSE[,3],code = 3, angle=90)
abline(a = 0, b=1)
