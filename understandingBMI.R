plot(YearPheno$BMI, YearPheno$RelativeJulian) 

plot(YearPheno$BMI[YearPheno$Age=="J"], x=YearPheno$RelativeJulian[YearPheno$Age=="J"], col="blue")
points(YearPheno$BMI[YearPheno$Age=="A"], x=YearPheno$RelativeJulian[YearPheno$Age=="A"], col="red")

plot(YearPheno$BMI[YearPheno$Age=="J"], YearPheno$RelativeJulian[YearPheno$Age=="J"], xlim=c(150,350))

cor.test(YearPheno$BMI[YearPheno$Age=="J"], YearPheno$RelativeJulian[YearPheno$Age=="J"], xlim=c(150,350))

summary(lmer(formula = BMI ~ 1 +Sex + RJst + (1|ID), data=YearPheno))

RawPheno$BMI <- 1000*RawPheno$Weight/RawPheno$Body_Length
RawPheno$RJst <- (RawPheno$RelativeJulian-mean(RawPheno$RelativeJulian, na.rm=TRUE))/sd(RawPheno$RelativeJulian,na.rm = TRUE)
RawPheno$meanRJst <- NA
for (i in 1:nrow(RawPheno))
{
  if (RawPheno$Age[i]=="J")
  { 
    RawPheno$meanRJst[i] <- mean(RawPheno$RJst[RawPheno$ID_Individual==RawPheno$ID_Individual[i] & RawPheno$Age=="J"], na.rm=TRUE)
  }
}

RawPheno$diffRJst <- RawPheno$RJst - RawPheno$meanRJst
head(RawPheno)
summary(lmer(formula = BMI ~ 1 +Sex +RJst + I(RJst^2)+ (1|ID_Individual), data=RawPheno[RawPheno$Age=="J",]))
summary(lmer(formula = BMI ~ 1 +Sex + diffRJst + meanRJst + I(diffRJst^2) + I(meanRJst^2) + (1|ID_Individual), data=RawPheno[RawPheno$Age=="J",]))

summary(lmer(formula = Weight ~ 1 +Sex + diffRJst + meanRJst + I(diffRJst^2) + I(meanRJst^2) + (1|ID_Individual), data=RawPheno[RawPheno$Age=="J",]))


summary(lmer(formula = BMI ~ 1 +Sex + Weight + (1|ID_Individual), data=RawPheno[RawPheno$Age=="J",]))
summary(lmer(formula = BMI ~ 1 +Sex + Weight + (1|ID_Individual), data=RawPheno[RawPheno$Age=="A",]))
