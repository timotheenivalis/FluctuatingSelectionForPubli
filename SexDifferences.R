YearPheno$BMI1 <- (YearPheno$BMI1 - mean(YearPheno$BMI1, na.rm=TRUE))/sd(YearPheno$BMI1,na.rm = TRUE)
YearPheno$BMI2 <- (YearPheno$BMI2 - mean(YearPheno$BMI2, na.rm=TRUE))/sd(YearPheno$BMI2,na.rm = TRUE)
YP <- YearPheno[,c("ID", "Year", "Sex", "Age", "BMI", "Mother", "animal", "BMI1", "BMI2", "RJst", "RJ2st", "BMIst", "PhiZ", "RhoZ", "FitnessZ", "GGImm")]

write.table(YP, file = "YearPhenoSimple.txt", quote = FALSE, row.names = FALSE)
YP$BMI1


mmARRfitness <- glmer(Fitness ~ 1 + BMIst + Sex *  Age  + Age*RJst +(1+BMIst|Year),
                      data=YearPheno, family=poisson, na.action = "na.omit")


mmARRfitnessSEX <- glmer(Fitness ~ 1 + BMIst + Sex *  Age  + Age*RJst +(1+BMIst|Sex:Year),data=YearPheno, family=poisson, na.action = "na.omit")
summary(mmARRfitnessSEX)                      

mmARRfitnessF <- glmer(FitnessZ ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                      data=YearPheno[YearPheno$Sex=="Female",], family=poisson, na.action = "na.omit")

mmARRfitnessM <- glmer(FitnessZ ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                       data=YearPheno[YearPheno$Sex=="Male",], family=poisson, na.action = "na.omit")
summary(mmARRfitnessF)
summary(mmARRfitnessM)


mmIfitnessF <- glmer(FitnessZ ~ 1 + BMIst  + Age*RJst +(1|Year),
                     data=YearPheno[YearPheno$Sex=="Female",], family=poisson, na.action = "na.omit")
mmIfitnessM <- glmer(FitnessZ ~ 1 + BMIst  + Age*RJst +(1|Year),
                     data=YearPheno[YearPheno$Sex=="Male",], family=poisson, na.action = "na.omit")

anova(mmARRfitnessF, mmIfitnessF)
anova(mmARRfitnessM, mmIfitnessM)



mmIfitnessFrho <- glmer(RhoZ ~ 1 + BMIst  + Age*RJst +(1|Year),
                        data=YearPheno[YearPheno$Sex=="Female",], family=poisson, na.action = "na.omit")
mmRfitnessFrho <- glmer(RhoZ ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                        data=YearPheno[YearPheno$Sex=="Female",], family=poisson, na.action = "na.omit")
summary(mmRfitnessFrho)

mmIfitnessMrho <- glmer(RhoZ ~ 1 + BMIst  + Age*RJst +(1|Year),
                        data=YearPheno[YearPheno$Sex=="Male",], family=poisson, na.action = "na.omit")
mmRfitnessMrho <- glmer(RhoZ ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                        data=YearPheno[YearPheno$Sex=="Male",], family=poisson, na.action = "na.omit")
summary(mmRfitnessMrho)
anova(mmRfitnessFrho, mmIfitnessFrho)
anova(mmRfitnessMrho, mmIfitnessMrho)


mmIfitnessFphi <- glmer(Phi ~ 1 + BMIst  + Age*RJst +(1|Year),
                        data=YearPheno[YearPheno$Sex=="Female",], family=binomial, na.action = "na.omit")
mmRfitnessFphi <- glmer(Phi ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                        data=YearPheno[YearPheno$Sex=="Female",], family=binomial, na.action = "na.omit")
summary(mmRfitnessFphi)

mmIfitnessMphi <- glmer(Phi ~ 1 + BMIst  + Age*RJst +(1|Year),
                        data=YearPheno[YearPheno$Sex=="Male",], family=binomial, na.action = "na.omit")
mmRfitnessMphi <- glmer(Phi ~ 1 + BMIst  + Age*RJst +(1+BMIst|Year),
                        data=YearPheno[YearPheno$Sex=="Male",], family=binomial, na.action = "na.omit")
summary(mmRfitnessMphi)

anova(mmRfitnessFphi, mmIfitnessFphi)
anova(mmRfitnessMphi, mmIfitnessMphi)
