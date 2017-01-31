setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")
YearPheno <- read.table(file = "YearPheno.txt", header=T)
head(YearPheno)

repetBMI_A <- lmer(formula =BMI ~ 1 + (1|ID) + (1|Year), data=YearPheno[YearPheno$Age=="A",])
summary(repetBMI_A)
576.11/(576.11+508.32)#R2 in adults

repetBMI <- lmer(formula =BMI ~ 1 + Age + (1|ID) + (1|Year), data=YearPheno)
summary(repetBMI)
505.86/(505.86+753.02)#R2 at all ages
