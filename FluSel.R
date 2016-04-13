YearPheno <- read.table(file = "YearPheno.txt", header=T)

SelByYear <- vector(length = 2015-2006)
SeSelByYear <- vector(length = 2015-2006)
for (t in 2006:2015)
  {
  m0 <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno[YearPheno$Year==t,], family=poisson)
  SelByYear[t-2005] <- coefficients(m0)[2]
  sm0<-summary(m0)
  SeSelByYear[t-2005] <- sm0$coefficients[2,2]
}
plot(SelByYear, x=2006:2015, ylim=c(-1.5,0.8), xlab="Year", ylab = "Selection gradient")
abline(h=0)
arrows(x0 = 2006:2015,x1 = 2006:2015,code = 3, y0 = SelByYear+1.96*SeSelByYear,
       y1 = SelByYear-1.96*SeSelByYear, angle = 90,length = 0.1)
m0all <- glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=poisson)
abline(h=coefficients(m0all)[2], lty=2)
sm0all <- summary(m0all)
lowm0all <- coefficients(m0all)[2]+1.96*sm0all$coefficients[2,2]
highm0all <- coefficients(m0all)[2]-1.96*sm0all$coefficients[2,2]
polygon(x=c(2005,2016,2016,2005),y=c(lowm0all,lowm0all, highm0all, highm0all),
        fillOddEven = TRUE, col=rgb(0.1,0.1,0.1,0.5), lty=2, )



summary(glm(Phi ~ 1 + Mass + Sex , data=YearPheno[YearPheno$Age=="A" & YearPheno$Year<max(YearPheno$Year),], family=binomial))
summary(glmer(Phi ~ 1 + as.factor(Year) + Mass + Sex +(0+Mass|Year), data=YearPheno[YearPheno$Age=="A" YearPheno$Year<max(YearPheno$Year),], family=binomial))

summary(glm(Rho ~ 1 + StMass + Sex , data=YearPheno, family=poisson))
summary(glmer(Rho ~ 1 + as.factor(Year) + StMass + Sex +(0+Mass|Year), data=YearPheno, family=poisson))


summary(glm(Fitness ~ 1 + StMass + Sex +Age , data=YearPheno, family=poisson))
summary(glmer(Fitness ~ 1 + as.factor(Year) + StMass + Sex +(0+Mass|Year), data=YearPheno, family=poisson))
summary(glmer(Fitness ~ 1 + StMass + Sex + Age +(1+Mass|Year), data=YearPheno[!is.na(YearPheno$StMass),], family=poisson))

prior0 <- list(G=list(G1= list(V=0.00001*diag(2), nu=2)), R=list(V=diag(1),nu=1))
m0 <- MCMCglmm(fixed = Fitness ~ 1 + Sex + Age + StMass,
               random = ~ us(1+StMass):Year, prior= prior0, family="poisson",
               data=YearPheno[!is.na(YearPheno$StMass),], nitt = 110000, burnin = 10000, thin = 100)
summary(m0)
autocorr(m0$VCV)
plot(m0)
