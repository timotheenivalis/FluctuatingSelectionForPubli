#### Pre instructions ####
#setwd(dir = "~/thesis/Mass/FluctuatingSelectionForPubli/")
setwd(dir = "/home/timothee/Documents/GitHub/FluctuatingSelectionForPubli/")

library("pedantics")
library('lme4')
library('tidyr')

#### Load data ####
RawPheno <- read.table(file = "ForFluctuatingSelectionRaw.txt",header=TRUE )
pedRaw <- read.table(file="Pedigree.txt",header=TRUE)

#### Edit data ####
names(pedRaw)<-c("animal","dam","sire")
ped<-orderPed(pedRaw)
ped[ped==0]<-NA

write.table(ped,file = "ped.txt", quote = FALSE,row.names = FALSE)
#add julian date
RawPheno$Origin=as.POSIXct(paste(paste(rep("01/01/",dim(RawPheno)[1])),RawPheno$Calendar_year,sep=""),format="%d/%m/%Y",tz="GMT")
RawPheno$Julian<-0
for (i in 1:nrow(RawPheno))
{
  RawPheno$Julian[i]<-round(julian(as.POSIXct(RawPheno$Date[i],format="%d/%m/%Y"),origin=as.POSIXct(RawPheno$Origin[i],format="%Y/%m/%d")),digit=0)
}
RawPheno$RelativeJulian<-RawPheno$Julian-min(RawPheno$Julian)

m0 <- lm(Weight ~ 1+RelativeJulian+I(RelativeJulian^2), data=RawPheno[which(RawPheno$Age=="A" & !is.na(RawPheno$Weight)),])
summary(m0)

plot(RawPheno$Weight[which(RawPheno$Age=="A"& !is.na(RawPheno$Weight))],residuals(m0)+m0$coefficients[1])

RawPheno$Weight[which(RawPheno$Age=="A" & !is.na(RawPheno$Weight))] <- residuals(m0)+m0$coefficients[1]

#### Year*Ind based data ####
Ind_unique <- unique(as.character(RawPheno$ID_Individual))
Year_unique <- unique(RawPheno$Calendar_year)

mBL_TL <- lm(RawPheno$Body_Length ~ RawPheno$Tail_Length)

YearPheno <- data.frame(ID=as.character(Ind_unique[1]), 
                        Year=Year_unique[1],
                        Sex="Female",
                        Age="A",
                        Mass=0,
                        BL = 0,
                        BMI = 0,
                        RelativeJulian=0,
                        Phi=0,
                        Rho=0,
                        Fitness=0,
                        stringsAsFactors = FALSE)
count <- 1
for (i in 1:length(Ind_unique))
  {
  yearsind <- unique (RawPheno$Calendar_year[RawPheno$ID_Individual==Ind_unique[i]])
  offsprings <- unique(RawPheno$ID_Individual[as.character(RawPheno$Mother) == as.character(Ind_unique[i]) | 
                                            as.character(RawPheno$Father) == as.character(Ind_unique[i])] )
  for (j in 1:length(yearsind))
    {
      YearPheno[count,"ID"] <- Ind_unique[i]
      YearPheno[count,"Year"] <-  yearsind[j]
      YearPheno[count,"Age"] <- as.character(RawPheno$Age[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                                  RawPheno$Calendar_year==  yearsind[j])][1])
      YearPheno[count,"Mass"] <- mean(RawPheno$Weight[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                              RawPheno$Calendar_year==  yearsind[j])], na.rm=TRUE)
      YearPheno[count,"BL"] <- 0.001*mean(RawPheno$Body_Length[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                              RawPheno$Calendar_year==  yearsind[j])], na.rm=TRUE)
      if (yearsind[j]==2006 & all(is.na(RawPheno$Body_Length[which(RawPheno$ID_Individual==Ind_unique[i] & RawPheno$Calendar_year==  yearsind[j])])))
      {
        YearPheno[count,"BMI"] <-  1000*mean(RawPheno$Weight[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                                     RawPheno$Calendar_year==  yearsind[j])] / (coefficients(mBL_TL)[1]+coefficients(mBL_TL)[2]*RawPheno$Tail_Length[which(RawPheno$ID_Individual==Ind_unique[i] & RawPheno$Calendar_year==  yearsind[j])]), na.rm=TRUE)
      }else{
        YearPheno[count,"BMI"] <-  1000*mean(RawPheno$Weight[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                              RawPheno$Calendar_year==  yearsind[j])] / RawPheno$Body_Length[which(RawPheno$ID_Individual==Ind_unique[i] & RawPheno$Calendar_year==  yearsind[j])], na.rm=TRUE)
        }
      YearPheno[count,"RelativeJulian"] <- mean(RawPheno$RelativeJulian[which(RawPheno$ID_Individual==Ind_unique[i] &
                                                             RawPheno$Calendar_year==  yearsind[j])], na.rm=TRUE)
      YearPheno[count,"Rho"] <- sum(as.character(offsprings) %in% as.character(RawPheno$ID_Individual[RawPheno$Cohort==yearsind[j]]))
      count <- count +1
    }#end for (j in 1:length(Year_unique))
  YearPheno[YearPheno$ID==Ind_unique[i],"Sex"] <- as.character(RawPheno$Sex[RawPheno$ID_Individual==Ind_unique[i]][1])
  
}# end for (i in 1:length(Ind_unique))

# correct adult mass for time


# Survival to the next year
YearPheno$Phi <- c(as.numeric(YearPheno$ID[-1] == YearPheno$ID[-length(YearPheno$ID)]),0)
#YearPheno$Phi[YearPheno$Year==max(YearPheno$Year)] <- NA
# but always exclude the last year from survival analyses

YearPheno$Fitness <- YearPheno$Rho + YearPheno$Phi * 2

YearPheno$StMass <- (YearPheno$Mass - mean(YearPheno$Mass, na.rm=TRUE))/sd(YearPheno$Mass, na.rm=TRUE)

YearPheno$Mother <- NA
for (i in as.character(ped$animal[!is.na(ped$dam)]))
{
  YearPheno$Mother[YearPheno$ID==i] <- as.character(ped$dam[ped$animal==i])
}

YearPhenoIntermediate <- YearPheno
################## Must include age-corrected mass

AtoMerge <- read.table(file="AtoMerge.txt", header = TRUE)

YearPheno<- merge(x=YearPhenoIntermediate, y=AtoMerge,by.x = "ID", by.y="id",all.x = TRUE)

YearPheno$A[YearPheno$Age=="A"] <- YearPheno$Mass[YearPheno$Age=="A"]

YearPheno$animal <- YearPheno$ID


YearPheno$M2006 <- NA
YearPheno$M2007 <- NA
YearPheno$M2008 <- NA
YearPheno$M2009 <- NA
YearPheno$M2010 <- NA
YearPheno$M2011 <- NA
YearPheno$M2012 <- NA
YearPheno$M2013 <- NA
YearPheno$M2014 <- NA
YearPheno$M2015 <- NA

for (i in 1:nrow(YearPheno))
{
  YearPheno[i, paste("M",YearPheno$Year[i],sep="")] <- YearPheno$Mass[i]
}

YearPheno$A2006 <- NA
YearPheno$A2007 <- NA
YearPheno$A2008 <- NA
YearPheno$A2009 <- NA
YearPheno$A2010 <- NA
YearPheno$A2011 <- NA
YearPheno$A2012 <- NA
YearPheno$A2013 <- NA
YearPheno$A2014 <- NA
YearPheno$A2015 <- NA

for (i in 1:nrow(YearPheno))
{
  YearPheno[i, paste("A",YearPheno$Year[i],sep="")] <- YearPheno$A[i]
}

YearPheno$A1 <- NA
YearPheno$A2 <- NA

YearPheno$BL1 <- NA
YearPheno$BL2 <- NA

YearPheno$BMI1 <- NA
YearPheno$BMI2 <- NA

for (i in 1:nrow(YearPheno))
{
  if (YearPheno$Year[i] %in% c(2006,2007,2009,2010,2013,2014))
  { YearPheno$A1[i] <- YearPheno$A[i]
  } else{YearPheno$A2[i] <- YearPheno$A[i]
  }
  
  if (YearPheno$Year[i] %in% c(2006, 2007, 2008, 2009, 2013, 2014))
    {YearPheno$BL1[i] <- YearPheno$BL[i]
    YearPheno$BMI1[i] <- YearPheno$BMI[i]
    }else{
      YearPheno$BL2[i] <- YearPheno$BL[i]
      YearPheno$BMI2[i] <- YearPheno$BMI[i]
      }
  }

YearPheno$RJst <- (YearPheno$RelativeJulian - mean(YearPheno$RelativeJulian))/sd(YearPheno$RelativeJulian)
YearPheno$RJ2st <- (YearPheno$RelativeJulian^2 - mean(YearPheno$RelativeJulian^2))/sd(YearPheno$RelativeJulian^2)

YearPheno$Ast <- (YearPheno$A - mean(YearPheno$A, na.rm=T))/sd(YearPheno$A, na.rm=T)
YearPheno$BLst <-  (YearPheno$BL - mean(YearPheno$BL, na.rm=T))/sd(YearPheno$BL, na.rm=T)
YearPheno$BMIst <-  (YearPheno$BMI - mean(YearPheno$BMI, na.rm=T))/sd(YearPheno$BMI, na.rm=T)

yearmeanfit <- tapply(YearPheno$Fitness,YearPheno$Year,mean)
yearmeanfit["2006"]
for (i in names(yearmeanfit))
{
  YearPheno$FitnessYear[YearPheno$Year==i] <- YearPheno$Fitness[YearPheno$Year==i]/ yearmeanfit[i]
}

YearPheno$OS <- NA
for(i in 1:nrow(YearPheno))
{
  offsp <- ped$animal[which(ped$dam==YearPheno$ID[i] | ped$sire==YearPheno$ID[i])]
  val <- mean(YearPheno$A[which(YearPheno$ID%in%offsp & YearPheno$Age=="J" & YearPheno$Year==YearPheno$Year[i])])
  if(!is.nan(val) & !is.na(val)){
  YearPheno$OS[i] <- val
  }
}


write.table(x = YearPheno, file = "YearPheno.txt",quote = FALSE, row.names = FALSE)

hist(YearPheno$BMI[YearPheno$Year==2006])
hist(YearPheno$BMI[YearPheno$Year!=2006], add=T)


plot(YearPheno$BMI[YearPheno$Year!=2006], YearPheno$FitnessYear[YearPheno$Year!=2006]);cor.test(YearPheno$BMI[YearPheno$Year!=2006], YearPheno$FitnessYear[YearPheno$Year!=2006])
plot(YearPheno$BMI[YearPheno$Year==2006], YearPheno$FitnessYear[YearPheno$Year==2006]);cor.test(YearPheno$BMI[YearPheno$Year==2006], YearPheno$FitnessYear[YearPheno$Year==2006])
plot(YearPheno$BMI, YearPheno$FitnessYear);cor.test(YearPheno$BMI, YearPheno$FitnessYear)

plot(YearPheno$BL, YearPheno$FitnessYear);cor.test(YearPheno$BL, YearPheno$FitnessYear)
