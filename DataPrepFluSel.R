#### Pre instructions ####
setwd(dir = "~/thesis/Mass/FluctuatingSelectionForPubli/")
setwd(dir = "D:/Documents/GitHub/FluctuatingSelectionForPubli/")

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



#### Year*Ind based data ####
Ind_unique <- unique(as.character(RawPheno$ID_Individual))
Year_unique <- unique(RawPheno$Calendar_year)

YearPheno <- data.frame(ID=as.character(Ind_unique[1]), 
                        Year=Year_unique[1],
                        Sex="Female",
                        Age="A",
                        Mass=0,
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
      YearPheno[count,"Rho"] <- sum(as.character(offsprings) %in% as.character(RawPheno$ID_Individual[RawPheno$Cohort==yearsind[j]]))
      count <- count +1
    }#end for (j in 1:length(Year_unique))
  YearPheno[YearPheno$ID==Ind_unique[i],"Sex"] <- as.character(RawPheno$Sex[RawPheno$ID_Individual==Ind_unique[i]][1])
  
}# end for (i in 1:length(Ind_unique))

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

################## Must include age-corrected mass




write.table(x = YearPheno, file = "YearPheno.txt",quote = FALSE, row.names = FALSE)

