#### Pre instructions ####
setwd(dir = "~/thesis/Mass/FluctuatingSelectionForPubli/")
library("pedantics")
library('lme4')

#### Load data ####
RawPheno <- read.table(file = "ForFluctuatingSelectionRaw.txt",header=TRUE )
pedRaw <- read.table(file="Pedigree.txt",header=TRUE)

#### Edit data ####
names(pedRaw)<-c("id","dam","sire")
ped<-orderPed(pedRaw)

#add julian date
RawPheno$Origin=as.POSIXct(paste(paste(rep("01/01/",dim(RawPheno)[1])),RawPheno$Calendar_year,sep=""),format="%d/%m/%Y",tz="GMT")
RawPheno$Julian<-0
for (i in 1:nrow(RawPheno))
{
  RawPheno$Julian[i]<-round(julian(as.POSIXct(RawPheno$Date[i],format="%d/%m/%Y"),origin=as.POSIXct(RawPheno$Origin[i],format="%Y/%m/%d")),digit=0)
}
RawPheno$RelativeJulian<-RawPheno$Julian-min(RawPheno$Julian)



