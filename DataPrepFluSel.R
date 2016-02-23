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
RawPheno$RelativeJulianMinus<-RawPheno$Julian-max(RawPheno$Julian)

ground<-read.table(file="C:/Users/Thimothee Admin/Documents/Thesis/field ecology/GroundCov.txt",header=T)

RawPheno$rock_coverage<-NA
RawPheno$herb_coverage<-NA
for (i in 1:dim(RawPheno)[1])
{
  if (length(ground$rock_coverage[which(ground$X==RawPheno$X_Position[i] & ground$Y==RawPheno$Y_Position[i])])>0)
  {
    RawPheno$rock_coverage[i]<-ground$rock_coverage[which(ground$X==RawPheno$X_Position[i] & ground$Y==RawPheno$Y_Position[i])]
    RawPheno$herb_coverage[i]<-ground$herb_coverage[which(ground$X==RawPheno$X_Position[i] & ground$Y==RawPheno$Y_Position[i])]
  }
}

RawPheno$rock_coverage[which(is.na(RawPheno$rock_coverage))]<-mean(RawPheno$rock_coverage,na.rm=T)
RawPheno$herb_coverage[which(is.na(RawPheno$herb_coverage))]<-mean(RawPheno$herb_coverage,na.rm=T)


