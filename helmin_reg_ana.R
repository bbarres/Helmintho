###############################################################################
###############################################################################
#Analizing Helmintosporiose strain resistance tests from Arvalis field trial
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#set the working directory
setwd("~/work/Rfichiers/Githuber/Helmintho_data")

#load the global dataset
helmdat<-read.table("helmindata.txt",header=T,sep="\t")


###############################################################################
#Analysis for the boscalid
###############################################################################

#subsetting the global dataset
bosc.dat<-helmdat[helmdat$active_substance=="boscalid",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
bosc_rez<-as.character(bosc.dat[bosc.dat$dose=="30" & bosc.dat$perc_croiss>50,
                                "sample_ID"])
REZbos<-data.frame("sample_ID"=bosc_rez,"ED50"=30)
#we limit the dataset to the sample that reach somehow a IC of 50%
bosc.dat<-bosc.dat[!(bosc.dat$sample_ID %in% bosc_rez),]
bosc.dat<-drop.levels(bosc.dat)
for (i in 1: dim(table(bosc.dat$sample_ID))[1]) {
  temp.m1<-drm(perc_croiss~dose,
      data=bosc.dat[bosc.dat$sample_ID==names(table(bosc.dat$sample_ID))[i],],
      fct=LL.4())
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(bosc.dat$sample_ID))[i],
                    "ED50"=temp[1])
  REZbos<-rbind(REZbos,tempx)
}

REZbos$ED50[REZbos$ED50>30]<-30
plot(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZbos,file="REZbos.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
abline(v=10,col="red",lwd=3)
#export to pdf 4.5 x 9 inches

###############################################################################
#Analysis for the bixafen
###############################################################################

#subsetting the global dataset
bixa.dat<-helmdat[helmdat$active_substance=="bixafen",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
bixa_rez<-as.character(bixa.dat[bixa.dat$dose=="30" & bixa.dat$perc_croiss>50,
                                "sample_ID"])
REZbix<-data.frame("sample_ID"=character(),"ED50"=numeric())
#we limit the dataset to the sample that reach somehow a IC of 50%
bixa.dat<-bixa.dat[!(bixa.dat$sample_ID %in% bixa_rez),]
bixa.dat<-drop.levels(bixa.dat)
for (i in 1: dim(table(bixa.dat$sample_ID))[1]) {
  temp.m1<-drm(perc_croiss~dose,
       data=bixa.dat[bixa.dat$sample_ID==names(table(bixa.dat$sample_ID))[i],],
       fct=LL.4())
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(bixa.dat$sample_ID))[i],
                    "ED50"=temp[1])
  REZbix<-rbind(REZbix,tempx)
}

REZbix$ED50[REZbix$ED50>30]<-30
plot(REZbix$ED50[order(REZbix$ED50)]/0.08,main="Bixafen",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.08/0.08,0,col="green4",lwd=2)
abline(0.8/0.08,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZbix,file="REZbix.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZbix$ED50[order(REZbix$ED50)]/0.08,main="Bixafen",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
abline(v=10,col="red",lwd=3)
#export to pdf 4.5 x 9 inches


###############################################################################
#Analysis for the fluopyram
###############################################################################

#subsetting the global dataset
fluo.dat<-helmdat[helmdat$active_substance=="fluopyram",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
fluo_rez<-as.character(fluo.dat[fluo.dat$dose=="30" & fluo.dat$perc_croiss>50,
                                "sample_ID"])
REZflo<-data.frame("sample_ID"=character(),"ED50"=numeric())
#we limit the dataset to the sample that reach somehow a IC of 50%
fluo.dat<-fluo.dat[!(fluo.dat$sample_ID %in% fluo_rez),]
fluo.dat<-drop.levels(fluo.dat)
for (i in 1: dim(table(fluo.dat$sample_ID))[1]) {
  temp.m1<-drm(perc_croiss~dose,
       data=fluo.dat[fluo.dat$sample_ID==names(table(fluo.dat$sample_ID))[i],],
       fct=LL.4())
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(fluo.dat$sample_ID))[i],
                    "ED50"=temp[1])
  REZflo<-rbind(REZflo,tempx)
}

plot(REZflo$ED50[order(REZflo$ED50)]/0.44,main="Fluopyram",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.44/0.44,0,col="green4",lwd=2)
abline(4.4/0.44,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZflo,file="REZflo.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZflo$ED50[order(REZflo$ED50)]/0.44,main="Fluopyram",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
abline(v=10,col="red",lwd=3)

#export to pdf 4.5 x 9 inches

###############################################################################
#Analysis for the fluxapyroxad
###############################################################################

#subsetting the global dataset
flux.dat<-helmdat[helmdat$active_substance=="fluxapyroxad",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
flux_rez<-as.character(flux.dat[flux.dat$dose=="30" & flux.dat$perc_croiss>50,
                                "sample_ID"])
REZflx<-data.frame("sample_ID"=flux_rez,"ED50"=30)
#we limit the dataset to the sample that reach somehow a IC of 50%
flux.dat<-flux.dat[!(flux.dat$sample_ID %in% flux_rez),]
flux.dat<-drop.levels(flux.dat)
for (i in 1: dim(table(flux.dat$sample_ID))[1]) {
  temp.m1<-drm(perc_croiss~dose,
               data=flux.dat[flux.dat$sample_ID==names(table(flux.dat$sample_ID))[i],],
               fct=LL.4())
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(flux.dat$sample_ID))[i],
                    "ED50"=temp[1])
  REZflx<-rbind(REZflx,tempx)
}

REZflx$ED50[REZflx$ED50>30]<-30
plot(REZflx$ED50[order(REZflx$ED50)]/0.21,main="Fluxapyroxade",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.21/0.21,0,col="green4",lwd=2)
abline(2.1/0.21,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZflx,file="REZflx.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZflx$ED50[order(REZflx$ED50)]/0.21,main="Fluxapyroxade",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     xlab="FR Classes",las=1,col=heat.colors(8)[c(8:1,rep(1,7))],ylim=c(0,450))
abline(v=10,col="red",lwd=3)

#export to pdf 4.5 x 9 inches


###############################################################################
#combined plot, figure for CIMA article
###############################################################################

op<-par(mfrow=c(2,2))

plot(REZbos$ED50[order(REZbos$ED50)]/0.39,main="boscalide",cex
     xlab="Identifiant souche", ylab="Facteur de résistance",las=1)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)

plot(REZbix$ED50[order(REZbix$ED50)]/0.08,main="bixafène",
     xlab="Identifiant souche", ylab="Facteur de résistance",las=1)
abline(0.08/0.08,0,col="green4",lwd=2)
abline(0.8/0.08,0,col="red",lwd=2)

plot(REZflo$ED50[order(REZflo$ED50)]/0.44,main="fluopyram",
     xlab="Identifiant souche", ylab="Facteur de résistance",las=1)
abline(0.44/0.44,0,col="green4",lwd=2)
abline(4.4/0.44,0,col="red",lwd=2)

plot(REZflx$ED50[order(REZflx$ED50)]/0.21,main="fluxapyroxade",
     xlab="Identifiant souche", ylab="Facteur de résistance",las=1)
abline(0.21/0.21,0,col="green4",lwd=2)
abline(2.1/0.21,0,col="red",lwd=2)

par(op)

#export to pdf 10 x 8 inches


###############################################################################
#plot with individual group by pop
###############################################################################

EC50_pop<-read.table("EC50_byPOP.txt",header=TRUE)

op<-par(mfrow=c(2,2),mar=c(2,2.5,3,1))
EC50bosc<-EC50_pop[EC50_pop$SA_ID=="boscalid",]
plot(EC50bosc$ED50,col=as.character(EC50bosc$pop_col),main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)

EC50bixa<-EC50_pop[EC50_pop$SA_ID=="bixafen",]
plot(EC50bixa$ED50,col=as.character(EC50bixa$pop_col),main="Bixafen")
abline(0.08,0,col="green3",lwd=2)
abline(0.8,0,col="orange3",lwd=2)

EC50fluo<-EC50_pop[EC50_pop$SA_ID=="fluopyram",]
plot(EC50fluo$ED50,col=as.character(EC50fluo$pop_col),main="Fluopyram")
abline(0.44,0,col="green3",lwd=2)
abline(4.4,0,col="orange3",lwd=2)

EC50flux<-EC50_pop[EC50_pop$SA_ID=="fluxapyroxade",]
plot(EC50flux$ED50,col=as.character(EC50flux$pop_col),main="Fluxapyroxade")
abline(0.21,0,col="green3",lwd=2)
abline(2.1,0,col="orange3",lwd=2)

par(op)


###############################################################################
#END
###############################################################################
