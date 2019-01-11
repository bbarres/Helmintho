###############################################################################
###############################################################################
#Analizing Helmintosporiose strain resistance tests from Arvalis field trial
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#load the global dataset
helmdat<-read.table("data/helmindata2018.txt",header=T,sep="\t")


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
CI50sens<-mean(REZbos$ED50[order(REZbos$ED50)][1:80])
REZbos<-merge(REZbos,helmdat[helmdat$active_substance=="boscalid" & 
                             helmdat$dose==0,],by="sample_ID")
REZbos<-cbind(REZbos,"FR"=REZbos$ED50/CI50sens)

op<-par(mfrow=c(2,2))
plot(REZbos$ED50[order(REZbos$ED50)]/CI50sens,main="Boscalid",xlab="Souches ID",
     ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)

hist(REZbos$ED50[order(REZbos$ED50)]/CI50sens,main="Boscalid",xlab="FR Classes",
     breaks=c(0,5,10,15,20,25,30,35,40,45,50),
     las=1,col=heat.colors(10)[10:1],ylim=c(0,110))
abline(v=10,col="red",lwd=3)

plot(REZbos$ED50[order(REZbos$site_ID, REZbos$ED50)]/CI50sens,
     main="Boscalid",xlab="Souches ID",ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)
par(op)
#export to pdf 8 x 9 inches

#export the results file
write.table(REZbos,file="output/REZbos18.txt",quote=FALSE,sep="\t",
            row.names=FALSE)


###############################################################################
#Analysis for the bixafen
###############################################################################

#subsetting the global dataset
bixa.dat<-helmdat[helmdat$active_substance=="bixafene",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
bixa_rez<-as.character(bixa.dat[bixa.dat$dose=="30" & bixa.dat$perc_croiss>50,
                                "sample_ID"])
REZbix<-data.frame("sample_ID"=bixa_rez,"ED50"=30)
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
CI50sens<-mean(REZbix$ED50[order(REZbix$ED50)][1:80])
REZbix<-merge(REZbix,helmdat[helmdat$active_substance=="bixafene" & 
                               helmdat$dose==0,],by="sample_ID")
REZbix<-cbind(REZbix,"FR"=REZbix$ED50/CI50sens)

op<-par(mfrow=c(2,2))
plot(REZbix$ED50[order(REZbix$ED50)]/CI50sens,main="Bixafen",xlab="Souches ID",
     ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)

hist(REZbix$ED50[order(REZbix$ED50)]/CI50sens,main="Bixafen",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,140,150,200,250,260,280,290),
     las=1,col=heat.colors(17)[17:1],ylim=c(0,200),freq=TRUE)
abline(v=10,col="red",lwd=3)

plot(REZbix$ED50[order(REZbix$site_ID,REZbix$ED50)]/CI50sens,
     main="Bixafen",xlab="Souches ID",ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)
par(op)
#export to pdf 10 x 8 inches

#export the results file
write.table(REZbix,file="output/REZbix18.txt",quote=FALSE,sep="\t",
            row.names=FALSE)


###############################################################################
#Analysis for the fluopyram
###############################################################################

#subsetting the global dataset
fluo.dat<-helmdat[helmdat$active_substance=="fluopyram",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
fluo_rez<-as.character(fluo.dat[fluo.dat$dose=="30" & fluo.dat$perc_croiss>50,
                                "sample_ID"])
REZflo<-data.frame("sample_ID"=fluo_rez,"ED50"=30)
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

CI50sens<-mean(REZflo$ED50[order(REZflo$ED50)][1:80])
REZflo<-merge(REZflo,helmdat[helmdat$active_substance=="fluopyram" & 
                               helmdat$dose==0,],by="sample_ID")
REZflo<-cbind(REZflo,"FR"=REZflo$ED50/CI50sens)

op<-par(mfrow=c(2,2))
plot(REZflo$ED50[order(REZflo$ED50)]/CI50sens,main="Fluopyram",
     xlab="Souches ID",ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)

hist(REZflo$ED50[order(REZflo$ED50)]/CI50sens,
     main="Fluopyram",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,190,200),
     freq=TRUE,las=1,col=heat.colors(16)[16:1],ylim=c(0,200))
abline(v=10,col="red",lwd=3)

plot(REZflo$ED50[order(REZflo$site_ID,REZflo$ED50)]/CI50sens,
     main="Fluopyram",xlab="Souches ID",ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)
par(op)
#export to pdf 10 x 8 inches

#export the results file
write.table(REZflo,file="output/REZflo18.txt",quote=FALSE,sep="\t",
            row.names=FALSE)


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
CI50sens<-mean(REZflx$ED50[order(REZflx$ED50)][1:80])
REZflx<-merge(REZflx,helmdat[helmdat$active_substance=="fluxapyroxad" & 
                               helmdat$dose==0,],by="sample_ID")
REZflx<-cbind(REZflx,"FR"=REZflx$ED50/CI50sens)

op<-par(mfrow=c(2,2))
plot(REZflx$ED50[order(REZflx$ED50)]/CI50sens,main="Fluxapyroxade",xlab="Souches ID",
     ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)

hist(REZflx$ED50[order(REZflx$ED50)]/CI50sens,main="Fluxapyroxade",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170),
     xlab="FR Classes",las=1,col=heat.colors(17)[17:1],ylim=c(0,170))
abline(v=10,col="red",lwd=3)

plot(REZflx$ED50[order(REZflx$site_ID,REZflx$ED50)]/CI50sens,
     main="Fluxapyroxade",xlab="Souches ID",ylab="FR",las=1)
abline(CI50sens/CI50sens,0,col="green4",lwd=2)
abline((CI50sens*10)/CI50sens,0,col="red",lwd=2)
par(op)
#export to pdf 10 x 8 inches

#export the results file
write.table(REZflx,file="output/REZflx18.txt",quote=FALSE,sep="\t",
            row.names=FALSE)


###############################################################################
#END
###############################################################################
