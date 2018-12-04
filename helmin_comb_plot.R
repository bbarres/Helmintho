###############################################################################
###############################################################################
#Comparison of the different level of resistance by strains
###############################################################################
###############################################################################

#loading the libraries
library(data.table)
library(psych)
library(RColorBrewer)
library(scales)

#loading the dataset of the combined results of the 4 SDHI resistance test
EC50_pop<-read.table("data/EC50_byPOP.txt",header=TRUE)
colovec<-brewer.pal(4,"Set1")

###############################################################################
#plot with individual group by pop
###############################################################################

op<-par(mfrow=c(2,2),mar=c(6.1,5.1,3,1))

EC50bosc<-EC50_pop[EC50_pop$SA_ID=="boscalid",]
plot(EC50bosc$ED50/0.39,col=alpha(colovec[as.numeric(EC50bosc$pop_col)],0.5),
     main="boscalide",xlab="",ylab="FR",las=1, 
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3,pch=19)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)
box(lwd=3,bty="o")

EC50bixa<-EC50_pop[EC50_pop$SA_ID=="bixafen",]
plot(EC50bixa$ED50/0.08,col=alpha(colovec[as.numeric(EC50bixa$pop_col)],0.5),
     main="bixafène",xlab="",ylab="",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3,pch=19)
abline(0.08/0.08,0,col="green4",lwd=2)
abline(0.8/0.08,0,col="red",lwd=2)
box(lwd=3,bty="o")

EC50fluo<-EC50_pop[EC50_pop$SA_ID=="fluopyram",]
plot(EC50fluo$ED50/0.44,col=alpha(colovec[as.numeric(EC50fluo$pop_col)],0.5),
     main="fluopyram",xlab="Identifiant souche",ylab="FR",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3,pch=19)
abline(0.44/0.44,0,col="green4",lwd=2)
abline(4.4/0.44,0,col="red",lwd=2)
box(lwd=3,bty="o")

EC50flux<-EC50_pop[EC50_pop$SA_ID=="fluxapyroxade",]
plot(EC50flux$ED50/0.21,col=alpha(colovec[as.numeric(EC50flux$pop_col)],0.5),
     main="fluxapyroxade",xlab="Identifiant souche",ylab="",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3,pch=19)
abline(0.21/0.21,0,col="green4",lwd=2)
abline(2.1/0.21,0,col="red",lwd=2)
box(lwd=3,bty="o")

par(op)

#export 15 x 10 for the figure in the CIMA poster


###############################################################################
#plot with individual with FR
###############################################################################

op<-par(mfrow=c(2,2),mar=c(6.1,5.1,3,1))

plot(EC50bosc$ED50[order(EC50bosc$ED50)]/0.39,main="boscalide",
     xlab="",ylab="FR",las=1, 
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)
box(lwd=3,bty="o")

plot(EC50bixa$ED50[order(EC50bixa$ED50)]/0.08,main="bixafène",
     xlab="",ylab="",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3)
abline(0.08/0.08,0,col="green4",lwd=2)
abline(0.8/0.08,0,col="red",lwd=2)
box(lwd=3,bty="o")

plot(EC50fluo$ED50[order(EC50fluo$ED50)]/0.44,main="fluopyram",
     xlab="Identifiant souche",ylab="FR",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3)
abline(0.44/0.44,0,col="green4",lwd=2)
abline(4.4/0.44,0,col="red",lwd=2)
box(lwd=3,bty="o")

plot(EC50flux$ED50[order(EC50flux$ED50)]/0.21,main="fluxapyroxade",
     xlab="Identifiant souche",ylab="",las=1,
     cex=1.5,cex.lab=3,cex.axis=2,cex.main=3)
abline(0.21/0.21,0,col="green4",lwd=2)
abline(2.1/0.21,0,col="red",lwd=2)
box(lwd=3,bty="o")

par(op)

#export 15 x 10 for the figure in the CIMA poster


###############################################################################
#covariation between the different FR for the different SDHI
###############################################################################

EC50_comb<-dcast(EC50_pop,sample_ID~SA_ID,value.var=c("ED50"),fun=mean)
plot(EC50_comb$boscalid[order(EC50_comb$boscalid)]~
       EC50_comb$bixafen[order(EC50_comb$boscalid)])
pairs(EC50_comb[order(EC50_comb$boscalid),-c(1)])
pairs(EC50_comb[order(EC50_comb$boscalid),-c(1)],upper.panel=NULL,
      diag.panel=NULL,log=TRUE)
pairs.panels(EC50_comb[,-c(1)],cex=0.3)

#same comparison without the strains that are very resistant 
#to boscalid (FR>30) or only with those strains
pairs.panels(EC50_comb[EC50_comb$boscalid<30,-c(1)] )
pairs.panels(EC50_comb[EC50_comb$boscalid==30,-c(1)] )

###############################################################################
#END
###############################################################################