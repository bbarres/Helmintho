###############################################################################
###############################################################################
#Analizing Helmintosporiose strain resistance tests from Arvalis field trial
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)

#set the working directory
setwd("~/work/Rfichiers/Githuber/Helmintho_data")

#load the global dataset
helmdat<-read.table("helmindata.txt",header=T,sep="\t")


###############################################################################
#Analysis for the boscalid
###############################################################################

bosc.dat<-helmdat[helmdat$active_substance=="boscalid",]

#modeling the dose response curve
bosc.m1<-drm(perc_croiss~dose,data=bosc.dat[bosc.dat$sample_ID=="17-001-00",],
             fct=LL.4())
summary(bosc.m1)
plot(bosc.m1)


bosc.m1.pop1<-drm(perc_croiss~dose,data=bosc.dat[bosc.dat$site_ID=="17-001",],
                  sample_ID,fct=LL.4())
plot(bosc.m1.pop1)

bosc.m1.pop1<-drm(perc_croiss~dose,data=bosc.dat[bosc.dat$site_ID=="17-002",],
                  sample_ID,fct=LL.3())
plot(bosc.m1.pop1)


###############################################################################
#Analysis for the bixafen
###############################################################################

bixa.dat<-helmdat[helmdat$active_substance=="bixafen",]

#modeling the dose response curve
bixa.m1<-drm(perc_croiss~dose,data=bixa.dat[bixa.dat$sample_ID=="17-001-01",],
             fct=LL.4())
summary(bixa.m1)
plot(bixa.m1)




###############################################################################
#END
###############################################################################


#comparison of different model
mselect(tavelure.m1,list(LL.3(),LN.3(),LN.4()))

#plot with 95% confidence interval
plot(tavelure.m1,broken=TRUE,type="confidence")
plot(tavelure.m1,broken=TRUE,add=TRUE)
#a simplier plot
plot(tavelure.m1,broken=TRUE)

#evaluating the ED50
ed50val<-ED(tavelure.m1,50,interval="delta")
#this is also working for a list of value for ED10, ED90...
ed_val<-ED(tavelure.m1,c(10,50,90),interval="delta")

