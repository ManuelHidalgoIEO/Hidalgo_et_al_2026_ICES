########################################################
# Portfolio analyses from worldwide surveys (WGCOMEDA) #
########################################################

quartz.options(dpi=75)

library(tidyverse)
library(vegan)
library(ecofolio)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(semTools)
library(patchwork)
library(cowplot)
library(car)
library(maps)
library(nlme)
library(ggrepel)
library(corrr)
library(ggcorrplot)

MyStd<-function (x) {(x-mean(x,na.rm=T))/sd(x,na.rm=T)}

m.world<-map_data("world")

setwd() # set working directory 

# ------------------------------------------------------ #
# ---------------------- STEP 1 ------------------------ #
# ------------------------------------------------------ #
#1# Compute metrics from all the species-abundance time series in each ecosystem

#1.1# Load species data

allspp<-read.csv2(file="sppdivALL2025.csv",header=T,dec=",",sep=",")

head(allspp)

allspp<-allspp[,-5]

allspp<-filter(allspp,year >= 1980) # avoid years before 1980 (n=7)

#1.2# Loop for calculating all required metrics

allAreasData<-matrix(data=NA,nrow=54,ncol=15)
# empty matrix with 54 rows (studied systems) and 15 columns (metrics to be calculated from each system)

ldf<-list() # make a list to store data.frames of time series

for (i in 1:54) {

  area1<-filter(allspp,map==i) # select an ecosystem

  min.yr<-min(area1$year) # start year
  max.yr<-max(area1$year) # end year

  sp.names<-unique(area1$species) # species present
  sp.num<-length(sp.names) # total number of species 

  template.yr<-data.frame(map=rep(unique(area1$map),times=length(seq(min.yr,max.yr,1))*sp.num),
                          region=rep(unique(area1$region),times=length(seq(min.yr,max.yr,1))*sp.num),
                          area=rep(unique(area1$area),times=length(seq(min.yr,max.yr,1))*sp.num),
                          species=rep(sp.names,each=length(seq(min.yr,max.yr,1))),
                          year=rep(seq(min.yr,max.yr,1),times=sp.num))
  # build a complete data.frame to account for the zero catch not included in the original data.frames

  area1<-merge(template.yr,area1,by.x=c("map","region","area","species","year"),
               by.y=c("map","region","area","species","year"),all=T) # merge data.frame with complete template

  area1$cpue<-ifelse(is.na(area1$cpue),0,area1$cpue) # NAs to zero catch

  area1<-area1[order(area1$species,area1$year),] # to be sure of the order by species and year 

  area1.wide<-data.frame(area1 %>% pivot_wider(names_from=species,values_from=cpue)) # reshape to wide format

  area1.wide$abundanceTot<-rowSums(area1.wide[,5:(sp.num+4)],na.rm=T) # calculate total abundance
  area1.wide$richness<-specnumber(area1.wide[,5:(sp.num+4)]) # richness
  area1.wide$shannon<-diversity(area1.wide[,5:(sp.num+4)],"shannon") # Shannon diversity
  area1.wide$evenness<-area1.wide$shannon/log(area1.wide$richness) # evenness

  area1.wide<-filter(area1.wide,abundanceTot > 0) # remove years with no survey data

  z.value<-fit_taylor(area1.wide[,5:(sp.num+4)],ci=T) # calculate Taylor Power Law z.value and CI
  avcvpe<-pe_avg_cv(area1.wide[,5:(sp.num+4)],detrending="not_detrended",boot_reps=500,ci=T) # Average_CV_PE and CI
  #avcvpe<-pe_avg_cv(area1.wide[,5:(sp.num+4)],detrending="linear_detrended",boot_reps=500,ci=T) # detrending included
  mvpe<-pe_mv(area1.wide[,5:(sp.num+4)],type="linear",ci=T) # MV_PE and CI
  #mvpe<-pe_mv(area1.wide[,5:(sp.num+4)],type="linear_detrended",ci=T) # detrending included

  # populate the matrix with all the metrics for each ecosystem
  allAreasData[i,1]<-mean(area1.wide$abundanceTot)
  allAreasData[i,2]<-sd(area1.wide$abundanceTot)
  allAreasData[i,3]<-mean(area1.wide$richness)
  allAreasData[i,4]<-mean(area1.wide$shannon)
  allAreasData[i,5]<-mean(area1.wide$evenness)
  allAreasData[i,6]<-z.value$z
  allAreasData[i,7]<-z.value$z.l
  allAreasData[i,8]<-z.value$z.u
  allAreasData[i,9]<-synchrony(area1.wide[,5:(sp.num+4)])
  allAreasData[i,10]<-avcvpe$pe
  allAreasData[i,11]<-avcvpe$ci[1]
  allAreasData[i,12]<-avcvpe$ci[2]
  allAreasData[i,13]<-mvpe$pe
  allAreasData[i,14]<-mvpe$ci[1]
  allAreasData[i,15]<-mvpe$ci[2]

  # extract the time series for abundance and diversity indices
  ldf[[i]]<-area1.wide[,c("map","region","area","year",
                          "abundanceTot","richness","shannon","evenness")]

}

#1.3# Set up ecosystem metrics data.frame

colnames(allAreasData)<-c("cpue","cpueSD","richness","shannon","evenness",
                          "z.value","z.value.l","z.value.u","synchrony",
                          "avcvpe","avcvpe.l","avcvpe.u",
                          "mvpe","mvpe.l","mvpe.u") # column names

allAreasData<-data.frame(allAreasData) # transform into a data.frame

allAreasData$area<-unique(allspp$area) # add ecosystem name
allAreasData$map<-seq(1,54,1) # add number to each ecosystem 

#1.4# Set up ecosystem time series data.frame

time.series.df<-bind_rows(ldf) # bind all time series files by row

#1.5# Example figure to calculate PE metrics for the Western Scotian Shelf

ex.pe.long<-filter(allspp,area=="Western Scotian Shelf") # select region

ex.pe.wide<-data.frame(ex.pe.long %>% # transform to wide format
                         pivot_wider(names_from=species,values_from=cpue))

ex.pe.wide<-ex.pe.wide %>% mutate_at(c(5:28),~replace_na(.,0)) # replace NA with 0 catch

mean.sel.sp<-filter(ex.pe.long,species=="Lophius americanus" |
                      species=="Hippoglossoides platessoides" |
                      species=="Melanogrammus aeglefinus") %>%
  group_by(species) %>% summarise(mean_pts=mean(log(cpue))) # mean values for each sp

ex1<-ggplot(data=filter(time.series.df,area=="Western Scotian Shelf"),
       aes(x=year,y=log(abundanceTot)))+
  geom_line(col="black")+
  geom_hline(aes(yintercept=mean(log(abundanceTot))),col="black")+
  scale_x_continuous("",breaks=c(1980,1985,1990,1995,2000,2005,2010,2015))+
  scale_y_continuous("log(Total Abundance)",limits=c(8.2,9.7),
                     breaks=c(8.2,8.5,8.8,9.1,9.4,9.7))+
  theme_bw()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

ex2<-ggplot(data=filter(ex.pe.long,species=="Lophius americanus" |
                     species=="Hippoglossoides platessoides" |
                     species=="Melanogrammus aeglefinus"),
       aes(x=year,y=log(cpue),col=species))+
  geom_line()+
  geom_hline(data=mean.sel.sp,aes(yintercept=mean_pts,col=species))+
  scale_x_continuous("Year",breaks=c(1980,1985,1990,1995,2000,2005,2010,2015))+
  scale_y_continuous("log(Abundance)",limits=c(1,9),breaks=c(1,3,5,7,9))+
  theme_bw()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

pdf(file="Ex_abund.pdf",width=6,height=6)
ex1/ex2
dev.off()

pdf(file="Ex_PE.pdf",width=6,height=6)
plot_mv(ex.pe.wide[,-c(1:4)],show="linear",cex.lab=1.4,cex.axis=1.2,cex=2)
mtext("(C)",adj=0,line=1,cex=1.2)
points(x=mean(ex.pe.wide$Lophius.americanus),
       y=var(ex.pe.wide$Lophius.americanus),col="#00BA38",pch=16,cex=2)
points(x=mean(ex.pe.wide$Hippoglossoides.platessoides),
       y=var(ex.pe.wide$Hippoglossoides.platessoides),col="#F8766D",pch=16,cex=2)
points(x=mean(ex.pe.wide$Melanogrammus.aeglefinus),
       y=var(ex.pe.wide$Melanogrammus.aeglefinus),col="#619CFF",pch=16,cex=2)
dev.off()


# ------------------------------------------------------ #
# ---------------------- STEP 2 ------------------------ #
# ------------------------------------------------------ #
#2# Q1: Do the environmental conditions relate to the average CPUE,
# and taxonomic and functional diversity across ecosystems?

#2.1# Load data and explore locations

portfolio<-read.table(file="portfolioDataForGEB.txt",header=T,dec=".",sep=",")

head(portfolio)

#2.2# Plot mid location of each system 

ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray60",col="black")+
  coord_cartesian(xlim=c(-95,40),ylim=c(15,70))+
  geom_point(size=4,col=2)+
  geom_text_repel(aes(x=Longitude,y=Latitude,label=System),col="blue")

#2.3# Examine simple correlations between CPUE, taxonomic and functional
# diversity with the environmental conditions 

portfolio %>% 
  dplyr::select(c(18,20,22,7,8,9,10:15)) %>% 
  correlate(use="pairwise.complete.obs") %>% 
  focus(cpue,richness,evenness,K,L50,TL)

#2.4# Get coefficients from linear models and plot for CPUE

cpueData<-portfolio[,c("cpue","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(cpueData)[2:7]

cpueCoefs<-matrix(data=NA,nrow=6,ncol=3)

cpueCoefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(log(cpue)~MyStd(i),list(i=as.name(x))),data=cpueData))[2]
}))

cpueCoefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(log(cpue)~MyStd(i),list(i=as.name(x))),data=cpueData))[2]
}))

cpueCoefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(log(cpue)~MyStd(i),list(i=as.name(x))),data=cpueData))[4]
}))

cpueCoefs<-as.data.frame(cpueCoefs)
colnames(cpueCoefs)<-c("coef","coefLow","coefUpp")
cpueCoefs$Predictor<-varlist

c1<-ggplot(data=cpueCoefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("red","black","black","black","black","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("red","black","black","black","black","black"))+
  scale_x_continuous("")+
  scale_y_continuous("Environmetal parameter",breaks=c(1,2,3,4,5,6),
                     labels=c("Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

#2.5# Get coefficients from linear models and plot for Richness

richData<-portfolio[,c("richness","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(richData)[2:7]

richCoefs<-matrix(data=NA,nrow=6,ncol=3)

richCoefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(richness~MyStd(i),list(i=as.name(x))),data=richData))[2]
}))

richCoefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(richness~MyStd(i),list(i=as.name(x))),data=richData))[2]
}))

richCoefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(richness~MyStd(i),list(i=as.name(x))),data=richData))[4]
}))

richCoefs<-as.data.frame(richCoefs)
colnames(richCoefs)<-c("coef","coefLow","coefUpp")
richCoefs$Predictor<-varlist

c2<-ggplot(data=richCoefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("black","black","black","black","red","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("black","black","black","black","red","black"))+
  scale_x_continuous("")+
  scale_y_continuous("",breaks=c(1,2,3,4,5,6),labels=c("","","","","",""))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

#2.6# Get coefficients from linear models and plot for Evenness

evenData<-portfolio[,c("evenness","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(evenData)[2:7]

evenCoefs<-matrix(data=NA,nrow=6,ncol=3)

evenCoefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(evenness~MyStd(i),list(i=as.name(x))),data=evenData))[2]
}))

evenCoefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(evenness~MyStd(i),list(i=as.name(x))),data=evenData))[2]
}))

evenCoefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(evenness~MyStd(i),list(i=as.name(x))),data=evenData))[4]
}))

evenCoefs<-as.data.frame(evenCoefs)
colnames(evenCoefs)<-c("coef","coefLow","coefUpp")
evenCoefs$Predictor<-varlist

c3<-ggplot(data=evenCoefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("black","black","black","black","red","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("black","black","black","black","red","black"))+
  scale_x_continuous("")+
  scale_y_continuous("",breaks=c(1,2,3,4,5,6),labels=c("","","","","",""))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#2.7# Get coefficients from linear models and plot for K

kData<-portfolio[,c("K","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(kData)[2:7]

kCoefs<-matrix(data=NA,nrow=6,ncol=3)

kCoefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(K~MyStd(i),list(i=as.name(x))),data=kData))[2]
}))

kCoefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(K~MyStd(i),list(i=as.name(x))),data=kData))[2]
}))

kCoefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(K~MyStd(i),list(i=as.name(x))),data=kData))[4]
}))

kCoefs<-as.data.frame(kCoefs)
colnames(kCoefs)<-c("coef","coefLow","coefUpp")
kCoefs$Predictor<-varlist

c4<-ggplot(data=kCoefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("black","black","black","black","red","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("black","black","black","black","red","black"))+
  scale_x_continuous("Standadized effect")+
  scale_y_continuous("Environmental parameter",breaks=c(1,2,3,4,5,6),
                     labels=c("Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

#2.8# Get coefficients from linear models and plot for L50

l50Data<-portfolio[,c("L50","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(l50Data)[2:7]

l50Coefs<-matrix(data=NA,nrow=6,ncol=3)

l50Coefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(L50~MyStd(i),list(i=as.name(x))),data=l50Data))[2]
}))

l50Coefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(L50~MyStd(i),list(i=as.name(x))),data=l50Data))[2]
}))

l50Coefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(L50~MyStd(i),list(i=as.name(x))),data=l50Data))[4]
}))

l50Coefs<-as.data.frame(l50Coefs)
colnames(l50Coefs)<-c("coef","coefLow","coefUpp")
l50Coefs$Predictor<-varlist

c5<-ggplot(data=l50Coefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("black","black","black","black","red","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("black","black","black","black","red","black"))+
  scale_x_continuous("Standadized effect")+
  scale_y_continuous("",breaks=c(1,2,3,4,5,6),labels=c("","","","","",""))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

#2.9# Get coefficients from linear models and plot for TL

tlData<-portfolio[,c("TL","Depth","DepthCV","Chla","ChlaCV","SBT","SBTrange")]

varlist<-names(tlData)[2:7]

tlCoefs<-matrix(data=NA,nrow=6,ncol=3)

tlCoefs[,1]<-unlist(lapply(varlist,function (x) {
  coef(lm(substitute(TL~MyStd(i),list(i=as.name(x))),data=tlData))[2]
}))

tlCoefs[,2]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(TL~MyStd(i),list(i=as.name(x))),data=tlData))[2]
}))

tlCoefs[,3]<-unlist(lapply(varlist,function (x) {
  confint(lm(substitute(TL~MyStd(i),list(i=as.name(x))),data=tlData))[4]
}))

tlCoefs<-as.data.frame(tlCoefs)
colnames(tlCoefs)<-c("coef","coefLow","coefUpp")
tlCoefs$Predictor<-varlist

c6<-ggplot(data=tlCoefs,aes(x=coef,y=seq(1,6,1)))+
  geom_vline(xintercept=0,col="gray",linetype="dashed")+
  geom_point(size=4,col=c("black","black","black","black","red","black"))+
  geom_pointrange(aes(y=seq(1,6,1),xmin=coefLow,xmax=coefUpp),
                  col=c("black","black","black","black","red","black"))+
  scale_x_continuous("Standadized effect")+
  scale_y_continuous("",breaks=c(1,2,3,4,5,6),labels=c("","","","","",""))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#2.10# All plots together

pdf(file="Abund_Div_Env_coefs.pdf",width=12,height=8)
c1+c2+c3+c4+c5+c6+plot_layout(ncol=3)
dev.off()

#2.11# Plot the most relevant relationships for the main text

p1<-ggplot(data=portfolio,aes(x=Depth,y=log(cpue)))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous(expression(paste("ln-Abundance (ind × ",km^-2,")")))+
  scale_x_continuous("Depth (m)")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position=c(0.95,0.05),legend.justification=c("right","bottom"),
        axis.text=element_text(size=14),axis.title=element_text(size=16),
        legend.text=element_text(size=8),legend.title=element_text(size=9),
        legend.key.spacing.y=unit(0.01,"pt"))+
  ggtitle("(A)")

p2<-ggplot(data=portfolio,aes(x=SBT,y=richness))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Richness")+
  scale_x_continuous("")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

p3<-ggplot(data=portfolio,aes(x=SBT,y=evenness))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Evenness")+
  scale_x_continuous("")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

p4<-ggplot(data=portfolio,aes(x=SBT,y=K))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous(expression(CWM[K]))+
  scale_x_continuous("SBT (ºC)")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

p5<-ggplot(data=portfolio,aes(x=SBT,y=L50))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous(expression(CWM[L50]))+
  scale_x_continuous("SBT (ºC)")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

p6<-ggplot(data=portfolio,aes(x=SBT,y=TL))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous(expression(CWM[TL]))+
  scale_x_continuous("SBT (ºC)")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

pdf(file="Abund_Div_Env_xy.pdf",width=12,height=8)
p1+p2+p3+p4+p5+p6+plot_layout(ncol=3)
dev.off()

#2.12# Summaries from the linear models and residuals check

#2.12.1# Run models

lm1<-lm(log(cpue)~Depth,data=portfolio)
summary(lm1)

lm2<-lm(richness~SBT,data=portfolio)
summary(lm2)

lm3<-lm(evenness~SBT,data=portfolio)
summary(lm3)

lm4<-lm(K~SBT,data=portfolio)
summary(lm4)

lm5<-lm(L50~SBT,data=portfolio)
summary(lm5)

lm6<-lm(TL~SBT,data=portfolio)
summary(lm6)

#2.12.2# Get residuals and fitted values

portfolio$residLM1<-residuals(lm1) # cpue~depth
portfolio$residLM2<-residuals(lm2) # richness~sbt
portfolio$residLM3<-residuals(lm3) # evenness~sbt
portfolio$residLM4<-residuals(lm4) # K~sbt
portfolio$residLM5<-residuals(lm5) # L50~sbt
portfolio$residLM6<-residuals(lm6) # TL~sbt

portfolio$fitLM1<-fitted(lm1)
portfolio$fitLM2<-fitted(lm2)
portfolio$fitLM3<-fitted(lm3)
portfolio$fitLM4<-fitted(lm4)
portfolio$fitLM5<-fitted(lm5)
portfolio$fitLM6<-fitted(lm6)

#2.12.3# Residuals plots

reLM1<-ggplot(data=portfolio,aes(x=fitLM1,y=residLM1))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM2<-ggplot(data=portfolio,aes(x=residLM1))+
  geom_histogram(binwidth=1)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM3<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM1),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM4<-ggplot(data=portfolio,aes(x=fitLM2,y=residLM2))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM5<-ggplot(data=portfolio,aes(x=residLM2))+
  geom_histogram(binwidth=5)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM6<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM2),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM7<-ggplot(data=portfolio,aes(x=fitLM3,y=residLM3))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM8<-ggplot(data=portfolio,aes(x=residLM3))+
  geom_histogram(binwidth=0.05)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM9<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM3),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM10<-ggplot(data=portfolio,aes(x=fitLM4,y=residLM4))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM11<-ggplot(data=portfolio,aes(x=residLM4))+
  geom_histogram(binwidth=0.05)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM12<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM4),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM13<-ggplot(data=portfolio,aes(x=fitLM5,y=residLM5))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM14<-ggplot(data=portfolio,aes(x=residLM5))+
  geom_histogram(binwidth=5)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM15<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM5),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM16<-ggplot(data=portfolio,aes(x=fitLM6,y=residLM6))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="Fitted")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM17<-ggplot(data=portfolio,aes(x=residLM6))+
  geom_histogram(binwidth=0.1)+
  labs(y="Count",x="Residuals")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

reLM18<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residLM6),size=3)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="Longitude (ºW-E)")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

pdf(file="resid_LMs.pdf",width=12,height=14)
(reLM1+reLM2+reLM3)/(reLM4+reLM5+reLM6)/(reLM7+reLM8+reLM9)/
  (reLM10+reLM11+reLM12)/(reLM13+reLM14+reLM15)/(reLM16+reLM17+reLM18)
dev.off()


# ------------------------------------------------------ #
# ---------------------- STEP 3 ------------------------ #
# ------------------------------------------------------ #
#3# Q2: Does Biodiversity relate to Ecosystem Functioning across ecosystems?
# Study BEF relationships using time series data for each survey

head(time.series.df)

summary(time.series.df[,5:8])

par(mfrow=c(2,3))
hist(time.series.df$abundanceTot)
hist(time.series.df$richness);abline(v=30,col=2) # roughly mean value
hist(time.series.df$evenness);abline(v=0.55,col=2) # roughly mean value
hist(log(time.series.df$abundanceTot)) # better distribution
hist(log(time.series.df$richness)) # not so much of an improvement
hist(log(time.series.df$evenness)) # worse distribution

time.series.df$logAbundTot<-log(time.series.df$abundanceTot) # log-transform abundance for a better distribution

#3.1# Model for Richness (centered at 30 species)

lmeRich4<-lme(logAbundTot~I(richness-30), # random intercept and slope, corAR1 because time series data and varExp for heteroskedasticity
              random=~I(richness-30)|area,
              correlation=corAR1(form=~1|area),
              weights=varExp(),
              data=time.series.df)

summary(lmeRich4)
intervals(lmeRich4)

coefRich<-as.data.frame(coef(lmeRich4))
coefRich$Area<-row.names(coefRich)
colnames(coefRich)<-c("InterceptRich","SlopeRich","Area")

#3.2# Model for Evenness (centered at 0.55 evenness units)

lmeEven5<-lme(logAbundTot~I(evenness-0.55),
              random=list(area=pdDiag(~I(evenness-0.55))),
              correlation=corAR1(form=~1|area),
              weights=varExp(),
              data=time.series.df)

summary(lmeEven5)
intervals(lmeEven5)

coefEven<-as.data.frame(coef(lmeEven5))
coefEven$Area<-row.names(coefEven)
colnames(coefEven)<-c("InterceptEven","SlopeEven","Area")

#3.3# Merge both LMEs random effects and with portfolio data

allDataAreas<-merge(coefRich,coefEven,by.x="Area",by.y="Area")

head(allDataAreas)
head(portfolio)

portfolio<-merge(portfolio,allDataAreas,by.x="System",by.y="Area",all=T)

#3.4# Residuals plots

time.series.df$residRichLME<-resid(lmeRich4,type="n")
time.series.df$fitRichLME<-fitted(lmeRich4) # predictions at the area level

time.series.df$residEvenLME<-resid(lmeEven5,type="n")
time.series.df$fitEvenLME<-fitted(lmeEven5) # predictions at the area level

reLME1<-ggplot(data=time.series.df,aes(x=fitRichLME,y=residRichLME))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

reLME2<-ggplot(data=time.series.df,aes(x=residRichLME))+
  geom_histogram()+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

reLME3<-ggplot(data=portfolio,aes(x=InterceptRich))+
  geom_histogram(bins=10)+
  labs(y="",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

reLME4<-ggplot(data=portfolio,aes(x=SlopeRich))+
  geom_histogram(bins=10)+
  labs(y="",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

reLME5<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=InterceptRich),size=2)+
  scale_colour_gradient2("Random \nintercept",low=("blue"),
                         mid=("white"),high=("red"),midpoint=9)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

reLME6<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=SlopeRich),size=2)+
  scale_colour_gradient2("Random \nslope",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0.045)+
  labs(y="",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

reLME7<-ggplot(data=time.series.df,aes(x=fitEvenLME,y=residEvenLME))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="Fitted")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

reLME8<-ggplot(data=time.series.df,aes(x=residEvenLME))+
  geom_histogram()+
  labs(y="Count",x="Residuals")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

reLME9<-ggplot(data=portfolio,aes(x=InterceptEven))+
  geom_histogram(bins=10)+
  labs(y="",x="Random intercept")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

reLME10<-ggplot(data=portfolio,aes(x=SlopeEven))+
  geom_histogram(bins=10)+
  labs(y="",x="Random slope")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

reLME11<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=InterceptEven),size=2)+
  scale_colour_gradient2("Random \nintercept",low=("blue"),
                         mid=("white"),high=("red"),midpoint=9)+
  labs(y="Latitude (ºN)",x="Longitude (ºW-E)")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

reLME12<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=SlopeEven),size=2)+
  scale_colour_gradient2("Random \nslope",low=("blue"),
                         mid=("white"),high=("red"),midpoint=-3.2)+
  labs(y="",x="Longitude (ºW-E)")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

pdf(file="resid_LMMs.pdf",width=22,height=6)
reLME1+reLME2+reLME3+reLME4+reLME5+reLME6+
  reLME7+reLME8+reLME9+reLME10+reLME11+reLME12+plot_layout(ncol=6)
dev.off()

#3.5# Plot predicted results from the LMEs

time.series.df$pred_popRich<-predict(lmeRich4,level=0) # get predictions at the population level
time.series.df$pred_popEven<-predict(lmeEven5,level=0)

plme1<-ggplot(data=time.series.df,aes(x=richness,y=fitRichLME,group=area))+
  geom_point(data=time.series.df,aes(x=richness,y=logAbundTot),
             alpha=0.4,col="gray60",size=0.1)+
  geom_line(color="gray50")+
  geom_line(aes(x=richness,y=pred_popRich),color="black",linewidth=2)+
  scale_y_continuous(expression(paste("ln-Abundance (ind × ",km^-2,")")))+
  scale_x_continuous("Richness",breaks=c(1,25,50,75,100,125))+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

plme2<-ggplot(data=time.series.df,aes(x=evenness,y=fitEvenLME,group=area))+
  geom_point(data=time.series.df,aes(x=evenness,y=logAbundTot),
             alpha=0.4,col="gray60",size=0.1)+
  geom_line(color="gray50")+
  geom_line(aes(x=evenness,y=pred_popEven),color="black",linewidth=2)+
  scale_y_continuous(expression(paste("ln-Abundance (ind × ",km^-2,")")))+
  scale_x_continuous("Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

#3.6# Explore relationships between random slopes and
# community properties and environment and plot relevant relationships

#3.6.1# All potential correlations

names(portfolio)

portfolio %>% 
  dplyr::select(c(7:18,20,22,23,26,27,30,46,48)) %>% 
  correlate(use="pairwise.complete.obs") %>% 
  focus(SlopeRich,SlopeEven)

#3.6.2# Relationships for richness slopes

mr1<-ggplot(data=portfolio,aes(x=log(cpue),y=SlopeRich))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Effect of richness \non abundance")+
  scale_x_continuous("ln-Abundance")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  guides(col=guide_legend(ncol=2,override.aes=list(size=1)))+
  ggtitle("(B)")

mr2<-ggplot(data=portfolio,aes(x=SBTrange,y=SlopeRich))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Effect of richness \non abundance")+
  scale_x_continuous("SBTrange")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

mr3<-ggplot(data=portfolio,aes(x=Chla,y=SlopeRich))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("")+
  scale_x_continuous("Chla")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

mr4<-ggplot(data=portfolio,aes(x=Depth,y=SlopeRich))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Effect of richness \non abundance")+
  scale_x_continuous("Depth")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

mr5<-ggplot(data=portfolio,aes(x=Trawl,y=SlopeRich))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("")+
  scale_x_continuous("Trawl")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#3.6.3# Relationships for evenness slopes

mr6<-ggplot(data=portfolio,aes(x=synchrony,y=SlopeEven))+
  geom_point(aes(col=Region),size=3)+
  geom_smooth(method="lm",se=F,col="gray50")+
  scale_y_continuous("Effect of evenness \non abundance")+
  scale_x_continuous("Synchrony")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

pdf(file="metareg_LMMs.pdf",width=8,height=12)
plme1+mr1+mr2+mr3+mr4+mr5+plme2+mr6+plot_layout(ncol=2)
dev.off()


# ------------------------------------------------------ #
# ---------------------- STEP 4 ------------------------ #
# ------------------------------------------------------ #
#4# Q3: What are the direct and indirect drivers of the Portfolio Effect
# across the studied ecosystems? Does it change with the type of PE metrics? 
# Use a SEM for MV_PE and then for AvCV_PE
# From Theory and Experiments: MV_PE may be affected by all variables
# (the four community metrics, environment and human impacts). Z-value may be
# affected by synchrony (see Fig. S12f in Anderson et al.) and richness
# (see Fig. S12e in Anderson et al.). Synchrony may be affected by evenness
# (see Fig. 5 in Anderson et al.) and richness (theoretical).
# Finally, evenness may be related to richness, and richness related to the environment.

#4.1# Plot distribution of MV_PE, AvCV_PE and z-value

names(portfolio)

portfolio<-portfolio[order(portfolio$plotOrd),] # order data for plotting purposes

pt1<-ggplot(data=portfolio,aes(x=z.value,xmin=z.value.l,xmax=z.value.u,y=plotOrd))+
  geom_point(aes(col=Region),size=2.5)+
  geom_pointrange(aes(col=Region),position=position_dodge(width=0.25))+
  geom_vline(xintercept=2)+
  scale_x_continuous("Taylor's power law z-value",breaks=c(1.7,2,2.3,2.6))+
  scale_y_continuous("",breaks=seq(1,54,1),labels=portfolio$System)+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

pt2<-ggplot(data=portfolio,aes(x=mvpe,xmin=mvpe.l,xmax=mvpe.u,y=plotOrd))+
  geom_point(aes(col=Region),size=2.5)+
  geom_pointrange(aes(col=Region),position=position_dodge(width=0.25))+
  geom_vline(xintercept=1)+
  scale_x_continuous("Portfolio Effect (MVPE)")+
  scale_y_continuous("",breaks=seq(1,54,1),labels=portfolio$System)+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",axis.text.y=element_blank(),
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

pt3<-ggplot(data=portfolio,aes(x=avcvpe,xmin=avcvpe.l,xmax=avcvpe.u,y=plotOrd))+
  geom_point(aes(col=Region),size=2.5)+
  geom_pointrange(aes(col=Region),position=position_dodge(width=0.25))+
  geom_vline(xintercept=1)+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),"Portfolio Effect (CVPE)")+
  scale_y_continuous("",breaks=seq(1,54,1),labels=portfolio$System)+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",axis.text.y=element_blank(),
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

pdf(file="PEandZ_values.pdf",width=12,height=6)
pt1+pt2+pt3+plot_layout(ncol=3)
dev.off()

#4.2# Explore correlations for each of the key variables to later build the SEM
# Start with theoretical complete models and run a stepAIC

names(portfolio)

m1.step<-lm(mvpe~z.value+synchrony+evenness+richness+
              SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m1.step)

m2.step<-lm(z.value~synchrony+evenness+richness+
              SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m2.step)

m3.step<-lm(synchrony~evenness+richness+
              SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m3.step)

m4.step<-lm(evenness~richness+
              SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m4.step)

m5.step<-lm(richness~SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m5.step)

#4.3# Standardize variables before fitting the SEM

portfolio$mvpe<-MyStd(portfolio$mvpe)
portfolio$avcvpe<-MyStd(portfolio$avcvpe)
portfolio$z.value<-MyStd(portfolio$z.value)
portfolio$synchrony<-MyStd(portfolio$synchrony)
portfolio$evenness<-MyStd(portfolio$evenness)
portfolio$richness<-MyStd(portfolio$richness)
portfolio$Depth<-MyStd(portfolio$Depth)
portfolio$DepthCV<-MyStd(portfolio$DepthCV)
portfolio$Chla<-MyStd(portfolio$Chla)
portfolio$ChlaCV<-MyStd(portfolio$ChlaCV)
portfolio$SBT<-MyStd(portfolio$SBT)
portfolio$SBTrange<-MyStd(portfolio$SBTrange)
portfolio$TAI<-MyStd(portfolio$TAI)
portfolio$Trawl<-MyStd(portfolio$Trawl)

#4.4# Fit SEM to MV_PE using lavaan

peSEM.mv<-'mvpe      ~ z.value + synchrony
           z.value   ~ evenness + richness + SBT + Depth
           synchrony ~ evenness + richness + SBTrange + Chla + TAI
           evenness  ~ SBTrange + Chla + DepthCV + Trawl
           richness  ~ Chla + ChlaCV + DepthCV'

defSEM<-sem(peSEM.mv,data=portfolio)
summary(defSEM,rsquare=T,standardized=T)
lavaanPlot(model=defSEM,coefs=T,stand=T,graph_options=list(layout="circo"),sig=0.05)

#4.5# Fit SEM using picewiseSEM

psemMVPE1<-lm(mvpe      ~ synchrony + z.value,data=portfolio)
psemMVPE2<-lm(z.value   ~ evenness + richness + SBT + Depth,data=portfolio)
psemMVPE3<-lm(synchrony ~ evenness + richness + SBTrange + Chla + TAI,data=portfolio)
psemMVPE4<-lm(evenness  ~ SBTrange + Chla + DepthCV + Trawl,data=portfolio)
psemMVPE5<-lm(richness  ~ Chla + ChlaCV + DepthCV,data=portfolio)

psemMVPE.def<-psem(psemMVPE1,psemMVPE2,psemMVPE3,psemMVPE4,psemMVPE5)

coefs(psemMVPE.def,standardize="scale")
rsquared(psemMVPE.def)
dSep(psemMVPE.def)
fisherC(psemMVPE.def)
plot(psemMVPE.def)

vif(psemMVPE1);vif(psemMVPE2);vif(psemMVPE3);vif(psemMVPE4);vif(psemMVPE5)

#4.6# Calculating direct and indirect effects for MV_PE using lavaan syntax

#4.6.1# Calculation within the model formulation

modelEffMVPE<-'mvpe      ~ mz*z.value + ms*synchrony
               z.value   ~ ze*evenness + zr*richness + zs*SBT + zd*Depth 
               synchrony ~ se*evenness + sr*richness + ss*SBTrange + sc*Chla + st*TAI 
               evenness  ~ es*SBTrange + ec*Chla + ed*DepthCV + et*Trawl
               richness  ~ ra*Chla + rc*ChlaCV + rd*DepthCV
           
               directZ          := mz
               directSyn        := ms
               
               indirectRich     := (zr*mz) + (sr*ms)
               indirectEven     := (ze*mz) + (se*ms)
               
               indirectSBT      := zs*mz
               indirectSBTrange := (ss*ms) + (es*se*ms) + (es*ze*mz)
               
               indirectChla     := (sc*ms) + (ec*se*ms) + (ec*ze*mz) + (ra*sr*ms) + (ra*zr*mz)
               indirectChlaCV   := (rc*sr*ms) + (rc*zr*mz)
               
               indirectDepth    := zd*mz
               indirectDepthCV  := (ed*se*ms) + (ed*ze*mz) + (rd*sr*ms) + (rd*zr*mz)
               
               indirectTAI      := st*ms
               indirectTrawl    := (et*se*ms) + (et*ze*mz)
               
               totalZ           := directZ
               totalSyn         := directSyn
               totalRich        := indirectRich
               totalEven        := indirectEven
               totalSBT         := indirectSBT
               totalSBTrange    := indirectSBTrange
               totalChla        := indirectChla
               totalChlaCV      := indirectChlaCV
               totalDepth       := indirectDepth
               totalDepthCV     := indirectDepthCV
               totalTAI         := indirectTAI
               totalTrawl       := indirectTrawl'

partial.defSEM.MVPEeffects<-sem(modelEffMVPE,data=portfolio)
summary(partial.defSEM.MVPEeffects,standardized=T)

#4.6.2# Plot direct and indirect effects

sem.MVPEeffects<-data.frame(summary(partial.defSEM.MVPEeffects,standardized=T)$pe[60:83,])

sem.MVPEeffects$varName<-c("z-value","Synchrony","Richness","Evenness","SBT",
                           "SBTrange","Chla","ChlaCV","Depth","DepthCV","TAI",
                           "Trawl","z-value","Synchrony","Richness","Evenness","SBT",
                           "SBTrange","Chla","ChlaCV","Depth","DepthCV","TAI","Trawl")

sem.MVPEeffects$effType<-c("Direct","Direct","Indirect","Indirect","Indirect","Indirect",
                           "Indirect","Indirect","Indirect","Indirect","Indirect","Indirect",
                           "Total","Total","Total","Total","Total","Total",
                           "Total","Total","Total","Total","Total","Total")

sem.MVPEeffects$varNum<-factor(c(1,2,3,4,5,6,7,8,9,10,11,12,
                                 1,2,3,4,5,6,7,8,9,10,11,12))

pdf(file="MVPE_SEM_dir_indir.pdf")

ggplot(data=filter(sem.MVPEeffects,effType != "Total"),aes(y=varNum,x=std.all))+
  geom_col(aes(fill=effType))+
  #geom_point(data=filter(sem.MVPEeffects,effType == "Total"),aes(y=varNum,x=std.all))+ # net effects
  #geom_point(data=filter(sem.MVPEeffects,effType == "Total" & pvalue < 0.05),
  #           aes(y=varNum,x=std.all),size=3,shape=8)+ # statistically significant effects
  #geom_point(data=filter(sem.MVPEeffects,effType == "Total" & pvalue > 0.05),
  #           aes(y=varNum,x=std.all),size=3,shape=4)+ # statistically non-significant effects
  scale_y_discrete("",labels=c("1"="z-value","2"="Synchrony","3"="Richness",
                               "4"="Evenness","5"="SBT","6"="SBTrange",
                               "7"="Chla","8"="ChlaCV","9"="Depth",
                               "10"="DepthCV","11"="TAI","12"="Trawl"))+
  scale_x_continuous("Standardized effect on MVPE",limits=c(-0.75,0.5))+
  scale_fill_manual("",values=c("#E7B800","#00AFBB"))+
  geom_vline(xintercept=0,linetype="dashed")+
  theme_classic()+
  background_grid()+
  theme(legend.position=c(0,1),legend.justification=c(-0.4,1.2),
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

dev.off()

#4.7# Residuals for each model included in the psem

portfolio$residMVPE1<-residuals(psemMVPE1) # MVPE model
portfolio$residMVPE2<-residuals(psemMVPE2) # z-value model
portfolio$residMVPE3<-residuals(psemMVPE3) # Synchrony model
portfolio$residMVPE4<-residuals(psemMVPE4) # Evenness model
portfolio$residMVPE5<-residuals(psemMVPE5) # Richness model

portfolio$fitMVPE1<-fitted(psemMVPE1)
portfolio$fitMVPE2<-fitted(psemMVPE2)
portfolio$fitMVPE3<-fitted(psemMVPE3)
portfolio$fitMVPE4<-fitted(psemMVPE4)
portfolio$fitMVPE5<-fitted(psemMVPE5)

reMVPE1<-ggplot(data=portfolio,aes(x=fitMVPE1,y=residMVPE1))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

reMVPE2<-ggplot(data=portfolio,aes(x=residMVPE1))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

reMVPE3<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residMVPE1),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

reMVPE4<-ggplot(data=portfolio,aes(x=fitMVPE2,y=residMVPE2))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

reMVPE5<-ggplot(data=portfolio,aes(x=residMVPE2))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

reMVPE6<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residMVPE2),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

reMVPE7<-ggplot(data=portfolio,aes(x=fitMVPE3,y=residMVPE3))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

reMVPE8<-ggplot(data=portfolio,aes(x=residMVPE3))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

reMVPE9<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residMVPE3),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

reMVPE10<-ggplot(data=portfolio,aes(x=fitMVPE4,y=residMVPE4))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

reMVPE11<-ggplot(data=portfolio,aes(x=residMVPE4))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

reMVPE12<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residMVPE4),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

reMVPE13<-ggplot(data=portfolio,aes(x=fitMVPE5,y=residMVPE5))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="Fitted")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(M)")

reMVPE14<-ggplot(data=portfolio,aes(x=residMVPE5))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="Residuals")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(N)")

reMVPE15<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residMVPE5),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="Longitude (ºW-E)")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(O)")

pdf(file="resids_MVPE.pdf",width=12,height=15)
(reMVPE1+reMVPE2+reMVPE3)/(reMVPE4+reMVPE5+reMVPE6)/
  (reMVPE7+reMVPE8+reMVPE9)/(reMVPE10+reMVPE11+reMVPE12)/
  (reMVPE13+reMVPE14+reMVPE15)
dev.off()

#4.8# Partial plots for MV_PE SEM

coef(defSEM)

#4.8.1# Partial z-value effect on MV_PE

portfolio$MV_PEpart_z<- portfolio$mvpe -
  (coef(defSEM)[2]*portfolio$synchrony)

x.line.pe.z<-seq(from=min(portfolio$z.value),to=max(portfolio$z.value),by=0.1)
y.line.pe.z<- coef(defSEM)[1]*x.line.pe.z
line.pe.z<-as.data.frame(cbind(x.line.pe.z,y.line.pe.z))

#4.8.2# Partial synchrony effect on MV_PE

portfolio$MV_PEpart_syn<- portfolio$mvpe -
  (coef(defSEM)[1]*portfolio$z.value)

x.line.pe.syn<-seq(from=min(portfolio$synchrony),to=max(portfolio$synchrony),by=0.1)
y.line.pe.syn<- coef(defSEM)[2]*x.line.pe.syn
line.pe.syn<-as.data.frame(cbind(x.line.pe.syn,y.line.pe.syn))

#4.8.3# Partial evenness effect on z-value

portfolio$Zpart_even<- portfolio$z.value -
  (coef(defSEM)[4]*portfolio$richness) -
  (coef(defSEM)[5]*portfolio$SBT) - 
  (coef(defSEM)[6]*portfolio$Depth)

x.line.z.even<-seq(from=min(portfolio$evenness),to=max(portfolio$evenness),by=0.1)
y.line.z.even<-coef(defSEM)[3]*x.line.z.even
line.z.even<-as.data.frame(cbind(x.line.z.even,y.line.z.even))

#4.8.4# Partial richness effect on z-value

portfolio$Zpart_rich<- portfolio$z.value -
  (coef(defSEM)[3]*portfolio$evenness) -
  (coef(defSEM)[5]*portfolio$SBT) - 
  (coef(defSEM)[6]*portfolio$Depth)

x.line.z.rich<-seq(from=min(portfolio$richness),to=max(portfolio$richness),by=0.1)
y.line.z.rich<-coef(defSEM)[4]*x.line.z.rich
line.z.rich<-as.data.frame(cbind(x.line.z.rich,y.line.z.rich))

#4.8.5# Partial SBT effect on z-value

portfolio$Zpart_sbt<- portfolio$z.value -
  (coef(defSEM)[3]*portfolio$evenness) -
  (coef(defSEM)[4]*portfolio$richness) -
  (coef(defSEM)[6]*portfolio$Depth)

x.line.z.sbt<-seq(from=min(portfolio$SBT),to=max(portfolio$SBT),by=0.1)
y.line.z.sbt<-coef(defSEM)[5]*x.line.z.sbt
line.z.sbt<-as.data.frame(cbind(x.line.z.sbt,y.line.z.sbt))

#4.8.6# Partial Depth effect on z-value

portfolio$Zpart_depth<- portfolio$z.value -
  (coef(defSEM)[3]*portfolio$evenness) -
  (coef(defSEM)[4]*portfolio$richness) -
  (coef(defSEM)[5]*portfolio$SBT)

x.line.z.depth<-seq(from=min(portfolio$Depth),to=max(portfolio$Depth),by=0.1)
y.line.z.depth<-coef(defSEM)[6]*x.line.z.depth
line.z.depth<-as.data.frame(cbind(x.line.z.depth,y.line.z.depth))

#4.8.7# Partial evenness effect on synchrony

portfolio$SYNpart_even<- portfolio$synchrony -
  (coef(defSEM)[8]*portfolio$richness) -
  (coef(defSEM)[9]*portfolio$SBTrange) - 
  (coef(defSEM)[10]*portfolio$Chla) - 
  (coef(defSEM)[11]*portfolio$TAI)

x.line.syn.even<-seq(from=min(portfolio$evenness),to=max(portfolio$evenness),by=0.1)
y.line.syn.even<-coef(defSEM)[7]*x.line.syn.even
line.syn.even<-as.data.frame(cbind(x.line.syn.even,y.line.syn.even))

#4.8.8# Partial richness effect on synchrony

portfolio$SYNpart_rich<- portfolio$synchrony -
  (coef(defSEM)[7]*portfolio$evenness) -
  (coef(defSEM)[9]*portfolio$SBTrange) - 
  (coef(defSEM)[10]*portfolio$Chla) - 
  (coef(defSEM)[11]*portfolio$TAI)

x.line.syn.rich<-seq(from=min(portfolio$richness),to=max(portfolio$richness),by=0.1)
y.line.syn.rich<-coef(defSEM)[8]*x.line.syn.rich
line.syn.rich<-as.data.frame(cbind(x.line.syn.rich,y.line.syn.rich))

#4.8.9# Partial SBTrange effect on synchrony

portfolio$SYNpart_sbt<- portfolio$synchrony -
  (coef(defSEM)[7]*portfolio$evenness) -
  (coef(defSEM)[8]*portfolio$richness) - 
  (coef(defSEM)[10]*portfolio$Chla) - 
  (coef(defSEM)[11]*portfolio$TAI)

x.line.syn.sbt<-seq(from=min(portfolio$SBTrange),to=max(portfolio$SBTrange),by=0.1)
y.line.syn.sbt<-coef(defSEM)[9]*x.line.syn.sbt
line.syn.sbt<-as.data.frame(cbind(x.line.syn.sbt,y.line.syn.sbt))

#4.8.10# Partial Chla effect on synchrony

portfolio$SYNpart_chla<- portfolio$synchrony -
  (coef(defSEM)[7]*portfolio$evenness) -
  (coef(defSEM)[8]*portfolio$richness) - 
  (coef(defSEM)[9]*portfolio$SBTrange) - 
  (coef(defSEM)[11]*portfolio$TAI)

x.line.syn.chla<-seq(from=min(portfolio$Chla),to=max(portfolio$Chla),by=0.1)
y.line.syn.chla<-coef(defSEM)[10]*x.line.syn.chla
line.syn.chla<-as.data.frame(cbind(x.line.syn.chla,y.line.syn.chla))

#4.8.11# Partial TAI effect on synchrony

portfolio$SYNpart_tai<- portfolio$synchrony -
  (coef(defSEM)[7]*portfolio$evenness) -
  (coef(defSEM)[8]*portfolio$richness) - 
  (coef(defSEM)[9]*portfolio$SBTrange) - 
  (coef(defSEM)[10]*portfolio$Chla)

x.line.syn.tai<-seq(from=min(portfolio$TAI),to=max(portfolio$TAI),by=0.1)
y.line.syn.tai<-coef(defSEM)[11]*x.line.syn.tai
line.syn.tai<-as.data.frame(cbind(x.line.syn.tai,y.line.syn.tai))

#4.8.12# Partial SBTrange effect on evenness

portfolio$EVENpart_sbt<- portfolio$evenness -
  (coef(defSEM)[13]*portfolio$Chla) -
  (coef(defSEM)[14]*portfolio$DepthCV) -
  (coef(defSEM)[15]*portfolio$Trawl)

x.line.even.sbt<-seq(from=min(portfolio$SBTrange),to=max(portfolio$SBTrange),by=0.1)
y.line.even.sbt<-coef(defSEM)[12]*x.line.even.sbt
line.even.sbt<-as.data.frame(cbind(x.line.even.sbt,y.line.even.sbt))

#4.8.13# Partial Chla effect on evenness

portfolio$EVENpart_chla<- portfolio$evenness -
  (coef(defSEM)[12]*portfolio$SBTrange) -
  (coef(defSEM)[14]*portfolio$DepthCV) -
  (coef(defSEM)[15]*portfolio$Trawl)

x.line.even.chla<-seq(from=min(portfolio$Chla),to=max(portfolio$Chla),by=0.1)
y.line.even.chla<-coef(defSEM)[13]*x.line.even.chla
line.even.chla<-as.data.frame(cbind(x.line.even.chla,y.line.even.chla))

#4.8.14# Partial DepthCV effect on evenness

portfolio$EVENpart_depth<- portfolio$evenness -
  (coef(defSEM)[12]*portfolio$SBTrange) -
  (coef(defSEM)[13]*portfolio$Chla) -
  (coef(defSEM)[15]*portfolio$Trawl)

x.line.even.depth<-seq(from=min(portfolio$DepthCV),to=max(portfolio$DepthCV),by=0.1)
y.line.even.depth<-coef(defSEM)[14]*x.line.even.depth
line.even.depth<-as.data.frame(cbind(x.line.even.depth,y.line.even.depth))

#4.8.15# Partial Trawl effect on evenness

portfolio$EVENpart_trawl<- portfolio$evenness -
  (coef(defSEM)[12]*portfolio$SBTrange) -
  (coef(defSEM)[13]*portfolio$Chla) -
  (coef(defSEM)[14]*portfolio$DepthCV)

x.line.even.trawl<-seq(from=min(portfolio$Trawl),to=max(portfolio$Trawl),by=0.1)
y.line.even.trawl<-coef(defSEM)[15]*x.line.even.trawl
line.even.trawl<-as.data.frame(cbind(x.line.even.trawl,y.line.even.trawl))

#4.8.16# Partial Chla effect on richness

portfolio$RICHpart_chla<- portfolio$richness -
  (coef(defSEM)[17]*portfolio$ChlaCV) - 
  (coef(defSEM)[18]*portfolio$DepthCV)

x.line.rich.chla<-seq(from=min(portfolio$Chla),to=max(portfolio$Chla),by=0.1)
y.line.rich.chla<-coef(defSEM)[16]*x.line.rich.chla
line.rich.chla<-as.data.frame(cbind(x.line.rich.chla,y.line.rich.chla))

#4.8.17# Partial ChlaCV effect on richness

portfolio$RICHpart_chlacv<- portfolio$richness -
  (coef(defSEM)[16]*portfolio$Chla) -
  (coef(defSEM)[18]*portfolio$DepthCV)

x.line.rich.chlacv<-seq(from=min(portfolio$ChlaCV),to=max(portfolio$ChlaCV),by=0.1)
y.line.rich.chlacv<-coef(defSEM)[17]*x.line.rich.chlacv
line.rich.chlacv<-as.data.frame(cbind(x.line.rich.chlacv,y.line.rich.chlacv))

#4.8.18# Partial DepthCV effect on richness

portfolio$RICHpart_depth<- portfolio$richness -
  (coef(defSEM)[16]*portfolio$Chla) - 
  (coef(defSEM)[17]*portfolio$ChlaCV)

x.line.rich.depth<-seq(from=min(portfolio$DepthCV),to=max(portfolio$DepthCV),by=0.1)
y.line.rich.depth<-coef(defSEM)[18]*x.line.rich.depth
line.rich.depth<-as.data.frame(cbind(x.line.rich.depth,y.line.rich.depth))

#4.8.19# All partial effects plots together

names(portfolio)

pe1<-ggplot()+
  geom_point(data=portfolio,aes(y=MV_PEpart_z,x=z.value,col=Region))+
  geom_line(data=line.pe.z,aes(x=x.line.pe.z,y=y.line.pe.z),linewidth=1,col="gray50")+
  labs(y="MVPE",x="z-value")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

pe2<-ggplot()+
  geom_point(data=portfolio,aes(y=MV_PEpart_syn,x=synchrony,col=Region))+
  geom_line(data=line.pe.syn,aes(x=x.line.pe.syn,y=y.line.pe.syn),linewidth=1,col="gray50")+
  labs(y="",x="Synchrony")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

pe3<-ggplot()+
  geom_point(data=portfolio,aes(y=Zpart_even,x=evenness,col=Region))+
  geom_line(data=line.z.even,aes(x=x.line.z.even,y=y.line.z.even),linewidth=1,col="gray50")+
  labs(y="z-value",x="Evenness")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

pe4<-ggplot()+
  geom_point(data=portfolio,aes(y=Zpart_rich,x=richness,col=Region))+
  geom_line(data=line.z.rich,aes(x=x.line.z.rich,y=y.line.z.rich),linewidth=1,col="gray50")+
  labs(y="",x="Richness")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

pe5<-ggplot()+
  geom_point(data=portfolio,aes(y=Zpart_sbt,x=SBT,col=Region))+
  geom_line(data=line.z.sbt,aes(x=x.line.z.sbt,y=y.line.z.sbt),linewidth=1,col="gray50")+
  labs(y="",x="SBT")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

pe6<-ggplot()+
  geom_point(data=portfolio,aes(y=Zpart_depth,x=Depth,col=Region))+
  geom_line(data=line.z.depth,aes(x=x.line.z.depth,y=y.line.z.depth),linewidth=1,col="gray50")+
  labs(y="",x="Depth")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

pe7<-ggplot()+
  geom_point(data=portfolio,aes(y=SYNpart_even,x=evenness,col=Region))+
  geom_line(data=line.syn.even,aes(x=x.line.syn.even,y=y.line.syn.even),linewidth=1,col="gray50")+
  labs(y="Synchrony",x="Evenness")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

pe8<-ggplot()+
  geom_point(data=portfolio,aes(y=SYNpart_rich,x=richness,col=Region))+
  geom_line(data=line.syn.rich,aes(x=x.line.syn.rich,y=y.line.syn.rich),linewidth=1,col="gray50")+
  labs(y="",x="Richness")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

pe9<-ggplot()+
  geom_point(data=portfolio,aes(y=SYNpart_sbt,x=SBTrange,col=Region))+
  geom_line(data=line.syn.sbt,aes(x=x.line.syn.sbt,y=y.line.syn.sbt),linewidth=1,col="gray50")+
  labs(y="",x="SBTrange")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

pe10<-ggplot()+
  geom_point(data=portfolio,aes(y=SYNpart_chla,x=Chla,col=Region))+
  geom_line(data=line.syn.chla,aes(x=x.line.syn.chla,y=y.line.syn.chla),linewidth=1,col="gray50")+
  labs(y="",x="Chla")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

pe11<-ggplot()+
  geom_point(data=portfolio,aes(y=SYNpart_tai,x=TAI,col=Region))+
  geom_line(data=line.syn.tai,aes(x=x.line.syn.tai,y=y.line.syn.tai),linewidth=1,col="gray50")+
  labs(y="",x="TAI")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

pe12<-ggplot()+
  geom_point(data=portfolio,aes(y=EVENpart_sbt,x=SBTrange,col=Region))+
  geom_line(data=line.even.sbt,aes(x=x.line.even.sbt,y=y.line.even.sbt),linewidth=1,col="gray50")+
  labs(y="Evenness",x="SBTrange")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

pe13<-ggplot()+
  geom_point(data=portfolio,aes(y=EVENpart_chla,x=Chla,col=Region))+
  geom_line(data=line.even.chla,aes(x=x.line.even.chla,y=y.line.even.chla),linewidth=1,col="gray50")+
  labs(y="",x="Chla")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(M)")

pe14<-ggplot()+
  geom_point(data=portfolio,aes(y=EVENpart_depth,x=DepthCV,col=Region))+
  geom_line(data=line.even.depth,aes(x=x.line.even.depth,y=y.line.even.depth),linewidth=1,col="gray50")+
  labs(y="",x="DepthCV")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(N)")

pe15<-ggplot()+
  geom_point(data=portfolio,aes(y=EVENpart_trawl,x=Trawl,col=Region))+
  geom_line(data=line.even.trawl,aes(x=x.line.even.trawl,y=y.line.even.trawl),linewidth=1,col="gray50")+
  labs(y="",x="Trawl")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(O)")

pe16<-ggplot()+
  geom_point(data=portfolio,aes(y=RICHpart_chla,x=Chla,col=Region))+
  geom_line(data=line.rich.chla,aes(x=x.line.rich.chla,y=y.line.rich.chla),linewidth=1,col="gray50")+
  labs(y="Richness",x="Chla")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(P)")

pe17<-ggplot()+
  geom_point(data=portfolio,aes(y=RICHpart_chlacv,x=ChlaCV,col=Region))+
  geom_line(data=line.rich.chlacv,aes(x=x.line.rich.chlacv,y=y.line.rich.chlacv),linewidth=1,col="gray50")+
  labs(y="",x="ChlaCV")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(Q)")

pe18<-ggplot()+
  geom_point(data=portfolio,aes(y=RICHpart_depth,x=DepthCV,col=Region))+
  geom_line(data=line.rich.depth,aes(x=x.line.rich.depth,y=y.line.rich.depth),linewidth=1,col="gray50")+
  labs(y="",x="DepthCV")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(R)")

pdf(file="MVPE_SEM_partial_plots.pdf",width=13,height=11)
pe1+pe2+plot_spacer()+plot_spacer()+plot_spacer()+
  pe3+pe4+pe5+pe6+plot_spacer()+
  pe7+pe8+pe9+pe10+pe11+
  pe12+pe13+pe14+pe15+plot_spacer()+
  pe16+pe17+pe18+plot_spacer()+plot_spacer()+
  plot_layout(ncol=5)
dev.off()


# ------------------------------------------------------ #
# ------------------------------------------------------ #
#5# Check if community composition based on the CWMs
# influences somehow the partial relationships from the MV_PE SEM

#5.1# Establish CWM categories (below and above the mean value for each CWM)

names(portfolio)

ggcorrplot(cor(portfolio[,c(7:9)],method="pearson",
               use="pairwise.complete.obs"),
           type="lower",colors=c("blue","white","red"),
           lab=T,lab_size=5,show.legend=T,tl.cex=10)

summary(lm(K~L50,data=portfolio))
summary(lm(TL~L50,data=portfolio))
summary(lm(K~TL,data=portfolio))

tr1<-ggplot(data=portfolio,aes(x=K))+
  geom_histogram(bins=10,col="white")+
  geom_vline(aes(xintercept=mean(K)),col="red")+
  labs(y="Count",x=expression(CWM[K]))+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

tr2<-ggplot(data=portfolio,aes(x=L50))+
  geom_histogram(bins=10,col="white")+
  geom_vline(aes(xintercept=mean(L50)),col="red")+
  labs(y="",x=expression(CWM[L50]))+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

tr3<-ggplot(data=portfolio,aes(x=TL))+
  geom_histogram(bins=10,col="white")+
  geom_vline(aes(xintercept=mean(TL)),col="red")+
  labs(y="",x=expression(CWM[TL]))+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

tr4<-ggplot(data=portfolio,aes(x=L50,y=K))+
  geom_smooth(method="lm",se=F,col="black")+
  geom_point(aes(size=TL,col=Region))+
  scale_color_brewer(palette="Set1")+
  scale_y_continuous(expression(CWM[K]))+
  scale_x_continuous(expression(CWM[L50]),breaks=c(15,20,25,30,35,40,45))+
  labs(size=expression(CWM[TL]))+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

pdf(file="CWMs.pdf",width=10,height=8)
(tr1+tr2+tr3)/tr4
dev.off()

portfolio$Kgroup<-factor(ifelse(portfolio$K < mean(portfolio$K,na.rm=T),"Slow K","Fast K"))
portfolio$L50group<-factor(ifelse(portfolio$L50 < mean(portfolio$L50,na.rm=T),"Early maturity","Late maturity"))
portfolio$TLgroup<-factor(ifelse(portfolio$TL < mean(portfolio$TL,na.rm=T),"Low TL","High TL"))

#5.2# Univariate relationships for MV_PE partial effects

#5.2.1# MVPE~z.value

i1<-ggplot(data=portfolio,aes(y=MV_PEpart_z,x=z.value))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="MVPE",x="z-value")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

i2<-ggplot(data=portfolio,aes(y=MV_PEpart_z,x=z.value))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="z-value")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

i3<-ggplot(data=portfolio,aes(y=MV_PEpart_z,x=z.value))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="z-value")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#5.2.2# MVPE~synchrony

i4<-ggplot(data=portfolio,aes(y=MV_PEpart_syn,x=synchrony))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="MVPE",x="Synchrony")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

i5<-ggplot(data=portfolio,aes(y=MV_PEpart_syn,x=synchrony))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Synchrony")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

i6<-ggplot(data=portfolio,aes(y=MV_PEpart_syn,x=synchrony))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Synchrony")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

pdf(file="MVPE_cwm_int.pdf",width=10,height=5)
(i1+i2+i3)/(i4+i5+i6)
dev.off()

#5.3# Univariate relationships for Z-value partial effects

#5.3.1# z.value~evenness

i7<-ggplot(data=portfolio,aes(y=Zpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="z-value",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

i8<-ggplot(data=portfolio,aes(y=Zpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

i9<-ggplot(data=portfolio,aes(y=Zpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#5.3.2# z.value~richness

i10<-ggplot(data=portfolio,aes(y=Zpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="z-value",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

i11<-ggplot(data=portfolio,aes(y=Zpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

i12<-ggplot(data=portfolio,aes(y=Zpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#5.3.3# z.value~SBT

i13<-ggplot(data=portfolio,aes(y=Zpart_sbt,x=SBT))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="z-value",x="SBT")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

i14<-ggplot(data=portfolio,aes(y=Zpart_sbt,x=SBT))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="SBT")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

i15<-ggplot(data=portfolio,aes(y=Zpart_sbt,x=SBT))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="SBT")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

#5.3.4# z.value~Depth

i16<-ggplot(data=portfolio,aes(y=Zpart_depth,x=Depth))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="z-value",x="Depth")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

i17<-ggplot(data=portfolio,aes(y=Zpart_depth,x=Depth))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Depth")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

i18<-ggplot(data=portfolio,aes(y=Zpart_depth,x=Depth))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Depth")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

pdf(file="Z_cwm_int.pdf",width=12,height=10)
(i7+i8+i9)/(i10+i11+i12)/(i13+i14+i15)/(i16+i17+i18)
dev.off()

#5.4# Univariate relationships for Synchrony partial effects

#5.4.1# synchrony~evenness

i19<-ggplot(data=portfolio,aes(y=SYNpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Synchrony",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

i20<-ggplot(data=portfolio,aes(y=SYNpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

i21<-ggplot(data=portfolio,aes(y=SYNpart_even,x=evenness))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Evenness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#5.4.2# synchrony~richness

i22<-ggplot(data=portfolio,aes(y=SYNpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Synchrony",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

i23<-ggplot(data=portfolio,aes(y=SYNpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

i24<-ggplot(data=portfolio,aes(y=SYNpart_rich,x=richness))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Richness")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#5.4.3# synchrony~SBTrange

i25<-ggplot(data=portfolio,aes(y=SYNpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Synchrony",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

i26<-ggplot(data=portfolio,aes(y=SYNpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

i27<-ggplot(data=portfolio,aes(y=SYNpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

#5.4.4# synchrony~Chla

i28<-ggplot(data=portfolio,aes(y=SYNpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Synchrony",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

i29<-ggplot(data=portfolio,aes(y=SYNpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

i30<-ggplot(data=portfolio,aes(y=SYNpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

#5.4.5# synchrony~TAI

i31<-ggplot(data=portfolio,aes(y=SYNpart_tai,x=TAI))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Synchrony",x="TAI")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(M)")

i32<-ggplot(data=portfolio,aes(y=SYNpart_tai,x=TAI))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="TAI")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(N)")

i33<-ggplot(data=portfolio,aes(y=SYNpart_tai,x=TAI))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="TAI")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(O)")

pdf(file="SYN_cwm_int.pdf",width=9,height=12)
(i19+i20+i21)/(i22+i23+i24)/(i25+i26+i27)/(i28+i29+i30)/(i31+i32+i33)
dev.off()

#5.5# Univariate relationships for Evenness partial effects

#5.5.1# evenness~SBTrange

i34<-ggplot(data=portfolio,aes(y=EVENpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Evenness",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

i35<-ggplot(data=portfolio,aes(y=EVENpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

i36<-ggplot(data=portfolio,aes(y=EVENpart_sbt,x=SBTrange))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="SBTrange")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#5.5.2# evenness~Chla

i37<-ggplot(data=portfolio,aes(y=EVENpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Evenness",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

i38<-ggplot(data=portfolio,aes(y=EVENpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

i39<-ggplot(data=portfolio,aes(y=EVENpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#5.5.3# evenness~DepthCV

i40<-ggplot(data=portfolio,aes(y=EVENpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Evenness",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

i41<-ggplot(data=portfolio,aes(y=EVENpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

i42<-ggplot(data=portfolio,aes(y=EVENpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

#5.5.4# evenness~Trawl

i43<-ggplot(data=portfolio,aes(y=EVENpart_trawl,x=Trawl))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Evenness",x="Trawl")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

i44<-ggplot(data=portfolio,aes(y=EVENpart_trawl,x=Trawl))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Trawl")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

i45<-ggplot(data=portfolio,aes(y=EVENpart_trawl,x=Trawl))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Trawl")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

pdf(file="EVEN_cwm_int.pdf",width=12,height=10)
(i34+i35+i36)/(i37+i38+i39)/(i40+i41+i42)/(i43+i44+i45)
dev.off()

#5.6# Univariate relationships for Richness partial effects

#5.6.1# richness~Chla

i46<-ggplot(data=portfolio,aes(y=RICHpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Richness",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

i47<-ggplot(data=portfolio,aes(y=RICHpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

i48<-ggplot(data=portfolio,aes(y=RICHpart_chla,x=Chla))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="Chla")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

#5.6.2# richness~ChlaCV

i49<-ggplot(data=portfolio,aes(y=RICHpart_chlacv,x=ChlaCV))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Richness",x="ChlaCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

i50<-ggplot(data=portfolio,aes(y=RICHpart_chlacv,x=ChlaCV))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="ChlaCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

i51<-ggplot(data=portfolio,aes(y=RICHpart_chlacv,x=ChlaCV))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="ChlaCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

#5.6.3# richness~DepthCV

i52<-ggplot(data=portfolio,aes(y=RICHpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=Kgroup))+
  geom_smooth(method="lm",se=F,aes(col=Kgroup))+
  facet_wrap(~Kgroup)+
  labs(y="Richness",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

i53<-ggplot(data=portfolio,aes(y=RICHpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=L50group))+
  geom_smooth(method="lm",se=F,aes(col=L50group))+
  facet_wrap(~L50group)+
  labs(y="",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

i54<-ggplot(data=portfolio,aes(y=RICHpart_depth,x=DepthCV))+
  geom_point(size=3,alpha=0.4,aes(col=TLgroup))+
  geom_smooth(method="lm",se=F,aes(col=TLgroup))+
  facet_wrap(~factor(TLgroup,levels=c("Low TL","High TL")))+
  scale_color_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="",x="DepthCV")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

pdf(file="RICH_cwm_int.pdf",width=10,height=8)
(i46+i47+i48)/(i49+i50+i51)/(i52+i53+i54)
dev.off()


# ------------------------------------------------------ #
# ------------------------------------------------------ #
#6# Fit SEM to AvCV_PE using lavaan

m6.step<-lm(avcvpe~z.value+synchrony+evenness+richness+
              SBT+SBTrange+Chla+ChlaCV+Depth+DepthCV+Trawl+TAI,data=portfolio)
step(m6.step)

#6.1# Same model structure as for MV_PE (except for the avcvpe submodel)

peSEM.cv<-'avcvpe    ~ z.value + synchrony + richness + DepthCV
           z.value   ~ evenness + richness + SBT + Depth
           synchrony ~ evenness + richness + SBTrange + Chla + TAI
           evenness  ~ SBTrange + Chla + DepthCV + Trawl
           richness  ~ Chla + ChlaCV + DepthCV'

partial.peSEMCV<-sem(peSEM.cv,data=portfolio)
summary(partial.peSEMCV,rsquare=T,standardized=T)
lavaanPlot(model=partial.peSEMCV,coefs=T,stand=T,graph_options=list(layout="circo"),sig=0.05)

#6.2# Calculating direct and indirect effects for AvCV_PE using lavaan syntax

#6.2.1# Calculation within the model formulation

modelEffCVPE<-'avcvpe    ~ cz*z.value + cs*synchrony + cr*richness + cd*DepthCV
               z.value   ~ ze*evenness + zr*richness + zs*SBT + zd*Depth 
               synchrony ~ se*evenness + sr*richness + ss*SBTrange + sc*Chla + st*TAI 
               evenness  ~ es*SBTrange + ec*Chla + ed*DepthCV + et*Trawl
               richness  ~ ra*Chla + rc*ChlaCV + rd*DepthCV
           
               directZ          := cz
               directSyn        := cs
               directRich       := cr
               directDepthCV    := cd
               
               indirectRich     := (zr*cz) + (sr*cs)
               indirectEven     := (ze*cz) + (se*cs)
               
               indirectSBT      := zs*cz
               indirectSBTrange := (ss*cs) + (es*se*cs) + (es*ze*cz)
               
               indirectChla     := (sc*cs) + (ec*se*cs) + (ec*ze*cz) + (ra*sr*cs) + (ra*zr*cz) + (ra*cr)
               indirectChlaCV   := (rc*sr*cs) + (rc*zr*cz) + (rc*cr)
               
               indirectDepth    := zd*cz
               indirectDepthCV  := (ed*se*cs) + (ed*ze*cz) + (rd*sr*cs) + (rd*zr*cz) + (rd*cr)
               
               indirectTAI      := st*cs
               indirectTrawl    := (et*se*cs) + (et*ze*cz)
               
               totalZ           := directZ
               totalSyn         := directSyn
               totalRich        := directRich + indirectRich
               totalEven        := indirectEven
               totalSBT         := indirectSBT
               totalSBTrange    := indirectSBTrange
               totalChla        := indirectChla
               totalChlaCV      := indirectChlaCV
               totalDepth       := indirectDepth
               totalDepthCV     := directDepthCV + indirectDepthCV
               totalTAI         := indirectTAI
               totalTrawl       := indirectTrawl'

partial.defSEM.CVPEeffects<-sem(modelEffCVPE,data=portfolio)
summary(partial.defSEM.CVPEeffects,standardized=T)

#6.2.2# Plot direct and indirect effects

sem.CVPEeffects<-data.frame(summary(partial.defSEM.CVPEeffects,standardized=T)$pe[62:87,])

sem.CVPEeffects$varName<-c("z-value","Synchrony","Richness","DepthCV",
                           "Richness","Evenness","SBT","SBTrange","Chla",
                           "ChlaCV","Depth","DepthCV","TAI","Trawl","z-value",
                           "Synchrony","Richness","Evenness","SBT","SBTrange",
                           "Chla","ChlaCV","Depth","DepthCV","TAI","Trawl")

sem.CVPEeffects$effType<-c("Direct","Direct","Direct","Direct",
                           "Indirect","Indirect","Indirect","Indirect","Indirect",
                           "Indirect","Indirect","Indirect","Indirect","Indirect",
                           "Total","Total","Total","Total","Total","Total",
                           "Total","Total","Total","Total","Total","Total")

sem.CVPEeffects$varNum<-factor(c(1,2,3,10,3,4,5,6,7,8,9,10,11,12,
                                 1,2,3,4,5,6,7,8,9,10,11,12))

pdf(file="CVPE_SEM_dir_indir.pdf")

ggplot(data=filter(sem.CVPEeffects,effType != "Total"),
       aes(y=varNum,x=std.all))+
  geom_col(aes(fill=effType))+
  #geom_point(data=filter(sem.CVPEeffects,effType == "Total"),aes(y=varNum,x=std.all))+ # net effects
  #geom_point(data=filter(sem.CVPEeffects,effType == "Total" & pvalue < 0.05),
  #           aes(y=varNum,x=std.all),size=3,shape=8)+ # statistically significant total effects
  #geom_point(data=filter(sem.CVPEeffects,effType == "Direct" & pvalue < 0.05 & varName == "Richness"),
  #           aes(y=varNum,x=std.all),size=3,shape=8)+ # statistically significant direct richness effect
  #geom_point(data=filter(sem.CVPEeffects,effType == "Total" & pvalue > 0.05 & varName == "Richness"),
  #           aes(y=varNum,x=std.all),size=3,shape=4)+ # statistically non-significant total richness effect
  #geom_point(data=filter(sem.CVPEeffects,effType == "Indirect" & pvalue < 0.05),
  #           aes(y=varNum,x=std.all),size=3,shape=8)+ # statistically significant indirect effects
  #geom_point(data=filter(sem.CVPEeffects,effType == "Indirect" & pvalue > 0.05),
  #           aes(y=varNum,x=std.all),size=3,shape=4)+ # statistically significant indirect effects
  scale_y_discrete("",labels=c("1"="z-value","2"="Synchrony","3"="Richness",
                               "4"="Evenness","5"="SBT","6"="SBTrange",
                               "7"="Chla","8"="ChlaCV","9"="Depth",
                               "10"="DepthCV","11"="TAI","12"="Trawl"))+
  scale_x_continuous("Standardized effect on CVPE",limits=c(-0.75,0.5))+
  scale_fill_manual("",values=c("#E7B800","#00AFBB"))+
  geom_vline(xintercept=0,linetype="dashed")+
  theme_classic()+
  background_grid()+
  theme(legend.position=c(0,1),legend.justification=c(-0.4,1.2),
        axis.text=element_text(size=14),axis.title=element_text(size=16))

dev.off()

#6.3# SEM with picewiseSEM for AvCV_PE

psemCVPE1<-lm(avcvpe    ~ z.value + synchrony + richness + DepthCV,data=portfolio)
psemCVPE2<-lm(z.value   ~ evenness + richness + SBT + Depth,data=portfolio)
psemCVPE3<-lm(synchrony ~ evenness + richness + SBTrange + Chla + TAI,data=portfolio)
psemCVPE4<-lm(evenness  ~ SBTrange + Chla + DepthCV + Trawl,data=portfolio)
psemCVPE5<-lm(richness  ~ Chla + ChlaCV + DepthCV,data=portfolio)

psemCVPE.def<-psem(psemCVPE1,psemCVPE2,psemCVPE3,psemCVPE4,psemCVPE5)

coefs(psemCVPE.def,standardize="scale")
rsquared(psemCVPE.def)
dSep(psemCVPE.def)
fisherC(psemCVPE.def)
plot(psemCVPE.def)

vif(psemCVPE1);vif(psemCVPE2);vif(psemCVPE3);vif(psemCVPE4);vif(psemCVPE5)

#6.4# Residuals for each model included in the psem

portfolio$residCVPE1<-residuals(psemCVPE1) # CVPE model
portfolio$residCVPE2<-residuals(psemCVPE2) # z-value model
portfolio$residCVPE3<-residuals(psemCVPE3) # Synchrony model
portfolio$residCVPE4<-residuals(psemCVPE4) # Evenness model
portfolio$residCVPE5<-residuals(psemCVPE5) # Richness model

portfolio$fitCVPE1<-fitted(psemCVPE1)
portfolio$fitCVPE2<-fitted(psemCVPE2)
portfolio$fitCVPE3<-fitted(psemCVPE3)
portfolio$fitCVPE4<-fitted(psemCVPE4)
portfolio$fitCVPE5<-fitted(psemCVPE5)

reCVPE1<-ggplot(data=portfolio,aes(x=fitCVPE1,y=residCVPE1))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

reCVPE2<-ggplot(data=portfolio,aes(x=residCVPE1))+
  geom_histogram(binwidth=0.2)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

reCVPE3<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residCVPE1),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

reCVPE4<-ggplot(data=portfolio,aes(x=fitCVPE2,y=residCVPE2))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

reCVPE5<-ggplot(data=portfolio,aes(x=residCVPE2))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(E)")

reCVPE6<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residCVPE2),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(F)")

reCVPE7<-ggplot(data=portfolio,aes(x=fitCVPE3,y=residCVPE3))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(G)")

reCVPE8<-ggplot(data=portfolio,aes(x=residCVPE3))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(H)")

reCVPE9<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residCVPE3),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(I)")

reCVPE10<-ggplot(data=portfolio,aes(x=fitCVPE4,y=residCVPE4))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(J)")

reCVPE11<-ggplot(data=portfolio,aes(x=residCVPE4))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(K)")

reCVPE12<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residCVPE4),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(L)")

reCVPE13<-ggplot(data=portfolio,aes(x=fitCVPE5,y=residCVPE5))+
  geom_point(alpha=0.4,size=4)+
  geom_hline(yintercept=0)+
  labs(y="Residuals",x="Fitted")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(M)")

reCVPE14<-ggplot(data=portfolio,aes(x=residCVPE5))+
  geom_histogram(binwidth=0.4)+
  labs(y="Count",x="Residuals")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(N)")

reCVPE15<-ggplot(data=portfolio,aes(x=Longitude,y=Latitude))+
  geom_polygon(data=m.world,aes(x=long,y=lat,group=group),
               fill="gray80",col="white")+
  coord_cartesian(xlim=c(-95,35),ylim=c(10,70))+
  geom_point(aes(col=residCVPE5),size=2)+
  scale_colour_gradient2("Residuals",low=("blue"),
                         mid=("white"),high=("red"),midpoint=0)+
  labs(y="Latitude (ºN)",x="Longitude (ºW-E)")+
  theme_classic()+
  background_grid()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(O)")

pdf(file="resids_CVPE.pdf",width=12,height=15)
(reCVPE1+reCVPE2+reCVPE3)/(reCVPE4+reCVPE5+reCVPE6)/
  (reCVPE7+reCVPE8+reCVPE9)/(reCVPE10+reCVPE11+reCVPE12)/
  (reCVPE13+reCVPE14+reCVPE15)
dev.off()

#6.5# Partial plots for AvCV_PE SEM

coef(partial.peSEMCV)

#6.5.1# Partial z-value effect on AvCV_PE

portfolio$CV_PEpart_z<- portfolio$avcvpe -
  (coef(partial.peSEMCV)[2]*portfolio$synchrony) -
  (coef(partial.peSEMCV)[3]*portfolio$richness) -
  (coef(partial.peSEMCV)[4]*portfolio$DepthCV)

x.line.cvpe.z<-seq(from=min(portfolio$z.value),to=max(portfolio$z.value),by=0.1)
y.line.cvpe.z<-coef(partial.peSEMCV)[1]*x.line.cvpe.z
line.cvpe.z<-as.data.frame(cbind(x.line.cvpe.z,y.line.cvpe.z))

#6.5.2# Partial synchrony effect on AvCV_PE

portfolio$CV_PEpart_syn<- portfolio$avcvpe -
  (coef(partial.peSEMCV)[1]*portfolio$z.value) -
  (coef(partial.peSEMCV)[3]*portfolio$richness) -
  (coef(partial.peSEMCV)[4]*portfolio$DepthCV)

x.line.cvpe.syn<-seq(from=min(portfolio$synchrony),to=max(portfolio$synchrony),by=0.1)
y.line.cvpe.syn<-coef(partial.peSEMCV)[2]*x.line.cvpe.syn
line.cvpe.syn<-as.data.frame(cbind(x.line.cvpe.syn,y.line.cvpe.syn))

#6.5.3# Partial richness effect on AvCV_PE

portfolio$CV_PEpart_rich<- portfolio$avcvpe -
  (coef(partial.peSEMCV)[1]*portfolio$z.value) -
  (coef(partial.peSEMCV)[2]*portfolio$synchrony) -
  (coef(partial.peSEMCV)[4]*portfolio$DepthCV)

x.line.cvpe.rich<-seq(from=min(portfolio$richness),to=max(portfolio$richness),by=0.1)
y.line.cvpe.rich<-coef(partial.peSEMCV)[3]*x.line.cvpe.rich
line.cvpe.rich<-as.data.frame(cbind(x.line.cvpe.rich,y.line.cvpe.rich))

#6.5.4# Partial DepthCV effect on AvCV_PE

portfolio$CV_PEpart_depthcv<- portfolio$avcvpe -
  (coef(partial.peSEMCV)[1]*portfolio$z.value) -
  (coef(partial.peSEMCV)[2]*portfolio$synchrony) -
  (coef(partial.peSEMCV)[3]*portfolio$richness)

x.line.cvpe.depthcv<-seq(from=min(portfolio$DepthCV),to=max(portfolio$DepthCV),by=0.1)
y.line.cvpe.depthcv<-coef(partial.peSEMCV)[4]*x.line.cvpe.depthcv
line.cvpe.depthcv<-as.data.frame(cbind(x.line.cvpe.depthcv,y.line.cvpe.depthcv))

#6.5.5# Plots

cvpe1<-ggplot()+
  geom_point(data=portfolio,aes(y=CV_PEpart_z,x=z.value,col=Region))+
  geom_line(data=line.cvpe.z,aes(x=x.line.cvpe.z,y=y.line.cvpe.z),linewidth=1,col="gray50")+
  labs(y="CVPE",x="z-value")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(A)")

cvpe2<-ggplot()+
  geom_point(data=portfolio,aes(y=CV_PEpart_syn,x=synchrony,col=Region))+
  geom_line(data=line.cvpe.syn,aes(x=x.line.cvpe.syn,y=y.line.cvpe.syn),linewidth=1,col="gray50")+
  labs(y="",x="Synchrony")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(B)")

cvpe3<-ggplot()+
  geom_point(data=portfolio,aes(y=CV_PEpart_rich,x=richness,col=Region))+
  geom_line(data=line.cvpe.rich,aes(x=x.line.cvpe.rich,y=y.line.cvpe.rich),linewidth=1,col="gray50")+
  labs(y="CVPE",x="Richness")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(C)")

cvpe4<-ggplot()+
  geom_point(data=portfolio,aes(y=CV_PEpart_depthcv,x=DepthCV,col=Region))+
  geom_line(data=line.cvpe.depthcv,aes(x=x.line.cvpe.depthcv,y=y.line.cvpe.depthcv),linewidth=1,col="gray50")+
  labs(y="",x="DepthCV")+
  scale_color_brewer(palette="Set1")+
  theme_classic()+
  background_grid()+
  theme(legend.position="none",
        axis.text=element_text(size=14),axis.title=element_text(size=16))+
  ggtitle("(D)")

pdf(file="CVPE_SEM_partial_plots.pdf",width=6,height=6)
cvpe1+cvpe2+cvpe3+cvpe4+plot_layout(ncol=2)
dev.off()

