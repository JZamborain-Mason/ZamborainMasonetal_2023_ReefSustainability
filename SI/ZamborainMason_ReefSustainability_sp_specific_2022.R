# Title: Reef sustainability_SI (species-specific analyses)
# Author: Jessica Zamborain Mason
# Description: This script Calculates total species richness based on sp-abundance distribution from serf data, then it gets individual species intrinsic growth rates, 
# and shows how it is not straightforward to know the relationship between community growth rates and individual species growth rates (community intrinsic growth rates are not necessary representative of individual species  (e.g., trophic dynamics and species interactions))
# R version 3.5.3 (2019-03-11)


#clear R
rm(list=ls())

#set working directory (where the data and code of this repository is stored)
setwd("")


#load required libraries
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library (poilog)#Vidar Grøtan and Steinar Engen (2008). poilog: Poisson lognormal and bivariate Poisson lognormal distribution. R package version 0.4.
library (rfishbase)#C. Boettiger, D. T. Lang and P. C. Wainwright. "rfishbase: exploring, manipulating and visualizing FishBase data from R". In: _Journal of Fish Biology_ 81.6 (Nov. 2012), pp. 2030-2039. DOI: 10.1111/j.1095-8649.2012.03464.x
library(FishLife)#James Thorson (2019). FishLife: Predict Life History Parameters For Any Fish. R package version 2.0.0.
library(tidyverse)#Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686


#upload data 
sp_specificdata<-read.csv("individualspscale_data_submitted.csv", head=T)
reefscale_data<-read.csv("reefscale_data_submitted.csv", head=T)

###########################################################################3
#Species richness

#abundance distributions per uniquesite
unique(sp_specificdata$Genus)
abundancedata<-sp_specificdata[!sp_specificdata$Genus=="unknown",]
abundancedata<-droplevels(abundancedata)
length(unique(abundancedata$UniqueSite)) #those sites that provider gave to sp level

#loop to go over each column and repeat each species the number of times as "number"
z<-list()#create a list
for (i in 1:length(abundancedata$FullSpecies)){ #from 1 to the length of the dataframe
  zz<-rep(abundancedata$FullSpecies[i],length.out=abundancedata$Number[i]) #repeat the species the number of times as the "number"
  zz<- paste(zz,abundancedata$UniqueSite[i], sep="_") #paste the species to its unique site
  z[[length(z)+1]]<-zz #store it after each species
}
z<- unlist(z, use.names=FALSE)#make it into  a vector
abundance<-colsplit(z,"\\_", names=c("FullSpecies","UniqueSite"))#make it into a data.frame splitting the columns
uniquesites<-unique(abundance$UniqueSite) #create a vector with the unique unique sites

#calculate the species richness for all uniquesites
A<-data.frame() #empty data frame to store the frequencies per species per unique site
p<-vector()
pars<-matrix(NA, ncol=2,nrow=length(uniquesites))
mle<-vector()
gof<-vector()
nbsp_us<-vector()
total_sp<-vector()
for (i in 1: length(uniquesites)){
  site_freq<-table(abundance$FullSpecies[abundance$UniqueSite==uniquesites[i]])  #frequency table for the first uniquesite and so forth
  B<-as.data.frame(site_freq)#make the frequency table into a data-frame
  B$UniqueSite<-rep(uniquesites[i],length.out=length(B$Freq))#put the uniquesite
  tri<-poilogMLE(sort(B$Freq), startVals = c(mu=2, sig=3),
                nboot = 10, zTrunc = TRUE,
                method = "BFGS", control = list(maxit=100000))
  p[i]<-tri$p
  pars[i,]<-tri$par
  mle[i]<-tri$logLval
  gof[i]<-tri$gof
  nbsp_us[i]<-length(B$Var1)
  total_sp[i]<-nbsp_us[i]/p[i]
  A<-rbind(A,B) #combine the for-loops
}

us_spdata<-as.data.frame(cbind(uniquesites,nbsp_us,total_sp))
colnames(us_spdata)<-c("UniqueSite","nbsp_us","total_sp")
colnames(pars)<-c("mu","sigma")
us_spdata<-cbind(us_spdata,pars,gof)
us_spdata$goodfit<-ifelse(us_spdata$gof>0.05 & us_spdata$gof<0.95,1,0 ) #following the poilog package include those with good fit
us_spdata$total_sp<-ifelse(us_spdata$goodfit==1,us_spdata$total_sp,NA)
us_spdata$total_sp2<-us_spdata$total_sp

#check it is equal to the one used
serfsites<-merge(reefscale_data,us_spdata[,c("UniqueSite","total_sp2")], by="UniqueSite",all.x=T)
ggplot(serfsites, aes(x=log(total_sp), y=log(total_sp2)))+geom_point()+geom_abline(intercept=0, slope=1)
serfsites$Larger[!is.na(serfsites$total_sp2)&serfsites$total_sp2>1000]

#to download sp richness
speciesrichness<-serfsites[,c("UniqueSite", "total_sp2")]
colnames(speciesrichness)<-c("UniqueSite", "total_sp")
#write.csv(speciesrichness, "sitespeciesrichness.csv", row.names=F)

###############################################################################
#individual species intrinsic growth rates vs. community growth rates

#benchmark data
benchmarkdata<-reefscale_data[reefscale_data$section=="benchmark" &(reefscale_data$reserves==1|(reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0)),]
sp_specificdata$dataset<-NULL
benchmarkdata_sp<-merge(benchmarkdata, sp_specificdata, by="UniqueSite", all.x=T)
benchmarkdata_sp$dataset[is.na(benchmarkdata_sp$FullSpecies)]
benchmarkdata_sp<-benchmarkdata_sp[!is.na(benchmarkdata_sp$FullSpecies),]

z<-list()
for (i in 1:length(benchmarkdata_sp$FullSpecies)){ 
  zz<-rep(benchmarkdata_sp$FullSpecies[i],length.out=benchmarkdata_sp$Number[i]) 
  zz<- paste(zz,benchmarkdata_sp$UniqueSite[i], sep="_") 
  z[[length(z)+1]]<-zz 
}
z<- unlist(z, use.names=FALSE)
benchmark_abundance<-colsplit(z,"\\_", names=c("FullSpecies","UniqueSite"))

#get mean intrinsic growth rate for each individual in the dataset to the lowest taxonomic level possible
uniquefamilies<-as.data.frame(unique(benchmarkdata_sp$Family[!is.na(benchmarkdata_sp$Family)]))
colnames(uniquefamilies)<-"Family"
uniquegenus<-as.data.frame(unique(benchmarkdata_sp$Genus[!(is.na(benchmarkdata_sp$Genus)|benchmarkdata_sp$Genus=="sp." |benchmarkdata_sp$Genus=="unknown")]))
colnames(uniquegenus)<-"Genus"
uniquesp<-as.data.frame(unique(benchmarkdata_sp$FullSpecies[!(is.na(benchmarkdata_sp$FullSpecies)| benchmarkdata_sp$Species=="sp." |benchmarkdata_sp$Species=="unknown")]))
colnames(uniquesp)<-"FullSpecies"
uniquesp_info<-ddply(benchmarkdata_sp,.(FullSpecies),summarize,Genus=Genus[1],Species=Species[1], Family=Family[1])
uniquesp<-merge(uniquesp,uniquesp_info,by="FullSpecies",all.x=T)
uniquefamilies$r<-rep(NA, length(uniquefamilies$Family))
uniquegenus$r<-rep(NA, length(uniquegenus$Genus))
uniquesp$r<-rep(NA, length(uniquesp$FullSpecies))

for (i in 1:length(uniquefamilies$Family)){
  trial<-Plot_taxa( Search_species(Family=as.character(uniquefamilies$Family[i]))$match_taxonomy[1])
  uniquefamilies$r[i]<-exp(trial[[1]]$Mean_pred[17])
}
summary(uniquefamilies$r)

for (i in 1:length(uniquegenus$Genus)){
  trial<-Plot_taxa( Search_species(Genus=as.character(uniquegenus$Genus[i]))$match_taxonomy[1])
  uniquegenus$r[i]<-exp(trial[[1]]$Mean_pred[17])
}
summary(uniquegenus$r)

#some species dont match so we correct them
uniquesp$Genus<-as.factor(ifelse(uniquesp$Genus=="Carangoides" &uniquesp$Species=="ruber", "Caranx", as.character(uniquesp$Genus)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="flavicauda", "fisheri", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="loricula", "loriculus", as.character(uniquesp$Species)))
uniquesp$Genus<-as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="multifasciata", "Paracentropyge", as.character(uniquesp$Genus)))
uniquesp$Genus[uniquesp$FullSpecies=="Diplodus cervinus omanensis"]
uniquesp$Species[uniquesp$FullSpecies=="Diplodus cervinus omanensis"]
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Diplodus" &uniquesp$Species=="cervinus omanensis", "cervinus", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Kyphosus" &uniquesp$Species=="sectator", "sectatrix", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Macropharyngodon" &uniquesp$Species=="bipartitus marisrubri", "bipartitus", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Plectorhinchus" &uniquesp$Species=="orientalis", "vittatus", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Plectorhinchus" &uniquesp$Species=="unicolor", "schotaf", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="bimaculatus", "bimaculata", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="frenatus", "frenata", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="taeniatus", "taeniata", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Thalassoma" &uniquesp$Species=="duperrey/lutescens", "duperrey", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Upeneus" &uniquesp$Species=="arge", "taeniopterus", as.character(uniquesp$Species)))
uniquesp$Species<-as.factor(ifelse(uniquesp$Genus=="Zebrasoma" &uniquesp$Species=="veliferum", "velifer", as.character(uniquesp$Species)))

for (i in 1:length(uniquesp$FullSpecies)){
  trial<-Plot_taxa( Search_species(Genus=as.character(uniquesp$Genus[i]),Species=as.character(uniquesp$Species[i]))$match_taxonomy[1])
  uniquesp$r[i]<-exp(trial[[1]]$Mean_pred[17])
}
summary(uniquesp$r)

#MERGE INDIVIDUAL sp, genus and family intrinsic growth rates TO THE DATASET
uniquesp$sp_r<-uniquesp$r
uniquegenus$genus_r<-uniquegenus$r
uniquefamilies$family_r<-uniquefamilies$r
benchmark_abundance<-merge(benchmark_abundance,uniquesp_info[,c("FullSpecies","Genus","Family")], by="FullSpecies",all.x=T)
benchmark_abundance<-merge(benchmark_abundance, uniquesp[,c("sp_r","FullSpecies")], by="FullSpecies",all.x=T)
summary(benchmark_abundance$sp_r)
benchmark_abundance<-merge(benchmark_abundance, uniquegenus[,c("genus_r","Genus")], by="Genus",all.x=T)
summary(benchmark_abundance$genus_r)
benchmark_abundance<-merge(benchmark_abundance, uniquefamilies[,c("family_r","Family")], by="Family",all.x=T)
summary(benchmark_abundance$family_r)
benchmark_abundance$good_r<-ifelse(is.na(benchmark_abundance$sp_r)&is.na(benchmark_abundance$genus_r), benchmark_abundance$family_r, ifelse(is.na(benchmark_abundance$sp_r)& !is.na(benchmark_abundance$genus_r), benchmark_abundance$genus_r, benchmark_abundance$sp_r))
summary(benchmark_abundance$good_r)

#plot intrinsic growth rate
fig1<-benchmark_abundance%>%
  ggplot() +
  geom_histogram(aes(good_r),alpha=0.5) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(x = "")+theme_classic()+
  geom_vline(xintercept =0.08,lwd=1.3,col="darkred")+geom_vline(xintercept = c(0.05, 0.16),lty=2,color="darkred")+
  geom_text(aes(x=0.5,y=3500),label="Community
growth rate",col="darkred",family="Helvetica")+ggtitle("This study:reef fish")
benchmark_abundance$r=benchmark_abundance$good_r


windows()
fig1
#######################################################################################################
#do the same for the Worm et al. paper (assuming Graham-Schaefer surplus dynamics; i.e., assuming community growth rate is equal sustainable exploitation rate at MMSY *2  )
#using reported species and ummsy reported in SI

wormcalifornia<-read.csv("wormetal2009_californiacurrentsp.csv", header=T)
summary(wormcalifornia)

wormcalifornia$FullSpecies=paste(wormcalifornia$genus, wormcalifornia$species, sep="_")
length(unique(wormcalifornia$FullSpecies))
reefquery2<-paste0(wormcalifornia$FullSpecies, collapse = '|')
int_rate2<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery2) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(int_rate2$rowname)

mean(int_rate2$r)
median(int_rate2$r)

#FMSY= r*2 (FAO) under Graham-Schaefer
a<-int_rate2 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.03*2,.07*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "California Current",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

#do for all the different ecosystems
wormdata<-read.csv("wormetal2009_spdata.csv", header=T)
summary(wormdata)

wormdata$FullSpecies<-paste(wormdata$genus, wormdata$species, sep="_")
length(unique(wormdata$FullSpecies))
wormdata$location<-paste(wormdata$name.1, wormdata$name.2, sep=" ")

#do it for each ecosystem
length(unique(wormdata$location))
unique(wormdata$location)
reefquery3<-paste0(wormdata$FullSpecies[wormdata$location=="Iceland Shelf"], collapse = '|')

int_rate3<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery3) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Iceland Shelf"]))
length(int_rate3$rowname)

b<-int_rate3 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.23*2,.34*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Iceland Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery4<-paste0(wormdata$FullSpecies[wormdata$location=="North Sea"], collapse = '|')
int_rate4<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery4) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="North Sea"]))
length(int_rate4$rowname)
c<-int_rate4 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.08*2,.16*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "North Sea",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery5<-paste0(wormdata$FullSpecies[wormdata$location=="Celtic-Biscay Shelf"], collapse = '|')
int_rate5<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery5) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Celtic-Biscay Shelf"]))
length(int_rate5$rowname)
d<-int_rate5 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.08*2,.17*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Celtic-Biscay Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery6<-paste0(wormdata$FullSpecies[wormdata$location=="Southern Australian"], collapse = '|')
int_rate6<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery6) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Southern Australian"]))
length(int_rate6$rowname)
e<-int_rate6 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.12*2,.18*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Southern Australian Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery7<-paste0(wormdata$FullSpecies[wormdata$location=="Northeast U.S."], collapse = '|')
int_rate7<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery7) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Northeast U.S."]))
length(int_rate7$rowname)
f<-int_rate7 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.20*2,.32*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Northeast U.S.Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery8<-paste0(wormdata$FullSpecies[wormdata$location=="Newfoundland-Labrador Shelf"], collapse = '|')
int_rate8<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery8) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Newfoundland-Labrador Shelf"]))
length(int_rate8$rowname)
g<-int_rate8 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.20*2,.26*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Newfoundland-Labrador Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery9<-paste0(wormdata$FullSpecies[wormdata$location=="Baltic Sea"], collapse = '|')
int_rate9<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery9) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Baltic Sea"]))
length(int_rate9$rowname)
h<-int_rate9 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.08*2,.12*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Baltic Sea",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery10<-paste0(wormdata$FullSpecies[wormdata$location=="Eastern Bering"], collapse = '|')
int_rate10<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery10) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Eastern Bering"]))
length(int_rate10$rowname)
i<-int_rate10 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.14*2,.21*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Eastern Bering",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()

reefquery11<-paste0(wormdata$FullSpecies[wormdata$location=="New Zealand"], collapse = '|')
int_rate11<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery11) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="New Zealand"]))
length(int_rate11$rowname)
j<-int_rate11 %>%
  ggplot() +
  geom_histogram(aes(r)) +
  geom_vline(xintercept = c(0.08*2,.11*2), color = "darkred") +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "New Zealand Shelf",caption = "")+labs(x = "sp intrinsic growth rates")+theme_classic()
windows()
ggarrange(fig1,a,b,c,d,e,f,g,h,i,j, nrow=2,ncol=6)


#worm using ln_r instead of r
a=int_rate2 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.03*2,.07*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "California Current",caption = "")+theme_classic()+labs(x = "")

#do ot for all the different ecosystems
wormdata=read.csv("wormetal2009_spdata.csv", header=T)
summary(wormdata)

wormdata$FullSpecies=paste(wormdata$genus, wormdata$species, sep="_")
length(unique(wormdata$FullSpecies))
wormdata$location=paste(wormdata$name.1, wormdata$name.2, sep=" ")

#do it for each ecosystem
length(unique(wormdata$location))
unique(wormdata$location)
reefquery3<-paste0(wormdata$FullSpecies[wormdata$location=="Iceland Shelf"], collapse = '|')

int_rate3<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery3) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Iceland Shelf"]))
length(int_rate3$rowname)

b=int_rate3 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.23*2,.34*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Iceland Shelf",caption = "")+theme_classic()+labs(x = "")

reefquery4<-paste0(wormdata$FullSpecies[wormdata$location=="North Sea"], collapse = '|')
int_rate4<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery4) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="North Sea"]))
length(int_rate4$rowname)
c=int_rate4 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.08*2,.16*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "North Sea",caption = "")+theme_classic()+labs(x = "")

reefquery5<-paste0(wormdata$FullSpecies[wormdata$location=="Celtic-Biscay Shelf"], collapse = '|')
int_rate5<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery5) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Celtic-Biscay Shelf"]))
length(int_rate5$rowname)
d=int_rate5 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.08*2,.17*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Celtic-Biscay Shelf",caption = "")+theme_classic()+labs(x = "")

reefquery6<-paste0(wormdata$FullSpecies[wormdata$location=="Southern Australian"], collapse = '|')
int_rate6<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery6) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Southern Australian"]))
length(int_rate6$rowname)
e=int_rate6 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.12*2,.18*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Southern Australian Shelf",caption = "")+theme_classic()+labs(x = "")

reefquery7<-paste0(wormdata$FullSpecies[wormdata$location=="Northeast U.S."], collapse = '|')
int_rate7<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery7) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Northeast U.S."]))
length(int_rate7$rowname)
f=int_rate7 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.20*2,.32*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Northeast U.S.Shelf",caption = "")+theme_classic()+labs(x = "")

reefquery8<-paste0(wormdata$FullSpecies[wormdata$location=="Newfoundland-Labrador Shelf"], collapse = '|')
int_rate8<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery8) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Newfoundland-Labrador Shelf"]))
length(int_rate8$rowname)
g=int_rate8 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.20*2,.26*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Newfoundland-Labrador Shelf",caption = "")+theme_classic()+labs(x = "")

reefquery9<-paste0(wormdata$FullSpecies[wormdata$location=="Baltic Sea"], collapse = '|')
int_rate9<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery9) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Baltic Sea"]))
length(int_rate9$rowname)
h=int_rate9 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.08*2,.12*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Baltic Sea",caption = "")+theme_classic()+labs(x = "")

reefquery10<-paste0(wormdata$FullSpecies[wormdata$location=="Eastern Bering"], collapse = '|')
int_rate10<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery10) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="Eastern Bering"]))
length(int_rate10$rowname)
i=int_rate10 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.14*2,.21*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "Eastern Bering",caption = "")+theme_classic()+labs(x = "")

reefquery11<-paste0(wormdata$FullSpecies[wormdata$location=="New Zealand"], collapse = '|')
int_rate11<- FishLife::FishBase_and_RAM$beta_gv %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(reef = str_detect(rowname, reefquery11) & !str_detect(rowname, "predictive")) %>%
  filter(reef)
length(unique(wormdata$FullSpecies[wormdata$location=="New Zealand"]))
length(int_rate11$rowname)
j=int_rate11 %>%
  ggplot() +
  geom_histogram(aes(exp(ln_r))) +
  geom_vline(xintercept = c(0.08*2,.11*2), color = "darkred",lty=2) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(title = "New Zealand Shelf",caption = "")+theme_classic()+labs(x = "")
windows()
worm_fig<-ggarrange(a,b,c,d,e,f,g,h,i,j, nrow=2,ncol=5)
spvscom<-ggarrange(fig1,worm_fig,labels=c("a","b"),nrow=1,ncol=2,widths=c(1,4))
annotate_figure(spvscom, bottom="Sp. intrinsic growth rates")

#######################################################################################################################################
#save.image(file='ZamborainMasonetal_ReefSustainability_spspecific.RData')
