#Code for Zamborain-Mason et al. 2020: Sustainability of the world's coral reef fisheries


#This script implements the individual-fish specific analyses: 
#1. Calculates total species richness ased on sp-abundance distribution
#2. Estimatea species-specific intrinsic growth rates for our data
#3. Runs simulation to prove individual sp intrinsic growth rates do not necceseraly translate into community ones

#clear R
rm(list=ls())

#set working directory (where the data and code of this repository is stored)
setwd("c:/Users/jzamb/Documents/AUSTRALIA/PhD/Chapter 1/ALL ANALYSES/FINAL CODE AND DATA/FINAL AFTER REVISIONS (GITHUB)")


#load required libraries
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library (poilog)#Vidar Grøtan and Steinar Engen (2008). poilog: Poisson lognormal and bivariate Poisson lognormal distribution. R package version 0.4.
library(FishLife)# James Thorson (2019). FishLife: Predict Life History Parameters For Any Fish. R package version 2.0.0. http://github.com/James-Thorson/FishLife
library (rfishbase)#C. Boettiger, D. T. Lang and P. C. Wainwright. "rfishbase: exploring, manipulating and visualizing FishBase data from R". In: _Journal of Fish Biology_ 81.6 (Nov. 2012), pp. 2030-2039. DOI: 10.1111/j.1095-8649.2012.03464.x

#upload data 
sp_specificdata=read.csv("individualspscale_data_submitted.csv", head=T)
reefscale_data=read.csv("reefscale_data_submitted.csv", head=T)

###########################################################################3
#Species richness

#abundance distributions per uniquesite
unique(sp_specificdata$Genus)
abundancedata=sp_specificdata[!sp_specificdata$Genus=="unknown",]
abundancedata=droplevels(abundancedata)
length(unique(abundancedata$UniqueSite)) #those sites that provider gave to sp level

#loop to go over each colums and repeat each species the number of times as "number"
z=list()#create a list
for (i in 1:length(abundancedata$FullSpecies)){ #from 1 to the length of the dataframe
  zz=rep(abundancedata$FullSpecies[i],length.out=abundancedata$Number[i]) #repeat the species the number of times as the "number"
  zz= paste(zz,abundancedata$UniqueSite[i], sep="_") #paste the species to its unique site
  z[[length(z)+1]]=zz #store it after each species
}
z= unlist(z, use.names=FALSE)#make it into  a vector

abundance=colsplit(z,"\\_", names=c("FullSpecies","UniqueSite"))#make it into a data.frame slitting the columns
length(unique(abundance$UniqueSite))
uniquesites=unique(abundance$UniqueSite) #create a vector with the unique unique sites

#calculate the species richness for all uniquesites
A=data.frame() #empty data frame to store the frequencies per species per unique site
p=vector()
pars=matrix(NA, ncol=2,nrow=length(uniquesites))
mle=vector()
gof=vector()
nbsp_us=vector()
total_sp=vector()
for (i in 1: length(uniquesites)){
  site_freq=table(abundance$FullSpecies[abundance$UniqueSite==uniquesites[i]])  #frequency table for the first unique site and so forth
  B=as.data.frame(site_freq)#make the frequency table into a data-frame
  B$UniqueSite=rep(uniquesites[i],length.out=length(B$Freq))#put the unique site
  tri=poilogMLE(sort(B$Freq), startVals = c(mu=2, sig=3),
                nboot = 10, zTrunc = TRUE,
                method = "BFGS", control = list(maxit=100000))
  p[i]=tri$p
  pars[i,]=tri$par
  mle[i]=tri$logLval
  gof[i]=tri$gof
  nbsp_us[i]=length(B$Var1)
  total_sp[i]=nbsp_us[i]/p[i]
  A=rbind(A,B) #combine the for loops
}

us_spdata=as.data.frame(cbind(uniquesites,nbsp_us,total_sp))
colnames(us_spdata)=c("UniqueSite","nbsp_us","total_sp")
colnames(pars)=c("mu","sigma")
summary(us_spdata$gof)
us_spdata=cbind(us_spdata,pars,gof)
us_spdata$goodfit=ifelse(us_spdata$gof>0.05 & us_spdata$gof<0.95,1,0 ) #following the poilog package include those with good fit
us_spdata$total_sp=ifelse(us_spdata$goodfit==1,us_spdata$total_sp,NA)
length(us_spdata$total_sp[!is.na(us_spdata$total_sp)])
us_spdata$total_sp2=us_spdata$total_sp
summary(us_spdata$total_sp2)
length(us_spdata$total_sp2[us_spdata$total_sp2>1000])

#check it is equal to the one used
serfsites=merge(reefscale_data,us_spdata[,c("UniqueSite","total_sp2")], by="UniqueSite",all.x=T)
ggplot(serfsites, aes(x=log(total_sp), y=log(total_sp2)))+geom_point()+geom_abline(intercept=0, slope=1)
serfsites$Larger[!is.na(serfsites$total_sp2)&serfsites$total_sp2>1000]

#to download sp richness
speciesrichness=serfsites[,c("UniqueSite", "total_sp2")]
colnames(speciesrichness)=c("UniqueSite", "total_sp")
#write.csv(speciesrichness, "sitespeciesrichness.csv", row.names=F)

###########################################################################3
#Intrinsic growth rates

#benchmark data
benchmarkdata=reefscale_data[reefscale_data$section=="benchmark" &(reefscale_data$reserves==1|(reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0)),]

colnames(sp_specificdata)
colnames(benchmarkdata)
sp_specificdata$dataset=NULL
benchmarkdata_sp=merge(benchmarkdata, sp_specificdata, by="UniqueSite", all.x=T)
benchmarkdata_sp$dataset[is.na(benchmarkdata_sp$FullSpecies)]
benchmarkdata_sp=benchmarkdata_sp[!is.na(benchmarkdata_sp$FullSpecies),]

z=list()#create a list
for (i in 1:length(benchmarkdata_sp$FullSpecies)){ #from 1 to the length of the dataframe
  zz=rep(benchmarkdata_sp$FullSpecies[i],length.out=benchmarkdata_sp$Number[i]) #repeat the species the number of times as the "number"
  zz= paste(zz,benchmarkdata_sp$UniqueSite[i], sep="_") #paste the species to its unique site
  z[[length(z)+1]]=zz #store it after each species
}
z= unlist(z, use.names=FALSE)#make it into  a vector

benchmark_abundance=colsplit(z,"\\_", names=c("FullSpecies","UniqueSite"))#make it into a data.frame slitting the columns

#get mean intrinsic growth rate for each individual in the dataset to the lowest taxonomic level possible
uniquefamilies=as.data.frame(unique(benchmarkdata_sp$Family[!is.na(benchmarkdata_sp$Family)]))
colnames(uniquefamilies)="Family"
uniquegenus=as.data.frame(unique(benchmarkdata_sp$Genus[!(is.na(benchmarkdata_sp$Genus)|benchmarkdata_sp$Genus=="sp." |benchmarkdata_sp$Genus=="unknown")]))
colnames(uniquegenus)="Genus"
uniquesp=as.data.frame(unique(benchmarkdata_sp$FullSpecies[!(is.na(benchmarkdata_sp$FullSpecies)| benchmarkdata_sp$Species=="sp." |benchmarkdata_sp$Species=="unknown")]))
colnames(uniquesp)="FullSpecies"
uniquesp_info=ddply(benchmarkdata_sp,.(FullSpecies),summarize,Genus=Genus[1],Species=Species[1], Family=Family[1])
uniquesp=merge(uniquesp,uniquesp_info,by="FullSpecies",all.x=T)
uniquefamilies$r=rep(NA, length(uniquefamilies$Family))
uniquegenus$r=rep(NA, length(uniquegenus$Genus))
uniquesp$r=rep(NA, length(uniquesp$FullSpecies))

for (i in 1:length(uniquefamilies$Family)){
  trial=Plot_taxa( Search_species(Family=as.character(uniquefamilies$Family[i]))$match_taxonomy[1])
  uniquefamilies$r[i]=exp(trial[[1]]$Mean_pred[17])
}
summary(uniquefamilies$r)

for (i in 1:length(uniquegenus$Genus)){
  trial=Plot_taxa( Search_species(Genus=as.character(uniquegenus$Genus[i]))$match_taxonomy[1])
  uniquegenus$r[i]=exp(trial[[1]]$Mean_pred[17])
}
summary(uniquegenus$r)

#some species dont match so we correct them
uniquesp$Genus=as.factor(ifelse(uniquesp$Genus=="Carangoides" &uniquesp$Species=="ruber", "Caranx", as.character(uniquesp$Genus)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="flavicauda", "fisheri", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="loricula", "loriculus", as.character(uniquesp$Species)))
uniquesp$Genus=as.factor(ifelse(uniquesp$Genus=="Centropyge" &uniquesp$Species=="multifasciata", "Paracentropyge", as.character(uniquesp$Genus)))
uniquesp$Genus[uniquesp$FullSpecies=="Diplodus cervinus omanensis"]
uniquesp$Species[uniquesp$FullSpecies=="Diplodus cervinus omanensis"]
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Diplodus" &uniquesp$Species=="cervinus omanensis", "cervinus", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Kyphosus" &uniquesp$Species=="sectator", "sectatrix", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Macropharyngodon" &uniquesp$Species=="bipartitus marisrubri", "bipartitus", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Plectorhinchus" &uniquesp$Species=="orientalis", "vittatus", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Plectorhinchus" &uniquesp$Species=="unicolor", "schotaf", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="bimaculatus", "bimaculata", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="frenatus", "frenata", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Scolopsis" &uniquesp$Species=="taeniatus", "taeniata", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Thalassoma" &uniquesp$Species=="duperrey/lutescens", "duperrey", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Upeneus" &uniquesp$Species=="arge", "taeniopterus", as.character(uniquesp$Species)))
uniquesp$Species=as.factor(ifelse(uniquesp$Genus=="Zebrasoma" &uniquesp$Species=="veliferum", "velifer", as.character(uniquesp$Species)))

for (i in 1:length(uniquesp$FullSpecies)){
  trial=Plot_taxa( Search_species(Genus=as.character(uniquesp$Genus[i]),Species=as.character(uniquesp$Species[i]))$match_taxonomy[1])
  uniquesp$r[i]=exp(trial[[1]]$Mean_pred[17])
}
summary(uniquesp$r)

#MERGE INDIVIDUAL sp, genus and family intrinsic growth rates TO THE DATASET
uniquesp$sp_r=uniquesp$r
uniquegenus$genus_r=uniquegenus$r
uniquefamilies$family_r=uniquefamilies$r
benchmark_abundance=merge(benchmark_abundance,uniquesp_info[,c("FullSpecies","Genus","Family")], by="FullSpecies",all.x=T)
benchmark_abundance=merge(benchmark_abundance, uniquesp[,c("sp_r","FullSpecies")], by="FullSpecies",all.x=T)
summary(benchmark_abundance$sp_r)
benchmark_abundance=merge(benchmark_abundance, uniquegenus[,c("genus_r","Genus")], by="Genus",all.x=T)
summary(benchmark_abundance$genus_r)
benchmark_abundance=merge(benchmark_abundance, uniquefamilies[,c("family_r","Family")], by="Family",all.x=T)
summary(benchmark_abundance$family_r)
benchmark_abundance$good_r=ifelse(is.na(benchmark_abundance$sp_r)&is.na(benchmark_abundance$genus_r), benchmark_abundance$family_r, ifelse(is.na(benchmark_abundance$sp_r)& !is.na(benchmark_abundance$genus_r), benchmark_abundance$genus_r, benchmark_abundance$sp_r))
summary(benchmark_abundance$good_r)
head(benchmark_abundance)

#plot intrinsic growth rate

fig1=benchmark_abundance%>%
  ggplot() +
  geom_histogram(aes(good_r)) +
  scale_y_continuous(expand = expand_scale()) +
  scale_x_continuous(expand = expand_scale()) +
  labs(x = "sp intrinsic growth rates")+theme_classic()
benchmark_abundance$r=benchmark_abundance$good_r

#relationship with biomass
bysite=ddply(benchmark_abundance,.(UniqueSite),summarize,median_r=median(good_r, na.rm=T),mean_r=mean(good_r, na.rm=T),countindividuals=length(good_r))
benchmarkdata=merge(benchmarkdata,bysite, by="UniqueSite", all.x=T)

fig2=ggplot(benchmarkdata, aes(y=mean_r,  x=log(FamBiomass_tkm2)) )+
  geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+xlab("log (Biomass (t/km2)")+ylab ("Mean sp. intrinsic growth rate")


#we are going to simulate the replacement of fast-growing sp to low-growing ones as biomass increases
#and tehn we are going to fit a logistic to their combined biomass

#assumptions
#each sp grows in biomass following a logistic and they grow independently
#sp with larger r's have lower k's and start at larger biomass values han those with lower k's
#starts with more biomass of higher r species

multisp.logistic<- function(init.pop,r, k,time,sp){
  # Create matrix with species (rows) by time(columns) 
  comunity=matrix(NA,nrow = sp, ncol = time)
  comunity[,1]=init.pop
  for (i in 2:time) {
    for (j in 1:sp) {
      comunity[j,i]=comunity[j,i-1]+r[j]*comunity[j,i-1]*(1-(comunity[j,i-1]/k[j]))
    }
  }
  # }
  return(comunity)
}

#set paramter and state variables
Number.Species=500
timesteps=50
Growth.rates2= rnorm(Number.Species, mean=mean(benchmark_abundance$r, na.rm=T),sd=sd(benchmark_abundance$r, na.rm=T))
Growth.rates2=ifelse(Growth.rates2<min(benchmark_abundance$r,na.rm=T),min(benchmark_abundance$r,na.rm=T),Growth.rates2)
carrying.capacities2=3/Growth.rates2
init.biomass2=Growth.rates2*1.5
init.biomass2=ifelse(Growth.rates2<0.4,0.1,init.biomass2)

#simulate community
simulatedcommunity=multisp.logistic(init.pop = init.biomass2, r = Growth.rates2, k = carrying.capacities2, time = timesteps, sp = Number.Species)
simulatedcommunity_df=as.data.frame(simulatedcommunity)
simulatedcommunity_df$r=Growth.rates2
simulatedcommunity_melt=melt(simulatedcommunity_df,id=c("r"))
simulatedcommunity_melt$time=as.numeric(simulatedcommunity_melt$variable)
fig3=ggplot(simulatedcommunity_melt, aes(x=time, y=value, col=r))+geom_point()+theme_classic()+xlab("Time")+ylab("Simulated sp Biomass")

communitybiomass=colSums(simulatedcommunity)
stanDat_null <- list(b = communitybiomass,ag=seq(1,50,by=1),
                      N = length(communitybiomass), newage=seq(1,50,by=1))

Fit_null<- stan(file = "Null_reservessimulated.stan", data = stanDat_null, chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
windows()
pairs(Fit_null, pars=c("r","bmin","B0"))

#model fit
Fit_null_summary <- summary(Fit_null,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_null=as.data.frame(Fit_null_summary$summary)

#plot median reserve trajectory with 50 and 90% credible intervals
row.names(output_Fit_null)
newpred_null=output_Fit_null[55:104,]
colnames(newpred_null)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
newpred_null$Time=seq(1,50,by=1)
fig4=ggplot(NULL)+geom_point(aes(x=seq(1,50,by=1), y=communitybiomass), col="grey",alpha=0.5)+theme_classic()+xlab("Time")+ylab("Community Biomass")+
  geom_line(aes(y=newpred_null$median, x=newpred_null$Time), col="navyblue",lwd=3)+geom_ribbon(aes(x=newpred_null$Time,ymin=newpred_null$`5%`, ymax=newpred_null$`95%`), alpha=0.3, fill="blue")+geom_ribbon(aes(x=newpred_null$Time,ymin=newpred_null$`25%`, ymax=newpred_null$`75%`), alpha=0.3, fill="navyblue")+
  geom_text(aes(x=mean(seq(2,50,by=1)),y=500,label=paste("r=",round(output_Fit_null[2,1],3),"(",round(output_Fit_null[2,4],3),",",round(output_Fit_null[2,8],3),")")))
windows()
ggarrange(fig1,fig2,fig3,fig4, nrow=1,ncol=4, labels=c("a","b","c","d"))

#save.image(file='Zamborain-Masonetal2020_Revised_sp.RData')
#load('Zamborain-Masonetal2020_Revised_sp.RData')
