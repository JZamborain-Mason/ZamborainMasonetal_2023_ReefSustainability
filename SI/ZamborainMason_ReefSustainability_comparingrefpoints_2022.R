# Title: Reef sustainability_SI1_A
# Author: Jessica Zamborain Mason
# Description: This script provides a comparison of  our estimates to those in McClanahan and Graham 2015 and MacNeil et al. 2015
# Please ask for data and permission to use the data to the relevant authors
# R version 3.5.3 (2019-03-11)


#clear R
rm(list=ls())


#set working directory
setwd("")

#load required libraries
library(dplyr)
library(reshape)
library(plyr)
library(ggplot2)
library(rstan)
library(ggpubr)

#my data
JZMdata<-read.csv("reefscale_data_submitted.csv",header=T)
JZMdata$model_component<- ifelse(JZMdata$section=="benchmark"&JZMdata$reserves==1,1,ifelse(JZMdata$section=="benchmark"&JZMdata$remote_20h==1 & JZMdata$remote_20h_inhabited==0,2,3))
myreserves<-JZMdata[JZMdata$model_component==1,]
myremote<-JZMdata[JZMdata$model_component==2,]

#import McClanahan and Graham's 2015 data shared to me by Tim McClanahan
 Tdata<-read.csv("PRSBrawdata.csv", head=T)
 colnames(Tdata)
 summary(as.factor(Tdata$Management))
 #filter MG data to only include high compliance marine reserves
 Tdata<-Tdata[Tdata$Management=="High compliance closure",]
 Tdata=droplevels(Tdata)
 #trends reported by them: 
 mg1<-ggplot(data=Tdata, aes(x=Closure.age, y=Total.biomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Total Biomass (t/km2)")+ggtitle("McClanahan & Graham 2015 Original")+ xlab ("Closure age")


#import MacNeil et al. 2015 data (and code used) shared to me by Aaron MacNeil
  Adata<-read.csv("B0_master.csv", head=T)
  # Create species/functional group table
  Adata$species<- paste(Adata$Genus,Adata$Sp,sep="_")
  Adata$fish_habitat<- rep("Unknown",length(Adata$species))
  Adata$FG<- "UNK"
  Adata$TL<- 0
  # Select only reef-associated species
  fishdata<- read.csv("B0_spp_list_master.csv")
  # Filter out sites with excessive schools
  out<-unique(as.vector(Adata$Site[Adata$Abundance>40000&Adata$Provider=="Newman"]))
  out<- append(out,"Rib Sheltered 3")
  Adata<- Adata[is.na(match(Adata$Site,out)),]
  # Species list
  spx<- unique(Adata$species)
  nspx<- length(spx)
  for (i in 1:nspx){
  print(spx[i])
  Adata$species[Adata$species==spx[i]]<- as.character(fishdata$true_spp[fishdata$species==spx[i]])
  print(length(Adata$species[Adata$species==spx[i]]))
  Adata$fish_habitat[Adata$species==spx[i]] = as.character(fishdata$habitat[fishdata$species==spx[i]])
  Adata$FG[Adata$species==spx[i]]<- as.character(fishdata$fg[fishdata$species==spx[i]])
  Adata$TL[Adata$species==spx[i]]<- as.character(fishdata$trophic_level[fishdata$species==spx[i]])
  }
  Adata<-Adata[Adata$fish_habitat=="Reef",]
  sum(Adata$Biomass==0)
  Adata$lBiomass<- log(Adata$Biomass)
  Adata$BodySize<- Adata$Biomass/Adata$Abundance

  #Aarons data to calculate recovery rates and unfished biomass
  aaronsites<-read.csv("B0data.csv", head=T)
  aaronsites$Site<-aaronsites$site
  aaronsites$Aaronlbiomass<-aaronsites$lbiomass
  #trends reported by them: 
  m1<-ggplot(data=aaronsites[aaronsites$protection=="Closed",], aes(x=mpage, y=Aaronlbiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Total Biomass (log(t/km2))")+ggtitle("MacNeil et al. 2015 Original")+ xlab ("Closure age")
  m1_exp<-ggplot(data=aaronsites[aaronsites$protection=="Closed",], aes(x=mpage, y=exp(Aaronlbiomass)/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Total Biomass (t/km2)")+ggtitle("MacNeil et al. 2015 Original (arithmetic)")+ xlab ("Closure age")

windows()
ggarrange(mg1,m1,m1_exp,nrow=1,ncol=3)

#what if we fit our null model (allowing for exports) to their data, regardless of all the other assumptions?
aaronsreserves<-aaronsites[aaronsites$protection=="Closed" &!is.na(aaronsites$lbiomass),]
aaronsremote<-aaronsites[aaronsites$protection=="Remote" &!is.na(aaronsites$lbiomass),]

#mcclanahan and graham 
#fit null model of only reserves to the McClanahan et al data only using reserves (note that we omit the fact that some reserves were sampled multiple times in time as McClanahan and Graham)
stan_mg<- list(b = log(Tdata$Total.biomass/10),ag=Tdata$Closure.age,
                            res=nrow(Tdata))

stancode_null_onlyreserves<-"data {
  int<lower=1> res; //number of reserve  data points
  real b[res]; //Biomass response variable (y axes) (log-transformed)
  //reserves
  int<lower=0> ag[res]; //age reserve (only for reserves)
  
}
parameters {
  real<lower=0> sigma_e; //error sd for biomass reserves
  real log_r;//community growth rate
  real log_bmin; //biomass at reserve age 0
  real log_B0; //unfished biomass 
  real<lower=0,upper=1> p;// proportion of biomass exported
}
transformed parameters {
  vector[res] mu;//mean log-biomass  reserves
  real<lower=0> r;
  real<lower=0> B0;
  real<lower=0> bmin;
  r=exp(log_r); 
  bmin=exp(log_bmin); 
  B0=exp(log_B0); 
  //reserve component
  for (i in 1:res){ 
    mu[i] = log((1-p)*(B0/(1+((B0-bmin)/bmin)*exp(-(r)*ag[i]))));
  }
}
model {
  //priors
  log_r ~ normal (-2,1); //weekly informative prior  biomass growth rate
  sigma_e ~ cauchy(0,1); //prior sd
  log_bmin ~ normal (log(40),1); //weakly informative prior reserve biomass at age 0
  log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
  p~ uniform (0,1); // uninformative prior for export proportion
  //likelihoods  
  for(n in 1:res){
    target += normal_lpdf(b[n] | mu[n], sigma_e);
  }
}
generated quantities {
  real <lower=0> BMMSY; 
  real <lower=0> MMSY;
  BMMSY=B0/2;
  MMSY=((r*B0)/4); 
}

"
writeLines(stancode_null_onlyreserves,"~/null_onlyreserves.stan")
Fit_null_p_mg<- stan("~/null_onlyreserves.stan", data = stan_mg, iter=20000,warmup=19000,chains = 4,control = list(adapt_delta = 0.99,max_treedepth=20))
windows()
ggarrange(plot(Fit_null_p_mg,pars="B0"),plot(Fit_null_p_mg,pars="r"),plot(Fit_null_p_mg,pars="p"),plot(Fit_null_p_mg,pars="MMSY"),ncol=4,nrow=1)

quantile(rstan::extract(Fit_null_p_mg,par="r")$r,c(0.10,0.5,0.9))
quantile(rstan::extract(Fit_null_p_mg,par="p")$p,c(0.10,0.5,0.9))
quantile(rstan::extract(Fit_null_p_mg,par="B0")$B0,c(0.10,0.5,0.90))
quantile(rstan::extract(Fit_null_p_mg,par="MMSY")$MMSY,c(0.10,0.5,0.90))
#posterior contraction for relevant parameters
var_prior_B0<-(1^2)
var_prior_bmin<-(1^2)
var_prior_logr<-(1^2)
var_prior_p<-(sd(runif(4000,0,1))^2)
#posterior contraction 
contraction_null_mg<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_null_p_mg,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_null_p_mg,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_null_p_mg,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_null_p_mg,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
colnames( contraction_null_mg)<-c("posterior contraction","parameter")
print(contraction_null_mg)

#macneil et al.2015
stan_m<- list(b = log(exp(aaronsreserves$lbiomass)/10),ag=aaronsreserves$mpage,
                            res=nrow(aaronsreserves),b2=log(exp(aaronsremote$lbiomass)/10),
                            rem=nrow(aaronsremote))
stancode_null_reservesrem<-"data {
 int<lower=1> res; //number of reserve  data points
 real b[res]; //Biomass response variable (y axes) (log-transformed)
 int<lower=1> rem; //number of remote  data points
 real b2[rem]; //Biomass response variable (y axes) (log-transformed)
//explanatory variables for each component
//reserves
 int<lower=0> ag[res]; //age reserve (only for reserves)
}
parameters {
 real<lower=0> sigma_e; //error sd for biomass reserves
 real<lower=0> sigma_r; //error sd for biomass remote
 real log_r;//community growth rate
 real log_bmin; //biomass at reserve age 0
 real log_B0; //unfished biomass 
 real<lower=0,upper=1> p;// proportion of biomass exported
}
transformed parameters {
 vector[res] mu;//mean log-biomass  reserves
 vector[rem] mu2;//mean log-biomass  remote
  real<lower=0> r;
  real<lower=0> B0;
  real<lower=0> bmin;
  r=exp(log_r); 
  bmin=exp(log_bmin); 
  B0=exp(log_B0);
//reserve component
for (i in 1:res){ 
       mu[i] = log((1-p)*(B0/(1+((B0-bmin)/bmin)*exp(-(r)*ag[i]))));
    }
//remote component
for (i in 1:rem){ 
       mu2[i] = log(B0);
    }
}
model {
//priors
 log_r ~ normal (-2,1); //weekly informative prior  biomass growth rate
 sigma_e ~ cauchy(0,1); //prior sd
 sigma_r ~ cauchy(0,1); //prior sd
 log_bmin ~ normal (log(40),1); //weakly informative prior reserve biomass at age 0
 log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
 p~ uniform (0,1); // uninformative prior for export proportion
//likelihood  
 for(n in 1:res){
            target += normal_lpdf(b[n] | mu[n], sigma_e);
 }
     for(n in 1:rem){
            target += normal_lpdf(b2[n] | mu2[n], sigma_r);
    }
}
generated quantities {
 real <lower=0> BMMSY; 
 real <lower=0> MMSY;
 BMMSY=B0/2;
 MMSY=((r*B0)/4); 
}

"
writeLines(stancode_null_reservesrem,"~/null_reservesrem.stan")
Fit_null_p_m<- stan("~/null_reservesrem.stan", data = stan_m, iter=20000,warmup=19000,chains = 4,control = list(adapt_delta = 0.99,max_treedepth=20))
windows()
ggarrange(plot(Fit_null_p_m,pars="B0"),plot(Fit_null_p_m,pars="r"),plot(Fit_null_p_m,pars="p"),plot(Fit_null_p_m,pars="MMSY"),ncol=4,nrow=1)

quantile(rstan::extract(Fit_null_p_m,par="r")$r,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_m,par="p")$p,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_m,par="B0")$B0,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_m,par="MMSY")$MMSY,c(0.05,0.5,0.95))
#posterior contraction 
contraction_null_m<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_null_p_m,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_null_p_m,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_null_p_m,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_null_p_m,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
colnames( contraction_null_m)<-c("posterior contraction","parameter")
print(contraction_null_m)


####################################checking additional differences in approach
#families included in our study ..............................................
fams<-read.csv("familiesincluded.csv", head=T)
fams<-fams[fams$Chapter.1_nosharks==1,]
fams$Family<-fams$Fish.family

#Exclude families in McClanahan and Graham's data that we did not include
 Tdata$famsBiomass<-Tdata$Total.biomass-(Tdata$Aulostomidae+Tdata$Carcarhinidae+Tdata$Fistularidae+Tdata$Ginglymostomatidae+Tdata$Holocentridae+Tdata$Muraenidae+Tdata$Others+Tdata$Pempheridae+Tdata$Pomacentridae+Tdata$Scombridae+Tdata$Scorpaenidae)
 mg2<-ggplot(data=Tdata, aes(x=Closure.age, y=famsBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Biomass(t/km2)")+ggtitle("McClanahan & Graham 2015 Families restricted")+ xlab ("Closure age")


#Exclude families in MacNeil et al. data that we did not include
 #make family terminology comparable
 Adata$Family<-as.factor(ifelse(Adata$Family=="ACANTHURIDAE","Acanthuridae", ifelse(Adata$Family=="CHAETODONTIDAE","Chaetodontidae", ifelse(Adata$Family=="CIRRHITIDAE","Cirrhitidae", ifelse(Adata$Family=="DIODONTIDAE","Diodontidae" ,ifelse(Adata$Family=="KYPHOSIDAE","Kyphosidae", ifelse(Adata$Family=="LABRIDAE","Labridae",ifelse(Adata$Family=="MONACANTHIDAE","Monacanthidae",ifelse(Adata$Family=="MULLIDAE","Mullidae",ifelse(Adata$Family=="POMACANTHIDAE","Pomacanthidae",ifelse(Adata$Family=="SERRANIDAE","Serranidae", ifelse(Adata$Family=="TETRAODONTIDAE" |Adata$Family=="Tetradontidae","Tetraodontidae",ifelse(Adata$Family=="ZANCLIDAE","Zanclidae", as.character(Adata$Family))))))))))))))
 fams_Adata<-merge( Adata, fams[,c("Family","Common.family.name")],by="Family", all.x=T)
 fams_Adata$siteyear<- paste(fams_Adata$Site,fams_Adata$Year,sep="::")
 notincludedbiomass<-fams_Adata[is.na(fams_Adata$Common.family.name),]
 data2<-ddply(notincludedbiomass, .(siteyear), summarize, Atoll = Atoll[1], Depth = Depth[1],Latitude=Latitude[1], Longitude=Longitude[1], TideRange=TideRange[1],Year=Year[1],Country=Country[1], Region=Region[1], Protection=Protection[1], Complexity=Complexity[1], Complex_method=Complex_method[1], HardCoral=HardCoral[1], AlgalCover=AlgalCover[1], ReefHabitat=ReefHabitat[1], Biomass=sum(Biomass), Abundance=sum(Abundance), SampMethod=SampMethod[1], SampArea=SampArea[1],  N_reps=N_reps[1], MPAge=MPAge[1], MPSize=MPSize[1], Provider=Provider[1],Site=Site[1]) 
 #merge with all data to elimante the biomass of non required families
 data2<-merge(data2,aaronsites[,c("siteyr", "productivity", "Aaronlbiomass")], by.x="siteyear", by.y="siteyr",all.x=T)
 data2$famBiomass=exp(data2$Aaronlbiomass)-data2$Biomass
 m2_exp<-ggplot(data=data2[data2$Protection=="Closed",], aes(x=MPAge, y=famBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Biomass (t/km2)")+ggtitle("MacNeil et al. 2015 Families restricted")+ xlab ("Closure age")

#comparison of remote loctaions included
ggplot()+geom_boxplot(data=data2[data2$Protection=="Remote",], aes(x=reorder(as.factor(Country),log(famBiomass/10),na.rm=T), y=log(famBiomass/10)),fill="blue",alpha=0.5)+
geom_boxplot(data=myremote, aes(x=as.factor(Locality), y=log(FamBiomass_tkm2)),fill="black",alpha=0.5)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")+ylab("log(Biomass(t/km2)")

#only tropical sites..........................................................
 Tdata<-filter(Tdata, abs(Tdata$LAT)< 23)
 mg3<-ggplot(data=Tdata, aes(x=Closure.age, y=famsBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Biomass(t/km2)")+ggtitle("McClanahan & Graham 2015 Tropical sites")+ xlab ("Closure age")

 data2<-filter(data2,abs(data2$Latitude)< 23)
 m3_exp<-ggplot(data=data2[data2$Protection=="Closed",], aes(x=MPAge, y=famBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Biomass (t/km2)")+ggtitle("MacNeil et al. 2015 Tropical sites")+ xlab ("Closure age")


#McClanahan and Graham's also:................................................
# (i)mixed time-series and space-for time substitution
#(ii)included sites without coral cover information 
#for the reference points we used data with coral cover (following MacNeil et al. approach)
 #note different locations respond differently to protection (not all assymptote after 10-15 years)
 mg4<-ggplot(data=Tdata, aes(x=Closure.age, y=famsBiomass/10,col=Tdata$Location))+guides(col=F)+geom_point()+theme_classic()+geom_smooth(method="lm",level = 0.5)+ylab("Biomass (t/km2)")+ggtitle("McClanahan & Graham 2015 potential process error")
 ggplot(data=Tdata[!is.na(Tdata$Coral.cover),], aes(x=Closure.age, y=famsBiomass,col=Tdata$Location[!is.na(Tdata$Coral.cover)]))+guides(col=F)+geom_point()+theme_classic()+geom_smooth(method="lm",level = 0.5)+ylab("Biomass (t/km2)")+ggtitle("McClanahan & Graham 2015 sites with coral cover")

#only one reserve in mcclanahan and Graham (longest information)
sites<- ddply(Tdata,.(Location),summarize, reefsites=length(Total.biomass)) 
mombasa<-Tdata[Tdata$Location=="Mombasa",]

mombasa_trend<-ggplot(data=mombasa, aes(x=Closure.age, y=famsBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Total Biomass (t/km2)")+ggtitle("Mombasa")+ xlab ("Closure age")
stan_mombasa<- list(b = log(mombasa$famsBiomass/10),ag=mombasa$Closure.age,
                            res=nrow(mombasa))
Fit_null_p_mombasa<- stan("~/null_onlyreserves.stan", data = stan_mombasa, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999))

windows()
ggarrange(mombasa_trend,plot(Fit_null_p_mombasa,pars="B0"),plot(Fit_null_p_mombasa,pars="r"),plot(Fit_null_p_mombasa,pars="p"),plot(Fit_null_p_mombasa,pars="MMSY"),ncol=5,nrow=1)
quantile(rstan::extract(Fit_null_p_mombasa,par="r")$r,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_mombasa,par="p")$p,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_mombasa,par="B0")$B0,c(0.05,0.5,0.95))
quantile(rstan::extract(Fit_null_p_mombasa,par="MMSY")$MMSY,c(0.05,0.5,0.95))

#posterior contraction 
contraction_null_mombasa<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_null_p_mombasa,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_null_p_mombasa,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_null_p_mombasa,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_null_p_mombasa,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
colnames( contraction_null_mombasa)<-c("posterior contraction","parameter")
print(contraction_null_mombasa)
#non-identified p because b0 not bounded


#now including our remote locations to bound B0................................
stan_mombasa_rem<- list(b = log(mombasa$famsBiomass/10),ag=mombasa$Closure.age,
                            res=nrow(mombasa),b2=log(JZMdata$FamBiomass_tkm2[JZMdata$model_component==2]),
               rem=length(log(JZMdata$FamBiomass_tkm2[JZMdata$model_component==2])))

Fit_null_p_mombasa_rem<- stan("~/null_reservesrem.stan", data = stan_mombasa_rem, iter=10000,warmup=9000,chains = 4)

mombasa_trend2<-ggplot(data=mombasa, aes(x=Closure.age, y=famsBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Total Biomass (t/km2)")+ggtitle("Mombasa + our remote")+ xlab ("Closure age")

windows()
ggarrange(mombasa_trend2,plot(Fit_null_p_mombasa_rem,pars="B0"),plot(Fit_null_p_mombasa_rem,pars="r"),plot(Fit_null_p_mombasa_rem,pars="p"),plot(Fit_null_p_mombasa_rem,pars="MMSY"),ncol=5,nrow=1)
median(rstan::extract(Fit_null_p_mombasa_rem,par="r")$r)
median(rstan::extract(Fit_null_p_mombasa_rem,par="p")$p)
median(rstan::extract(Fit_null_p_mombasa_rem,par="B0")$B0)
median(rstan::extract(Fit_null_p_mombasa_rem,par="MMSY")$MMSY)

#posterior contraction 
contraction_null_mombasa_rem<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_null_p_mombasa_rem,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_null_p_mombasa_rem,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_null_p_mombasa_rem,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_null_p_mombasa_rem,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
colnames( contraction_null_mombasa_rem)<-c("posterior contraction","parameter")
print(contraction_null_mombasa_rem)
#now it is identified better


#due to the potential process error: we only used space for time-substitution 
#example (with mcclanahan and Graham)
ggplot(data=Tdata, aes(x=Closure.age, y=famsBiomass,col=Tdata$Location))+geom_point()+theme_classic()+geom_smooth(method="lm")+ylim(c(0,2000))+ylab("Biomass (kg/ha")
Tdata2<-Tdata[!is.na(Tdata$Coral.cover),]

subsampled_data<- ddply(Tdata2,.(Location),
                        function(x) {
                          x[sample(nrow(x),size=1),]
                        })

mg5<-ggplot(data=subsampled_data, aes(x=Closure.age, y=famsBiomass/10))+geom_point()+theme_classic()+geom_smooth(method="gam")+ylab("Biomass(t/km2)")+ggtitle("McClanahan & Graham 2015 methods as our study")+ xlab ("Closure age")

#This all this is without other model assumptions!
windows()
ggarrange(mg1, mg2,mg3, mg4,mg5,mombasa_trend,nrow=1,ncol=6)
windows()
ggarrange(m1_exp, m2_exp,m3_exp, nrow=1,ncol=3)


#all together
ggplot()+geom_point(data=subsampled_data, aes(x=Closure.age, y=famsBiomass/10),pch=21,fill="red",alpha=0.5)+geom_smooth(data=subsampled_data, aes(x=Closure.age, y=famsBiomass/10),col="red",fill="red",method="gam")+
geom_point(data=data2[data2$Protection=="Closed",], aes(x=MPAge, y=famBiomass/10),pch=21,fill="blue",alpha=0.5)+geom_smooth(data=data2[data2$Protection=="Closed",], aes(x=MPAge, y=famBiomass/10),col="blue",fill="blue",method="gam")+
geom_point(data=myreserves, aes(x=Closure.age, y=FamBiomass_tkm2),pch=21,fill="black",alpha=0.5)+geom_smooth(data=myreserves, aes(x=Closure.age, y=FamBiomass_tkm2),col="black",fill="grey",method="gam")+ylab("Biomass(t/km2)")+ xlab ("Closure age")+theme_classic()

################################################################################
#save.image("comparing_refpoints.RData")
#load(file="comparing_refpoints.RData")