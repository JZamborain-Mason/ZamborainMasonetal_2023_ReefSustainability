# Title: Reef sustainability_SI2
# Author: Jessica Zamborain Mason
# Description: This code implements the Bayesian workflow for "Sustainable reference points for multispecies coral reef fisheries"
# R version 3.5.3 (2019-03-11)

#the first part of the code is the same as the original analyses (setting up data for model). However, then I implement the Bayesian workflow analyses to show:
 #(i)  Models produce non-biased results
 #(ii) Model used produces identifiable parameters

#Note that overall this simulation takes a very long time to run!

#clear R
 rm(list=ls())


#set working directory (where the data and code of this repository is stored)
 setwd("")
 

#load required libraries
 library(ggrepel)  #Kamil Slowikowski (2018). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.8.0. https://CRAN.R-project.org/package=ggrepel
 library(ggplot2) #H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
 library(ggpubr) #  Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.2.https://CRAN.R-project.org/package=ggpubr
 library(plyr) #Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
 library(rstan) #Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. http://mc-stan.org/.
 library(lme4) #  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
 library(rstan)#Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. http://mc-stan.org/.
 library(dplyr)#Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2019). dplyr: A Grammar of Data Manipulation. R package version 0.8.0.1. https://CRAN.R-project.org/package=dplyr
 library(rworldmap)#South, Andy 2011 rworldmap: A New R package for Mapping Global Data. The R Journal Vol. 3/1 : 35-43.
 library(rworldxtra)#  Andy South (2012). rworldxtra: Country boundaries at high resolution.. R package version 1.01. https://CRAN.R-project.org/package=rworldxtra
 library(mgcv)#Wood, S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman and Hall/CRC
 library(mice)#Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. URL https://www.jstatsoft.org/v45/i03/
 library(VIM)#Alexander Kowarik, Matthias Templ (2016). Imputation with the R Package VIM. Journal of Statistical Software, 74(7), 1-16. doi:10.18637/jss.v074.i07
 library(coefplot)#Jared P. Lander (2018). coefplot: Plots Coefficients from Fitted Models. R package version 1.2.6. https://CRAN.R-project.org/package=coefplot
 library(bayesplot)#Jonah Gabry and Tristan Mahr (2018). bayesplot: Plotting for Bayesian Models. R package version 1.6.0. https://CRAN.R-project.org/package=bayesplot
 library(reshape)#H. Wickham. Reshaping data with the reshape package. Journal of Statistical Software, 21(12), 2007.
 library(tidyr)#Hadley Wickham and Lionel Henry (2019). tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.3. https://CRAN.R-project.org/package=tidy
 library(car)#John Fox and Sanford Weisberg (2011). An {R} Companion to Applied Regression, Second Edition. Thousand Oaks CA: Sage. URL:http://socserv.socsci.mcmaster.ca/jfox/Books/Companion
 library(tidyverse) #Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. R package version 1.2.1.https://CRAN.R-project.org/package=tidyverse
 library(stringr)#Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0.https://CRAN.R-project.org/package=string
 library(ggridges)#  Claus O. Wilke (2018). ggridges: Ridgeline Plots in 'ggplot2'. R package version 0.5.1. https://CRAN.R-project.org/package=ggridges
 library(rstanarm)#Goodrich B, Gabry J, Ali I & Brilleman S. (2018). rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17.4. http://mc-stan.org/.
 library(loo)#Vehtari A, Gabry J, Yao Y, Gelman A (2019). "loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models." R package version 2.1.0, <URL: https://CRAN.R-project.org/package=loo>.
 library(geoR)#Paulo J. Ribeiro Jr and Peter J. Diggle (2018). geoR: Analysis of Geostatistical Data. R package version 1.7-5.2.1. https://CRAN.R-project.org/package=geoR
 library(DHARMa)#Florian Hartig (2019). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.2.4. https://CRAN.R-project.org/package=DHARMa
 library(brms)# Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of StatisticalSoftware, 80(1), 1-28. doi:10.18637/jss.v080.i01
 library(ggeffects)#Lüdecke D (2018). "ggeffects: Tidy Data Frames of Marginal Effects from Regression Models." _Journal of Open SourceSoftware_, *3*(26), 772. doi: 10.21105/joss.00772 (URL: http://doi.org/10.21105/joss.00772).
 library(extraDistr)#  Tymoteusz Wolodzko (2020). extraDistr: Additional Univariate and Multivariate Distributions. R package version 1.9.1.
 library(Matching)#Jasjeet S. Sekhon (2011). Multivariate and Propensity Score Matching Software with Automated Balance Optimization: The Matching Package for R.Journal of Statistical Software, 42(7), 1-52. <doi:10.18637/jss.v042.i07>
 
#to run chains in parallel
 rstan_options(auto_write = T)
 options(mc.cores = parallel::detectCores())
 #backup_options <- options()
 #options(backup_options)

#upload data (data is available in the paper's supporting information)
 reefscale_data<- read.csv("reefscale_data_submitted.csv", head=T)
 jurisdictionscale_data<- read.csv("jurisdictionscale_data_submitted.csv", head=T)


#functions used in this script.........................................................................................................

 #mean center covariates: standardise (Following Gelman and Hill 2007)
 standardise <- function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 
 #Function to check missingness
 pMiss <- function(x){sum(is.na(x))/length(x)*50}

 #functions from Vehtari et al  2019 to monitor convergence
 #if you use these documents please see what is required for the authors (https://github.com/avehtari/rhat_ess)
 #and cite the paper: Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, Paul-Christian Bürkner (2019): Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
 source('monitornew_Vehtarietal2019.R')
 source('monitorplot_Vehtarietal2019.R')

#######################################################################################################################################
#clean reefscale and jurisdiction scale data for analyses..............................................................................
#standardize and relevel categorical variables
 reefscale_data$indexr<- ifelse(reefscale_data$Region=="C",1,ifelse(reefscale_data$Region=="I",2,3))
 reefscale_data$Region<- as.factor(reefscale_data$Region)
 reefscale_data$Atoll<- ifelse(reefscale_data$Atoll=="1",1,0)
 reefscale_data<- reefscale_data[!reefscale_data$ReefHabitat=="Bank",]
 reefscale_data<- droplevels(reefscale_data)
 reefscale_data<- reefscale_data %>%
  mutate(SampMethod=relevel(as.factor(SampMethod), ref="Standard belt transect"),
         ReefHabitat=relevel(as.factor(ReefHabitat),ref="Slope"),
         sClosure.size=standardise(log(Closure.size)),
         sSampArea=standardise(log(reefscale_data$SampArea)),
         sDepth=standardise(sqrt(Depth)),
         sOcean_prod=standardise(log(Prod_mgCm2d_finer)),
         sHardCoral=standardise(sqrt(HardCoral)),
         sSST=standardise(SST_dC),
         sSSTanom=standardise(SSTanom_dC),
         sgravtot=log(Grav_tot+1),
         stt_market=standardise(log(TT_market_h)))

#create dummy varibles for the different categorical variables
 reefscale_data$pointintercept<- ifelse(reefscale_data$SampMethod=="Point intercept", 1,0)
 reefscale_data$backreef<- ifelse(reefscale_data$ReefHabitat=="Lagoon_Back reef", 1,0)
 reefscale_data$crest<- ifelse(reefscale_data$ReefHabitat=="Crest", 1,0)
 reefscale_data$flat<- ifelse(reefscale_data$ReefHabitat=="Flat", 1,0)
 reefscale_data$distancesampling<- ifelse(reefscale_data$SampMethod=="Distance sampling",1,0)
 reefscale_data$model_component<- ifelse(reefscale_data$section=="benchmark"&reefscale_data$reserves==1,1,ifelse(reefscale_data$section=="benchmark"&reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,2,3))
 
#first read Aarons data to reclaify those he had defined as remote but that were <20h from human settlements
 Adata<-read.csv("B0data_MacNeiletal.csv",header=T)
 reefscale_data$siteyear<-paste(reefscale_data$Site, reefscale_data$Year,sep="::")
 reefscale_data<-merge(reefscale_data,Adata[,c("siteyr","MGMT_type")],by.x="siteyear",by.y="siteyr",all.x=T)
 reefscale_data$definedprotection<- ifelse(reefscale_data$section=="benchmark"& reefscale_data$reserves==1,"HC Marine reserves", ifelse(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,"Remote",ifelse(reefscale_data$Management=="Restricted"|(reefscale_data$Management=="Remote" & !(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0)&reefscale_data$MGMT_type=="Access restricted"),"Restricted","Openly fished")))
 
 #MacNeil et al. 2015 data defined as remote that does not satisfy our remote is classified as restricted (overlaps with some data providor classifications) 
 reefscale_data$status_restricted<- ifelse(reefscale_data$model_component==3 &reefscale_data$definedprotection=="Restricted",1,0)
 reefscale_data$status_fished<- ifelse(reefscale_data$model_component==3 &reefscale_data$definedprotection=="Fished",1,0)
 reefscale_data$status_management<- ifelse(reefscale_data$status_restricted==1,"Restricted",ifelse(reefscale_data$status_fished==1,"Fished", reefscale_data$definedprotection))
 
 #we do not have reserve size/coral cover for most sites so we  create variable to add 0 to those components (coral cover for the post-model estimation: average conditions)
 reefscale_data$sHardCoral2<- ifelse(is.na(reefscale_data$sHardCoral),0,reefscale_data$sHardCoral)
 reefscale_data$sClosure.size2<- ifelse(is.na(reefscale_data$sClosure.size),0,reefscale_data$sClosure.size)
 reefscale_data$Closure.age2<- ifelse(is.na(reefscale_data$Closure.age),0,reefscale_data$Closure.age)
 
#MacNeil et al. 2015 data defined as remote that does not satisfy our remote is classified as restricted (overlaps with some data providor classifications) 
 reefscale_data$status_restricted<- ifelse(reefscale_data$model_component==3 &(reefscale_data$Management=="Restricted"|reefscale_data$Management=="Remote"),1,0)
 reefscale_data$status_fished<- ifelse(reefscale_data$model_component==3 &reefscale_data$Management=="Fished",1,0)
 reefscale_data$status_management<- ifelse(reefscale_data$status_restricted==1,"Restricted",ifelse(reefscale_data$status_fished==1,"Fished", reefscale_data$definedprotection))

#we do not have reserve size/coral cover for most sites so we  create variable to add 0 to those components (coral cover for the post-model estimation: average conditions)
 reefscale_data$sHardCoral2<- ifelse(is.na(reefscale_data$sHardCoral),0,reefscale_data$sHardCoral)
 reefscale_data$sClosure.size2<- ifelse(is.na(reefscale_data$sClosure.size),0,reefscale_data$sClosure.size)
 reefscale_data$Closure.age2<- ifelse(is.na(reefscale_data$Closure.age),0,reefscale_data$Closure.age)


#separate reserves, remote and fished data and add index for random effects
 fished_data<- reefscale_data[reefscale_data$model_component==3,]
 fished_data<- droplevels(fished_data)
 reserves_complete<- reefscale_data[reefscale_data$model_component==1,]
 reserves_complete<- droplevels(reserves_complete)
 remote_complete<- reefscale_data[reefscale_data$model_component==2,]
 remote_complete<- droplevels(remote_complete)

#index for random effects 
 fished_data$indexj<- as.numeric(as.factor(fished_data$Larger))
 remote_complete$indexj<- as.numeric(as.factor(remote_complete$Larger))
 reserves_complete$indexj<- base::as.numeric(as.factor(reserves_complete$Larger))

#create new data to plot model predictions
#As we have standardized all continuous variables, we create newdata  based only on reserve age (i.e., at average conditions)
 newdata<- with(reserves_complete, data.frame(Closure.age=seq(min(Closure.age), max(Closure.age), len=length(reserves_complete$Closure.age)))) 
 Bio<- seq(0,350,1) #biomass data for surplus production curve
 B<- length(Bio)

#######################################################################################################################################
## Run reference point and status model  ..............................................................................................

#data
 stanDat_full_country<- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                             res=nrow(reserves_complete), d=reserves_complete$sDepth, hc=reserves_complete$sHardCoral, si=reserves_complete$sClosure.size,
                             at=reserves_complete$Atoll,sst=reserves_complete$sSST,
                             op=reserves_complete$sOcean_prod,rh_c=reserves_complete$crest,rh_b=reserves_complete$backreef,rh_f=reserves_complete$flat,
                             cm_pc=reserves_complete$pointintercept, sa=reserves_complete$sSampArea,
                             gr=reserves_complete$sgravtot,
                             b2 = log(remote_complete$FamBiomass_tkm2),rem=nrow(remote_complete), d2=remote_complete$sDepth, hc2=remote_complete$sHardCoral, 
                             at2=remote_complete$Atoll,sst2=remote_complete$sSST,
                             op2=remote_complete$sOcean_prod,rh_c2=remote_complete$crest,rh_b2=remote_complete$backreef,
                             sa2=remote_complete$sSampArea,
                             cm_ds2=remote_complete$distancesampling,
                             b3 = log(fished_data$FamBiomass_tkm2),fis=nrow(fished_data), d3=fished_data$sDepth, 
                             at3=fished_data$Atoll,sst3=fished_data$sSST,
                             op3=fished_data$sOcean_prod,rh_c3=fished_data$crest,rh_b3=fished_data$backreef,rh_f3=fished_data$flat,
                             cm_pc3=fished_data$pointintercept, sa3=fished_data$sSampArea,cm_ds3=fished_data$distancesampling,hc3=fished_data$sHardCoral2,
                             m_r3=fished_data$status_restricted,gr3=fished_data$sgravtot,
                             R=nlevels(as.factor(reserves_complete$Larger)),
                             pr=reserves_complete$indexj,
                             R2=nlevels(as.factor(remote_complete$Larger)),
                             pr2=remote_complete$indexj,
                             R3=nlevels(as.factor(fished_data$Larger)),
                             pr3=fished_data$indexj,
                             newage=newdata$Closure.age,
                             Bio=Bio, B=B)
 
 
 #null model 
 stancode_null="data {
 int<lower=1> res; //number of reserve  data points
 real b[res]; //Biomass response variable (y axes) (log-transformed)
 int<lower=1> rem; //number of remote  data points
 real b2[rem]; //Biomass response variable (y axes) (log-transformed)
 int<lower=1> fis; //number of fished  data points
 real b3[fis]; //Biomass response variable (y axes) (log-transformed)
//explanatory variables for each component
//reserves
 int<lower=0> ag[res]; //age reserve (only for reserves)
 
}
parameters {
 real<lower=0> sigma_e; //error sd for biomass reserves
 real<lower=0> sigma_r; //error sd for biomass remote
 real<lower=0> sigma_f; //error sd for biomass fished
 real log_r;//community biomass intrinsic growth rate (log)
 real log_bmin; //biomass at reserve age 0(log)
 real log_B0; //unfished biomass (log)
 real I_fished; //intercept for fished reefs
}
transformed parameters {
 vector[res] mu;//mean log-biomass  reserves
 vector[rem] mu2;//mean log-biomass  remote
 vector[fis] mu3;//mean log-biomass  fished
 real<lower=0> r;//community biomass intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass 
 r=exp(log_r); 
 bmin=exp(log_bmin); 
 B0=exp(log_B0); 
//reserve component
for (i in 1:res){ 
       mu[i] = log(B0*exp(log(bmin/B0)*exp(-r*ag[i])));
    }
//remote component
for (i in 1:rem){ 
       mu2[i] = log(B0);
    }
//fished component
for (i in 1:fis){ 
     mu3[i] =I_fished;
 }
}
model {
//priors
 log_r ~ normal (-2,1); //weekly informative prior  biomass intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 sigma_r ~ cauchy(0,1); //uninformative prior sd
 sigma_f ~ cauchy(0,1); //uninformative prior sd
 log_bmin ~ normal (log(10),1); //weakly informative prior reserve biomass at age 0
 log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
 I_fished~ normal(5,5); // prior for intercept in fished reefs
//likelihoods  
 for(n in 1:res){
            target += normal_lpdf(b[n] | mu[n], sigma_e);
    }
 for(n in 1:rem){
            target += normal_lpdf(b2[n] | mu2[n], sigma_r);
    }
 for(n in 1:fis){
            target += normal_lpdf(b3[n] | mu3[n], sigma_f);
   }
}
generated quantities {
 real <lower=0> BMMSY; //BMMSY for averge and most common environemental conditions
 real <lower=0> MMSY;//MMSY for average and most common environmental conditions
//average and common conditions
 BMMSY=B0/2.718281828;
 MMSY=((r*B0)/2.718281828);
}


"
 writeLines(stancode_null,"~/null.stan")
 Fit_null<- stan("~/null.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
 windows()
 pairs(Fit_null,pars=c("log_B0","log_r","log_bmin"))
 #posterior contraction for our model parameters ref points
 #posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
 var_prior_logB0<-(1^2)
 var_prior_logbmin<-(1^2)
 var_prior_logr<-(1^2)
var_prior_MMSY<-(sd(exp(rnorm(4000,log(120),1))*exp(rnorm(4000,-2,1))/2.718281828)^2)
var_prior_BMMSY<-(sd(exp(rnorm(4000,log(120),1))/2.718281828)^2)

 #posterior contraction for full model reference points
  (var_prior_MMSY-(sd(rstan::extract(Fit_null,pars=c("MMSY"))$MMSY)^2))/var_prior_MMSY
 (var_prior_BMMSY-(sd(rstan::extract(Fit_null,pars=c("BMMSY"))$BMMSY)^2))/var_prior_BMMSY

 #posterior contraction for our data
 contraction_null<-as.data.frame(cbind(round(c((var_prior_logB0-(sd(rstan::extract(Fit_null,pars=c("log_B0"))$log_B0)^2))/var_prior_logB0,(var_prior_logr-(sd(rstan::extract(Fit_null,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_logbmin-(sd(rstan::extract(Fit_null,pars=c("log_bmin"))$log_bmin)^2))/var_prior_logbmin),2),c("log_B0","log_r","log_bmin")))
 colnames( contraction_null)<-c("posterior contraction","parameter")
 print(contraction_null)
 
### BAYESIAN WORKFLOW ###  
#1. prior predictive checks.....................................................
 nsims<-50
 #sample from priors
 B0_sim<-exp(rnorm(nsims,log(120),1))
 bmin_sim<-exp(rnorm(nsims,log(10),1))
 r_sim<-exp(rnorm(nsims,-2,1))
 BMMSY_sim<-B0_sim/2.718281828
 MMSY_sim<-(B0_sim*r_sim)/2.718281828
 sigma_e_sim<-rhcauchy(nsims,1)
 sigma_r_sim<-rhcauchy(nsims,1)
 sigma_f_sim<-rhcauchy(nsims,1)
 I_fished_sim<-rnorm(nsims,5,5) 
#generate simulated data biomass
 res_sim<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 rem_sim<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 fis_sim<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)

 for (n in 1:nsims){
   for (i in 1:nrow(reserves_complete)){
  res_sim[i,n]<-rnorm(1,log(B0_sim[n]*exp(log(bmin_sim[n]/B0_sim[n])*exp(-r_sim[n]*reserves_complete$Closure.age[i]))),sigma_e_sim[n])
  }
   rem_sim[,n]<-rnorm(nrow(remote_complete),log(B0_sim[n]),sigma_r_sim[n])
   fis_sim[,n]<-rnorm(nrow(fished_data),I_fished_sim[n],sigma_f_sim[n])
 }
 res_sim2<-as.data.frame(res_sim)
 res_sim2$Closure.age<-reserves_complete$Closure.age
 res_sim2 <- melt(res_sim2 ,  id.vars = 'Closure.age')
 
 windows()
 a<-ggplot()+geom_density(aes(x=exp(melt(res_sim)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in reserves (t/km2)")
 b<-ggplot() + geom_point(data=res_sim2,aes(Closure.age,exp(value),group = variable),alpha=0.1,col="grey")+ geom_smooth(data=res_sim2,method="gam",aes(Closure.age,exp(value),group = variable),alpha=0.1,lty=2,lwd=0.1,col="grey")+
    geom_smooth(aes(x=reserves_complete$Closure.age,y=exp(matrixStats::rowMeans2(res_sim))),col="black",lwd=2,method="gam")+ylim(c(0,500))+ggtitle("")+guides(col=F)+ylab("Simulated biomass in reserves (t/km2)")+theme_classic()
 c<- ggplot()+geom_density(aes(x=exp(melt(rem_sim)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+ggtitle("")+xlab("Simulated biomass in remote (t/km2)")
 d<- ggplot()+geom_density(aes(x=exp(melt(fis_sim)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in fished (t/km2)")
 prior_pred_fig<-ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d")) 
 annotate_figure(prior_pred_fig,top="Prior predictive check")
 
#2.Computational faithfulness: simulation-based calibration.....................
 sbc_null_B0<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_bmin<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_MMSY<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_BMMSY<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_I_fished<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_sigma_e<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_sigma_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_null_sigma_f<-matrix(NA,nrow=4000,ncol=nsims)
 
 for (n in 1:nsims){
    print(paste("Simulation number ", n, "out of ",max(nsims)))
    stanDat_sbc<- list(b = res_sim[,n],ag=reserves_complete$Closure.age,
                                res=nrow(reserves_complete),
                                b2 = rem_sim[,n],rem=nrow(remote_complete), 
                                b3 = fis_sim[,n],fis=nrow(fished_data))
    fit_sbc<-stan( "~/null.stan", data = stanDat_sbc, chains = 4,control = list(adapt_delta = 0.999))
    sbc_null_B0[,n]<-rstan::extract(fit_sbc,pars=c("B0"))$B0
    sbc_null_bmin[,n]<-rstan::extract(fit_sbc,pars=c("bmin"))$bmin
    sbc_null_r[,n]<-rstan::extract(fit_sbc,pars=c("r"))$r
    sbc_null_MMSY[,n]<-rstan::extract(fit_sbc,pars=c("MMSY"))$MMSY
    sbc_null_BMMSY[,n]<-rstan::extract(fit_sbc,pars=c("BMMSY"))$BMMSY
    sbc_null_I_fished[,n]<-rstan::extract(fit_sbc,pars=c("I_fished"))$I_fished
    sbc_null_sigma_e[,n]<-rstan::extract(fit_sbc,pars=c("sigma_e"))$sigma_e
    sbc_null_sigma_r[,n]<-rstan::extract(fit_sbc,pars=c("sigma_r"))$sigma_r
    sbc_null_sigma_f[,n]<-rstan::extract(fit_sbc,pars=c("sigma_f"))$sigma_f
    
 }

 
a<-ggplot()+geom_histogram(aes(x=melt(log(sbc_null_B0))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(120),1)),fill="red",alpha=0.5)+xlab("log (Unfished biomass, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+
   geom_text(aes(6.2,29000),label="Simulated-data
posterior",col="darkgrey")+
   geom_text(aes(6.5,10000),label="Prior",col="darkred")
b<-ggplot()+geom_histogram(aes(x=melt(log(sbc_null_r))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,-2,1)),fill="red",alpha=0.5)+xlab("log (intrinsic growth rate, 1/t)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
c<-ggplot()+geom_histogram(aes(x=melt(log(sbc_null_bmin))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(10),1)),fill="red",alpha=0.5)+xlab("log (biomass reserve age 0, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
e<-ggplot()+geom_histogram(aes(x=melt(sbc_null_I_fished)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,5,5)),fill="red",alpha=0.5)+xlab("log(intercept biomass in fished reefs, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
f<-ggplot()+geom_histogram(aes(x=log(melt(sbc_null_sigma_e)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd reserve biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
g<-ggplot()+geom_histogram(aes(x=log(melt(sbc_null_sigma_r)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd remote biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
h<-ggplot()+geom_histogram(aes(x=log(melt(sbc_null_sigma_f)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd fished biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
i<-ggplot()+geom_histogram(aes(x=log(melt(sbc_null_BMMSY)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(exp(rnorm(4000*50,log(120),1))/2.718281828)),fill="red",alpha=0.5)+xlab("log BMMSY, t/km2")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
j<-ggplot()+geom_histogram(aes(x=log(melt(sbc_null_MMSY)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(exp(rnorm(4000*50,log(120),1))*exp(rnorm(4000*50,-2,1))/2.718281828)),fill="red",alpha=0.5)+xlab("log MMSY (t/km2/y)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
windows()
sbc_null_fig<-ggarrange(a,b,c,e,f,g,h,i,j,nrow=2,ncol=5)
annotate_figure(sbc_null_fig, top="Simulation-based calibration")


#rank plot for MMSY and BMMSY 
 sbc_rank_MMSY_null<- NA
 sbc_rank_BMMSY_null<- NA
 for (i in 1:nsims) { # Compute SBC rank
    post_MMSY_i <- sbc_null_MMSY[,i]
    post_BMMSY_i <- sbc_null_BMMSY[,i]
    idx_so <- seq(1,length(sbc_null_MMSY[,i]),5) # thin to remove autocorrelation (if wanted)
    sbc_rank_MMSY_null[i] <- sum( MMSY_sim[i] < post_MMSY_i [idx_so] )
    sbc_rank_BMMSY_null[i] <- sum( BMMSY_sim[i] < post_BMMSY_i [idx_so] )
 }
 est_parsMMSY_null <- data.frame(sbc_rank_MMSY_null)
 est_parsBMMSY_null <- data.frame(sbc_rank_BMMSY_null)
 ggarrange(ggplot()+geom_histogram(aes(x=est_parsMMSY_null$sbc_rank_MMSY_null),bins=15,col="black",fill="grey")+theme_classic()+xlab("Rank for MMSY"),ggplot()+geom_histogram(aes(x=est_parsBMMSY_null$sbc_rank_BMMSY_null),bins=15,col="black",fill="grey")+theme_classic()+xlab("Rank for BMMSY"),labels=c("a","b"))

 #test distribution is uniform
 mmsy_null_ks<-ks.boot(sbc_rank_MMSY_null,runif(length(sbc_rank_MMSY_null),min(sbc_rank_MMSY_null),max(sbc_rank_MMSY_null)),alternative = c("two.sided")) 
 bmmsy_null_ks<-ks.boot(sbc_rank_BMMSY_null,runif(length(sbc_rank_BMMSY_null),min(sbc_rank_BMMSY_null),max(sbc_rank_BMMSY_null))) 
#low p value indicates that we reject the null hypothesis that both samples were draw from the same distribution. As it is high, we do not reject it
 
#3.Model sensitivity: z-score vs posterior contraction.........................

#zscore: This measure estimates how close the posterior mean is to the truth relative to the posterior
#standard deviation, in other words how close the entire posterior distribution is to the
#truth value: posterior (posterior mean - prior)/sd posterior
z_score<-matrix(NA,ncol=9,nrow=nsims)
z_score[,1]<-(matrixStats::colMeans2(sbc_null_B0)-B0_sim)/matrixStats::colSds(sbc_null_B0)
z_score[,2]<-(matrixStats::colMeans2(sbc_null_r)-r_sim)/matrixStats::colSds(sbc_null_r)
z_score[,3]<-(matrixStats::colMeans2(sbc_null_bmin)-bmin_sim)/matrixStats::colSds(sbc_null_bmin)
z_score[,4]<-(matrixStats::colMeans2(sbc_null_I_fished)-I_fished_sim)/matrixStats::colSds(sbc_null_I_fished)
z_score[,5]<-(matrixStats::colMeans2(sbc_null_sigma_e)-sigma_e_sim)/matrixStats::colSds(sbc_null_sigma_e)
z_score[,6]<-(matrixStats::colMeans2(sbc_null_sigma_r)-sigma_r_sim)/matrixStats::colSds(sbc_null_sigma_r)
z_score[,7]<-(matrixStats::colMeans2(sbc_null_sigma_f)-sigma_f_sim)/matrixStats::colSds(sbc_null_sigma_f)
z_score[,8]<-(matrixStats::colMeans2(sbc_null_MMSY)-MMSY_sim)/matrixStats::colSds(sbc_null_MMSY)
z_score[,9]<-(matrixStats::colMeans2(sbc_null_BMMSY)-BMMSY_sim)/matrixStats::colSds(sbc_null_BMMSY)
colnames(z_score)<-c("B0","r","bmin","I_fished","sigma_e","sigma_r","sigma_f","MMSY","BMMSY")
#posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
pc<-matrix(NA,ncol=9,nrow=nsims)
var_prior_logB0<-(1^2)
var_prior_logbmin<-(1^2)
var_prior_logr<-(1^2)
var_prior_I_fished<-(5^2)
var_prior_sigma_e<-(sd(rhcauchy(4000,1))^2)
var_prior_sigma_r<-(sd(rhcauchy(4000,1))^2)
var_prior_sigma_f<-(sd(rhcauchy(4000,1))^2)
var_prior_MMSY<-(sd(exp(rnorm(4000,log(120),1))*exp(rnorm(4000,-2,1))/2.718281828)^2)
var_prior_BMMSY<-(sd(exp(rnorm(4000,log(120),1))/2.718281828)^2)


pc[,1]<-(var_prior_logB0-(matrixStats::colSds(log(sbc_null_B0))^2))/var_prior_logB0
pc[,2]<-(var_prior_logr-(matrixStats::colSds(log(sbc_null_r))^2))/var_prior_logr
pc[,3]<-(var_prior_logbmin-(matrixStats::colSds(log(sbc_null_bmin))^2))/var_prior_logbmin
pc[,4]<-(var_prior_I_fished-(matrixStats::colSds(sbc_null_I_fished)^2))/var_prior_I_fished
pc[,5]<-(var_prior_sigma_e-(matrixStats::colSds(sbc_null_sigma_e)^2))/var_prior_sigma_e
pc[,6]<-(var_prior_sigma_r-(matrixStats::colSds(sbc_null_sigma_r)^2))/var_prior_sigma_r
pc[,7]<-(var_prior_sigma_f-(matrixStats::colSds(sbc_null_sigma_f)^2))/var_prior_sigma_f
pc[,8]<-(var_prior_MMSY-(matrixStats::colSds(sbc_null_MMSY)^2))/var_prior_MMSY
pc[,9]<-(var_prior_BMMSY-(matrixStats::colSds(sbc_null_BMMSY)^2))/var_prior_BMMSY

colnames(pc)<-c("B0","r","bmin","I_fished","sigma_e","sigma_r","sigma_f","MMSY","BMMSY")

pc_zscore<-cbind(melt(pc),melt(z_score)$value)
ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
   geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3)+
   geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+
   xlab("Posterior contraction (Nsims=50)")+ylab("Z-score (Nsims=50)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)
 windows()
 pc_zscore_fig2<-ggplot()+geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Median posterior contraction (Nsims=50)")+ylab("Median z-score (Nsims=50)")+theme(panel.background = element_rect(fill="white",colour="black"))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")

 pc_zscore_fig<-ggplot()+geom_point(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Mean posterior contraction (Nsims=50)")+ylab("Mean z-score (Nsims=50)")+theme(panel.background = element_rect(fill="white",colour="black"))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")
 
 zscore_fig<-ggplot()+geom_violin(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2,draw_quantiles = 0.5)+theme(axis.text.x=element_blank(),panel.background = element_rect(fill="white",colour="black"))+ylim(c(5,-5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score")+geom_jitter(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2,pch=21)
 abszscore_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(x=X2,y=abs(melt(z_score)$value),fill=X2),alpha=0.2)+theme(axis.text.x=element_text(angle=90),panel.background = element_rect(fill="white",colour="black"))+ylim(c(0,5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score
(absolute)")+geom_jitter(data=pc_zscore,aes(x=X2,y=abs(melt(z_score)$value),fill=X2),alpha=0.2,pch=21)+
    geom_hline(yintercept = c(3,4),lty=3,col="red")
  zscore_fig<-ggarrange( zscore_fig, abszscore_fig,nrow=2,ncol=1,heights=c(1,0.8))
 pc_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme(panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction")+geom_jitter(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2,pch=21)

 pc_fig_zoom<- ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme(axis.text.y=element_blank(), panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction (zoom)")+geom_jitter(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2,pch=21)+coord_cartesian(xlim=c(0, 1))
  
ggarrange(zscore_fig,pc_fig, pc_fig_zoom, nrow=1,ncol=3,labels = c("a","b",""),widths=c(1,0.8,0.5))
 ggarrange( pc_zscore_fig2,pc_zscore_fig,nrow=1,ncol=2,labels = c("a","b"),widths=c(1,1))
 
#4. Posterior predictive checks.................................................
#for the last simulation
joined_sim <- rstan::extract(fit_sbc)
n_sims <- length(joined_sim $lp__)
y_rep_reserves <- array(NA, c(n_sims, nrow(reserves_complete)))
y_rep_remote <- array(NA, c(n_sims, nrow(remote_complete)))
y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

for (s in 1:n_sims){
   y_rep_reserves[s,] <- rnorm(nrow(reserves_complete), joined_sim$mu[s,], joined_sim$sigma_e[s])
   y_rep_remote[s,] <- rnorm(nrow(remote_complete), joined_sim$mu2[s,], joined_sim$sigma_r[s])
   y_rep_fished[s,] <- rnorm(nrow(fished_data), joined_sim$mu3[s,], joined_sim$sigma_f[s])
}
bayesplot::color_scheme_set(scheme = "gray")
a<- bayesplot::ppc_dens_overlay(res_sim[,50],y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<- bayesplot::ppc_dens_overlay(rem_sim[,50],y_rep_remote[1:4000,])+ggtitle("Remote")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))+guides(col=F)
c<- bayesplot::ppc_dens_overlay(fis_sim[,50],y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))
windows()
ggarrange(a,b,c,nrow=1,ncol=3)


#for all simulations
mu_reserves_sim<-array(NA,dim=c(4000,nrow(reserves_complete),nsims))
mu_remote_sim<-array(NA,dim=c(4000,nrow(remote_complete),nsims))
mu_fished_sim<-array(NA,dim=c(4000,nrow(fished_data),nsims))
y_rep_reserves <- array(NA, c(4000, nrow(reserves_complete),nsims))
y_rep_remote <- array(NA, c(4000, nrow(remote_complete),nsims))
y_rep_fished <- array(NA, c(4000, nrow(fished_data),nsims))


for(j in 1: nsims){
 for (i in 1:nrow(reserves_complete)){
  mu_reserves_sim[,i,j]<-log(sbc_null_B0[,j]*exp(log(sbc_null_bmin[,j]/sbc_null_B0[,j])*exp(-sbc_null_r[,j]*reserves_complete$Closure.age[i]))) 
 }  
 for(i in 1:nrow(remote_complete))  {
  mu_remote_sim [,i,j]<-log(sbc_null_B0[,j])  
 }
 for(i in 1:nrow(fished_data))  {
  mu_fished_sim [,i,j]<-sbc_null_I_fished[,j]  
 }
 for (s in 1:4000){
      y_rep_reserves[s,,j] <- rnorm(nrow(reserves_complete), mu_reserves_sim[s,,j], sbc_null_sigma_e[s,j])
      y_rep_remote[s,,j] <- rnorm(nrow(remote_complete), mu_remote_sim[s,,j], sbc_null_sigma_r[s,j])
      y_rep_fished[s,,j] <- rnorm(nrow(fished_data), mu_fished_sim[s,,j], sbc_null_sigma_f[s,j])
 }
}
windows()
a<-bayesplot::ppc_dens_overlay(res_sim[,50],y_rep_reserves[1:500,,50])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<-bayesplot::ppc_dens_overlay(rem_sim[,50],y_rep_remote[1:500,,50])+ggtitle("remote ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
c<-bayesplot::ppc_dens_overlay(fis_sim[,50],y_rep_fished[1:500,,50])+ggtitle("fished ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,nrow=1,ncol=3)


a2<-ppc_stat(res_sim[,50],y_rep_reserves[,,50],stat="mean")+ggtitle("Reserves")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
b2<-ppc_stat(rem_sim[,50],y_rep_remote[,,50],stat="mean")+ggtitle("Remote")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
c2<-ppc_stat(fis_sim[,50],y_rep_fished[,,50],stat="mean")+ggtitle("Fished")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,a2,b2,c2,nrow=2,ncol=3)

################################################################################


#run full model 

 stancode_full="data {
   int<lower=1> res; //number of reserve  data points
   real b[res]; //Biomass response variable (y axes) (log-transformed)
   int<lower=1> rem; //number of remote  data points
   real b2[rem]; //Biomass response variable (y axes) (log-transformed)
   int<lower=1> fis; //number of fished  data points
   real b3[fis]; //Biomass response variable (y axes) (log-transformed)
   //explanatory variables for each component
   //reserves
   int<lower=0> ag[res]; //age reserve (only for reserves)
   real hc[res]; //predictor hard coral (only for reserve/remote model components due to data)
   real si[res]; //predictor reserve size (only for reserves)
   real d[res]; //predictor depth 
   real sa[res]; //predictor sampling area
   real op[res]; //predictor ocean productivity 
   real sst[res]; //predictor sea surface temperature
   real gr[res]; //predictor gravity(not standardized)
   int at[res];//predictor atoll
   int rh_b[res];//predictor backreef/lagoon habitat
   int rh_f[res];//predictor flat habitat
   int rh_c[res];//predictor crest habitat
   int cm_pc[res];//predictor census method point count
   //remote
   real hc2[rem]; //predictor hard coral 
   real d2[rem]; //predictor depth 
   real sa2[rem]; //predictor sampling area
   real op2[rem]; //predictor ocean productivity 
   real sst2[rem]; //predictor sea surface temperature
   int at2[rem];//predictor atoll
   int rh_b2[rem];//predictor backreef/lagoon habitat
   int rh_c2[rem];//predictor crest habitat
   int cm_ds2[rem];//predictor census method distance sampling
   //fished
   real d3[fis]; //predictor depth 
   real sst3[fis]; //predictor sea surface temperature
   real sa3[fis]; //predictor sampling area
   real op3[fis]; //predictor ocean productivity 
   int at3[fis];//predictor atoll
   int rh_b3[fis];//predictor backreef/lagoon habitat
   int rh_f3[fis];//predictor flat habitat
   int rh_c3[fis];//predictor crest habitat
   int cm_ds3[fis];//predictor census method distance sampling
   int cm_pc3[fis];//predictor census method point count
   int m_r3[fis]; //predictor management restricted (for status models only)
   real gr3[fis]; //predictor gravity(not standardized)
  real hc3[fis]; //predictor hard coral (not used in model because majority fished sites dont have that info but did use in generated quantities (av conditions)))
   int<lower=1> R3; //number of data regions fished(groups)
   int<lower=1, upper=R3> pr3[fis]; //region id 
 }
 parameters {
   vector[13] beta; //slopes for sampling/environmental variables 
   real<lower=0> sigma_e; //error sd for biomass reserves
   real<lower=0> sigma_r; //error sd for biomass remote
   real<lower=0> sigma_f; //error sd for biomass fished
   vector[R3] u3; // random effects for fished regions
   real<lower=0> sigma_u; //region sd
   real log_r;// community biomass intrinsic growth rate
   real log_bmin; //biomass at reserve age 0
   real log_B0; //unfished biomass 
   real I_fished; //intercept for fished reefs
 }
 transformed parameters {
   vector[res] mu;//mean log-biomass  reserves
   vector[rem] mu2;//mean log-biomass  remote
   vector[fis] mu3;//mean log-biomass  fished
   vector[res] k; //site-specific reserve carrying capacity given enevironmental factors
   vector[rem] k2; //site-specific remote carrying capacity given enevironmental factors
   real<lower=0> r;
   real<lower=0> bmin; 
   real<lower=0> B0;
   r=exp(log_r); 
   bmin=exp(log_bmin); 
   B0=exp(log_B0);  
   //reserve component: without exports
   for (i in 1:res){ 
     k[i] = exp(log(B0) +beta[1] * op[i]+beta[2] * sst[i] +beta[3]*hc[i]+beta[4]*at[i]); 
         mu[i] = log(exp(log(k[i]*exp(log(bmin/k[i])*exp(-r*ag[i])))+beta[12]*d[i]+beta[5]*rh_c[i]+beta[6]*rh_f[i]+beta[7]*rh_b[i]+beta[8]*cm_pc[i]+beta[9]*sa[i]+ beta[10]*si[i]+ beta[13]*gr[i]));
  }
   //remote component
   for (i in 1:rem){ 
     k2[i]=exp(log(B0)+beta[1] * op2[i] +beta[2] * sst2[i]+beta[3]*hc2[i]+beta[4]*at2[i]);
     mu2[i] = log(k2[i])+beta[12]*d2[i]+beta[5]*rh_c2[i]+beta[7]*rh_b2[i]+beta[11]*cm_ds2[i]+beta[9]*sa2[i];
   }
   //fished component
   for (i in 1:fis){ 
     mu3[i] =I_fished+beta[12]*d3[i]+beta[5]*rh_c3[i]+beta[6]*rh_f3[i]+beta[7]*rh_b3[i]+beta[8]*cm_pc3[i]+beta[9]*sa3[i]+beta[11]*cm_ds3[i]+ beta[13]*gr3[i]+ u3[pr3[i]];
   }
 }
 model {
   //priors
   beta[1] ~ normal (0,2); //prior slope
   beta[2] ~ normal (0,2); //prior slope
   beta[3] ~ normal (0,2); //prior slope
   beta[4] ~ normal (0,2); //prior slope
   beta[5] ~ normal (0,2); //prior slope
   beta[6] ~ normal (0,2); //prior slope
   beta[7] ~ normal (0,2); //prior slope
   beta[8] ~ normal (0,2); //prior slope
   beta[9] ~ normal (0,2); //prior slope
   beta[10] ~ normal (0,2); //prior slope
   beta[11] ~ normal (0,2); //prior slope
    beta[12] ~ normal (0,2); //prior slope
  beta[13] ~ normal (0,2); //prior slope
      sigma_u ~ cauchy (0,1); //prior sd for group varying intercept
   u3 ~ normal(0,sigma_u); //prior re fished
   log_r ~ normal (-2,1); //weekly informative prior  biomass intrinsic growth rate
   sigma_e ~ cauchy(0,1); //uninformative prior sd
   sigma_r ~ cauchy(0,1); //uninformative prior sd
   sigma_f ~ cauchy(0,1); //uninformative prior sd
   log_bmin ~ normal (log(10),1); //weakly informative prior reserve biomass at age 0
   log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
   I_fished~ normal(5,5); // prior for intercept in fished reefs
    //likelihoods  
   for(n in 1:res){
     target += normal_lpdf(b[n] | mu[n], sigma_e);
   }
   for(n in 1:rem){
     target += normal_lpdf(b2[n] | mu2[n], sigma_r);
   }
   for(n in 1:fis){
     target += normal_lpdf(b3[n] | mu3[n], sigma_f);
   }
 }
generated quantities {
   real <lower=0> BMMSY; //BMMSY for averge and most common environemental conditions
   real <lower=0> MMSY;//MMSY for average and most common environmental conditions
     //average and common conditions
   BMMSY=B0/2.718281828;
   MMSY=((r*B0)/2.718281828); 
}
"
 writeLines(stancode_full,"~/full.stan") 
 
 
 Fit_full<- stan("~/full.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
 windows()
 pairs( Fit_full,pars=c("log_bmin","log_r","log_B0"))
 
 #posterior contraction for full model reference points
  (var_prior_MMSY-(sd(rstan::extract(Fit_full,pars=c("MMSY"))$MMSY)^2))/var_prior_MMSY
 (var_prior_BMMSY-(sd(rstan::extract(Fit_full,pars=c("BMMSY"))$BMMSY)^2))/var_prior_BMMSY
 #posterior contraction for our data
 contraction_full<-as.data.frame(cbind(round(c((var_prior_logB0-(sd(rstan::extract(Fit_full,pars=c("log_B0"))$log_B0)^2))/var_prior_logB0,(var_prior_logr-(sd(rstan::extract(Fit_full,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_logbmin-(sd(rstan::extract(Fit_full,pars=c("log_bmin"))$log_bmin)^2))/var_prior_logbmin),2),c("log_B0","log_r","log_bmin")))
 colnames( contraction_full)<-c("posterior contraction","parameter")
 print(contraction_full)
 
 
 #now we follow the bayesian workflow with simulated data
 ### BAYESIAN WORKFLOW ###  
 #1. prior predictive checks.....................................................
 nsims<-50
 #sample from priors
 B0_sim4<-exp(rnorm(nsims,log(120),1))
 bmin_sim4<-exp(rnorm(nsims,log(10),1))
 r_sim4<-exp(rnorm(nsims,-2,1))
 beta13_sim4<-rnorm(nsims,0,2)
 MMSY_sim4<-(B0_sim4*r_sim4)/2.718281828
 BMMSY_sim4<-(B0_sim4)/2.718281828
 beta1_sim4<-rnorm(nsims,0,2)
 beta2_sim4<-rnorm(nsims,0,2)
 beta3_sim4<-rnorm(nsims,0,2)
 beta4_sim4<-rnorm(nsims,0,2)
 beta5_sim4<-rnorm(nsims,0,2)
 beta6_sim4<-rnorm(nsims,0,2)
 beta7_sim4<-rnorm(nsims,0,2)
 beta8_sim4<-rnorm(nsims,0,2)
 beta9_sim4<-rnorm(nsims,0,2)
 beta10_sim4<-rnorm(nsims,0,2)
 beta11_sim4<-rnorm(nsims,0,2)
 beta12_sim4<-rnorm(nsims,0,2)
 sigma_u_sim4<-rhcauchy(nsims,1)
 u3_sim4<-matrix(NA,nrow=nlevels(as.factor(fished_data$Larger)),ncol=nsims)
 for (i in 1: length( sigma_u_sim4)){
    u3_sim4[,i]<-  rnorm(nlevels(as.factor(fished_data$Larger)),0,sigma_u_sim4[i])
 }
 sigma_e_sim4<-rhcauchy(nsims,1)
 sigma_r_sim4<-rhcauchy(nsims,1)
 sigma_f_sim4<-rhcauchy(nsims,1)
 I_fished_sim4<-rnorm(nsims,5,5) 
 #generate simulated data biomass
 k_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 mu_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 mu_obs_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 k_rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 mu_rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 mu_fis_sim4<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)
 res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 fis_sim4<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)

 for (n in 1:nsims){
    for (i in 1:nrow(reserves_complete)){
       k_res_sim4[i,n]<-exp(log(B0_sim4[n]) +beta1_sim4[n] * reserves_complete$sOcean_prod[i]+beta2_sim4[n]* reserves_complete$sSST[i] +beta3_sim4[n]*reserves_complete$sHardCoral[i]+beta4_sim4[n]*reserves_complete$Atoll[i])
       mu_res_sim4[i,n]<-  log(k_res_sim4[i,n]*exp(log(bmin_sim4[n]/k_res_sim4[i,n])*exp(-r_sim4[n]*reserves_complete$Closure.age[i])))
       mu_obs_res_sim4[i,n]<-mu_res_sim4[i,n]+beta12_sim4[n]*reserves_complete$sDepth[i]+beta5_sim4[n]*reserves_complete$crest[i]+beta6_sim4[n]*reserves_complete$flat[i]+beta7_sim4[n]*reserves_complete$backreef[i]+beta8_sim4[n]*reserves_complete$pointintercept[i]+beta9_sim4[n]*reserves_complete$sSampArea[i]+ beta10_sim4[n]*reserves_complete$sClosure.size[i]+beta13_sim4[n]*reserves_complete$sgravtot[i]
       res_sim4[i,n]<-rnorm(1,mu_obs_res_sim4[i,n],sigma_e_sim4[n])
    }
    for (i in 1:nrow(remote_complete)){
       k_rem_sim4[i,n]<-exp(log(B0_sim4[n])+beta1_sim4[n] *remote_complete$sOcean_prod[i] +beta2_sim4[n] * remote_complete$sSST[i]+beta3_sim4[n] *remote_complete$sHardCoral[i]+beta4_sim4[n] *remote_complete$Atoll[i]);
       mu_rem_sim4[i,n]<- log(k_rem_sim4[i,n])+beta12_sim4[n]*remote_complete$sDepth[i]+beta5_sim4[n]*remote_complete$crest[i]+beta7_sim4[n]*remote_complete$backreef[i]+beta11_sim4[n]*remote_complete$distancesampling[i]+beta9_sim4[n]*remote_complete$sSampArea[i]
       rem_sim4[i,n]<-rnorm(1,mu_rem_sim4[i,n],sigma_r_sim4[n])
    }   
    for (i in 1:nrow(fished_data)){
       mu_fis_sim4[i,n]<- I_fished_sim4[n]+beta12_sim4[n]*fished_data$sDepth[i]+beta5_sim4[n]*fished_data$crest[i]+beta6_sim4[n]*fished_data$flat[i]+beta7_sim4[n]*fished_data$backreef[i]+beta8_sim4[n]*fished_data$pointintercept[i]+beta11_sim4[n]*fished_data$distancesampling[i]+beta9_sim4[n]*fished_data$sSampArea[i]+beta13_sim4[n]*fished_data$sgravtot[i]+u3_sim4[fished_data$indexj[i],n]
       fis_sim4[i,n]<-rnorm(1,mu_fis_sim4[i,n],sigma_f_sim4[n])
    }
 }
 
 res_sim5<-as.data.frame(res_sim4)
 res_sim5$Closure.age<-reserves_complete$Closure.age
 res_sim6 <- melt(res_sim5 ,  id.vars = 'Closure.age')
 
 windows()
 a<-ggplot()+geom_density(aes(x=exp(melt(res_sim6)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in reserves (t/km2)")
 b<-ggplot() + geom_point(data=res_sim6,aes(Closure.age,exp(value),group = variable),alpha=0.1,col="grey")+ geom_smooth(data=res_sim6,method="gam",formula=y~s(x,k=3),aes(Closure.age,exp(value),group = variable),alpha=0.1,lty=2,lwd=0.1,col="grey")+
    geom_smooth(aes(x=reserves_complete$Closure.age,y=exp(matrixStats::rowMeans2(res_sim4))),col="black",lwd=2,method="gam",formula=y~s(x,k=3))+ylim(c(0,700))+ggtitle("")+guides(col=F)+ylab("Simulated biomass in reserves (t/km2)")+theme_classic()
 c<- ggplot()+geom_density(aes(x=exp(melt(rem_sim4)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+ggtitle("")+xlab("Simulated biomass in remote (t/km2)")
 d<- ggplot()+geom_density(aes(x=exp(melt(fis_sim4)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in fished (t/km2)")
 prior_pred_fig<-ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d")) 
 annotate_figure(prior_pred_fig,top="Prior predictive check")
 
 #2.Computational faithfullness: simulation-based calibration....................
 sbc_B0<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_bmin<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_MMSY<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_BMMSY<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta13<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_I_fished<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_e<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_f<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta1<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta2<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta3<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta4<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta5<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta6<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta7<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta8<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta9<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta10<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta11<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_beta12<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_u<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_u3<-array(NA,dim=c(4000,nsims,nlevels(as.factor(fished_data$Larger))))
 
 
 for (n in 1:nsims){
    print(paste("Starting simulation",n,sep = " "))
    
    stanDat_sbc<- list(b = res_sim4[,n], ag=reserves_complete$Closure.age,
                           res=nrow(reserves_complete), d=reserves_complete$sDepth, hc=reserves_complete$sHardCoral, si=reserves_complete$sClosure.size,
                           at=reserves_complete$Atoll,sst=reserves_complete$sSST,
                           op=reserves_complete$sOcean_prod,rh_c=reserves_complete$crest,rh_b=reserves_complete$backreef,rh_f=reserves_complete$flat,
                           cm_pc=reserves_complete$pointintercept, sa=reserves_complete$sSampArea,
                           gr=reserves_complete$sgravtot,
                           b2 = rem_sim4[,n],rem=nrow(remote_complete), d2=remote_complete$sDepth, hc2=remote_complete$sHardCoral, 
                           at2=remote_complete$Atoll,sst2=remote_complete$sSST,
                           op2=remote_complete$sOcean_prod,rh_c2=remote_complete$crest,rh_b2=remote_complete$backreef,
                           sa2=remote_complete$sSampArea,
                           cm_ds2=remote_complete$distancesampling,
                           b3 = fis_sim4[,n],fis=nrow(fished_data), d3=fished_data$sDepth, 
                           at3=fished_data$Atoll,sst3=fished_data$sSST,
                           op3=fished_data$sOcean_prod,rh_c3=fished_data$crest,rh_b3=fished_data$backreef,rh_f3=fished_data$flat,
                           cm_pc3=fished_data$pointintercept, sa3=fished_data$sSampArea,cm_ds3=fished_data$distancesampling,hc3=fished_data$sHardCoral2,
                           m_r3=fished_data$status_restricted,gr3=fished_data$sgravtot,
                           R2=nlevels(as.factor(remote_complete$Larger)),
                           pr2=remote_complete$indexj,
                           R3=nlevels(as.factor(fished_data$Larger)),
                           pr3=fished_data$indexj)
    fit_sbc_full<-stan("~/full.stan", data = stanDat_sbc,chains = 4,control = list(adapt_delta = 0.999) )
    sbc_B0[,n]<-rstan::extract(fit_sbc_full,pars=c("B0"))$B0
    sbc_bmin[,n]<-rstan::extract(fit_sbc_full,pars=c("bmin"))$bmin
    sbc_r[,n]<-rstan::extract(fit_sbc_full,pars=c("r"))$r
    sbc_beta13[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,13]
    sbc_MMSY[,n]<-rstan::extract(fit_sbc_full,pars=c("MMSY"))$MMSY
    sbc_BMMSY[,n]<-rstan::extract(fit_sbc_full,pars=c("BMMSY"))$BMMSY
    sbc_I_fished[,n]<-rstan::extract(fit_sbc_full,pars=c("I_fished"))$I_fished
    sbc_sigma_e[,n]<-rstan::extract(fit_sbc_full,pars=c("sigma_e"))$sigma_e
    sbc_sigma_r[,n]<-rstan::extract(fit_sbc_full,pars=c("sigma_r"))$sigma_r
    sbc_sigma_f[,n]<-rstan::extract(fit_sbc_full,pars=c("sigma_f"))$sigma_f
    sbc_beta1[,n]<- rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,1]
    sbc_beta2[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,2]
    sbc_beta3[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,3]
    sbc_beta4[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,4]
    sbc_beta5[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,5]
    sbc_beta6[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,6]
    sbc_beta7[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,7]
    sbc_beta8[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,8]
    sbc_beta9[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,9]
    sbc_beta10[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,10]
    sbc_beta11[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,11]
    sbc_beta12[,n]<-rstan::extract(fit_sbc_full,pars=c("beta"))$beta[,12]
    sbc_sigma_u[,n]<-rstan::extract(fit_sbc_full,pars=c("sigma_u"))$sigma_u
    sbc_u3[,n,]<-rstan::extract(fit_sbc_full,pars=c("u3"))$u3
 }
 
 windows()
 a<-ggplot()+geom_histogram(aes(x=melt(log(sbc_B0))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(120),1)),fill="red",alpha=0.5)+xlab("log (Unfished biomass, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+
    geom_text(aes(6.7,30000),label="Simulated-data
posterior",col="darkgrey")+
    geom_text(aes(6.5,18000),label="Prior",col="darkred")
 b<-ggplot()+geom_histogram(aes(x=melt(log(sbc_r))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,-2,1)),fill="red",alpha=0.5)+xlab("log (intrinsic growth rate, 1/t)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 c<-ggplot()+geom_histogram(aes(x=melt(log(sbc_bmin))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(10),1)),fill="red",alpha=0.5)+xlab("log (biomass reserve age 0, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 d<-ggplot()+geom_histogram(aes(x=melt(sbc_beta13)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size gravityd")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 e<-ggplot()+geom_histogram(aes(x=melt(sbc_I_fished)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,5,5)),fill="red",alpha=0.5)+xlab("log(intercept biomass in fished reefs, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 f<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_e)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd reserve biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 g<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_r)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd remote biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 h<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_f)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd fished biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 i<-ggplot()+geom_histogram(aes(x=melt(sbc_beta1)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size ocean productivity")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 j<-ggplot()+geom_histogram(aes(x=melt(sbc_beta2)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size sst")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 k<-ggplot()+geom_histogram(aes(x=melt(sbc_beta3)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size coral cover")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 l<-ggplot()+geom_histogram(aes(x=melt(sbc_beta4)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size atoll")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 m<-ggplot()+geom_histogram(aes(x=melt(sbc_beta5)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size crest")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 n<-ggplot()+geom_histogram(aes(x=melt(sbc_beta6)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size flat")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 o<-ggplot()+geom_histogram(aes(x=melt(sbc_beta7)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size backreef")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 p<-ggplot()+geom_histogram(aes(x=melt(sbc_beta8)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size point count meth.")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 q<-ggplot()+geom_histogram(aes(x=melt(sbc_beta9)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size sampling area")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 r<-ggplot()+geom_histogram(aes(x=melt(sbc_beta10)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size reserve size")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 s<-ggplot()+geom_histogram(aes(x=melt(sbc_beta11)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size distance sampling meth.")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 t<-ggplot()+geom_histogram(aes(x=melt(sbc_beta12)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size depth")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 u<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_u)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd random effects)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 w<-ggplot()+geom_histogram(aes(x=melt(log(sbc_BMMSY))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=log(exp(rnorm(4000*50,log(120),1))/2.718281828)),fill="red",alpha=0.5)+xlab("log (BMMSY, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 x<-ggplot()+geom_histogram(aes(x=melt(log(sbc_MMSY))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=log((exp(rnorm(4000*50,log(120),1))*exp(rnorm(4000*50,-2,1)))/2.718281828)),fill="red",alpha=0.5)+xlab("log (MMSY, t/km2/y)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 
ov_prioru<-matrix(NA,nrow=4000*50,ncol=nlevels(as.factor(fished_data$Larger)))
ov_sigmau<-rhcauchy(4000*50,1)
for (i in 1: length(ov_sigmau)){
    ov_prioru[i,]<-  rnorm(nlevels(as.factor(fished_data$Larger)),0,ov_sigmau[i])
 }
length(melt(sbc_u3)$value)
length(melt(ov_prioru)$value)
y<-ggplot()+geom_histogram(aes(x=melt(sbc_u3)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=melt(ov_prioru)$value),fill="red",alpha=0.5)+xlab("Random effects fished")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+xlim(c(-50,50))
 windows()
 ggarrange(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,y,w,x,nrow=4,ncol=6)


 
#rank plot for MMSY and BMMSY 
 sbc_rank_MMSY<- NA
 sbc_rank_BMMSY<- NA
 for (i in 1:nsims) { # Compute SBC rank
    post_MMSY_i <- sbc_MMSY[,i]
    post_BMMSY_i <- sbc_BMMSY[,i]
    idx_so <- seq(1,length(sbc_MMSY[,i]),5) # thin to remove autocorrelation (if wanted)
    sbc_rank_MMSY[i] <- sum( MMSY_sim4[i] < post_MMSY_i [idx_so] )
     sbc_rank_BMMSY[i] <- sum( BMMSY_sim4[i] < post_BMMSY_i [idx_so] )
   }
 est_parsMMSY <- data.frame(sbc_rank_MMSY)
 est_parsBMMSY <- data.frame(sbc_rank_BMMSY)
 ggarrange(ggplot()+geom_histogram(aes(x=est_parsMMSY$sbc_rank_MMSY),bins=15,col="black",fill="grey")+theme_classic()+xlab("Rank for MMSY"),ggplot()+geom_histogram(aes(x=est_parsBMMSY$sbc_rank_BMMSY),bins=15,col="black",fill="grey")+theme_classic()+xlab("Rank for BMMSY"),labels=c("a","b"))
#test distribution is uniform
  mmsy_full_ks<-ks.boot(sbc_rank_MMSY,runif(length(sbc_rank_MMSY),min(sbc_rank_MMSY),max(sbc_rank_MMSY))) 
 #low p value indicates that we reject the null hypothesis that both samples were draw from the same distribution. As it is high, we do not reject it
  bmmsy_full_ks<-ks.boot(sbc_rank_BMMSY,runif(length(sbc_rank_BMMSY),min(sbc_rank_BMMSY),max(sbc_rank_BMMSY))) 
 


 #3.Model sensitivity: z-score vs posterior contraction.........................
 
 #zscore: This measure estimates how close the posterior mean is to the truth relative to the posterior
 #standard deviation, in other words how close the entire posterior distribution is to the
 #truth value: posterior (posterior mean - prior)/sd posterior
 z_score<-matrix(NA,ncol=23,nrow=nsims)
 z_score[,1]<-(matrixStats::colMeans2(sbc_B0)-B0_sim4)/matrixStats::colSds(sbc_B0)
 z_score[,2]<-(matrixStats::colMeans2(sbc_r)-r_sim4)/matrixStats::colSds(sbc_r)
 z_score[,3]<-(matrixStats::colMeans2(sbc_bmin)-bmin_sim4)/matrixStats::colSds(sbc_bmin)
 z_score[,4]<-(matrixStats::colMeans2(sbc_beta13)-beta13_sim4)/matrixStats::colSds(sbc_beta13)
 z_score[,5]<-(matrixStats::colMeans2(sbc_I_fished)-I_fished_sim4)/matrixStats::colSds(sbc_I_fished)
 z_score[,6]<-(matrixStats::colMeans2(sbc_sigma_e)-sigma_e_sim4)/matrixStats::colSds(sbc_sigma_e)
 z_score[,7]<-(matrixStats::colMeans2(sbc_sigma_r)-sigma_r_sim4)/matrixStats::colSds(sbc_sigma_r)
 z_score[,8]<-(matrixStats::colMeans2(sbc_sigma_f)-sigma_f_sim4)/matrixStats::colSds(sbc_sigma_f)
 z_score[,9]<-(matrixStats::colMeans2(sbc_beta1)-beta1_sim4)/matrixStats::colSds(sbc_beta1)
 z_score[,10]<-(matrixStats::colMeans2(sbc_beta2)-beta2_sim4)/matrixStats::colSds(sbc_beta2)
 z_score[,11]<-(matrixStats::colMeans2(sbc_beta3)-beta3_sim4)/matrixStats::colSds(sbc_beta3)
 z_score[,12]<-(matrixStats::colMeans2(sbc_beta4)-beta4_sim4)/matrixStats::colSds(sbc_beta4)
 z_score[,13]<-(matrixStats::colMeans2(sbc_beta5)-beta5_sim4)/matrixStats::colSds(sbc_beta5)
 z_score[,14]<-(matrixStats::colMeans2(sbc_beta6)-beta6_sim4)/matrixStats::colSds(sbc_beta6)
 z_score[,15]<-(matrixStats::colMeans2(sbc_beta7)-beta7_sim4)/matrixStats::colSds(sbc_beta7)
 z_score[,16]<-(matrixStats::colMeans2(sbc_beta8)-beta8_sim4)/matrixStats::colSds(sbc_beta8)
 z_score[,17]<-(matrixStats::colMeans2(sbc_beta9)-beta9_sim4)/matrixStats::colSds(sbc_beta9)
 z_score[,18]<-(matrixStats::colMeans2(sbc_beta10)-beta10_sim4)/matrixStats::colSds(sbc_beta10)
 z_score[,19]<-(matrixStats::colMeans2(sbc_beta11)-beta11_sim4)/matrixStats::colSds(sbc_beta11)
 z_score[,20]<-(matrixStats::colMeans2(sbc_beta12)-beta12_sim4)/matrixStats::colSds(sbc_beta12)
 z_score[,21]<-(matrixStats::colMeans2(sbc_sigma_u)-sigma_u_sim4)/matrixStats::colSds(sbc_sigma_u)
 z_score[,22]<-(matrixStats::colMeans2(sbc_MMSY)-MMSY_sim4)/matrixStats::colSds(sbc_MMSY)
 z_score[,23]<-(matrixStats::colMeans2(sbc_BMMSY)-BMMSY_sim4)/matrixStats::colSds(sbc_BMMSY)
 
 colnames(z_score)<-c("B0","r","bmin","beta13","I_fished","sigma_e","sigma_r","sigma_f","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","beta9","beta10","beta11","beta12","sigma_u","MMSY","BMMSY")
 #posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
 pc<-matrix(NA,ncol=23,nrow=nsims)
  var_prior_beta<-(2^2)
 var_prior_sigma_u<-(sd(rhcauchy(4000,1))^2)

 pc[,1]<-(var_prior_logB0-(matrixStats::colSds(log(sbc_B0))^2))/var_prior_logB0
 pc[,2]<-(var_prior_logr-(matrixStats::colSds(log(sbc_r))^2))/var_prior_logr
 pc[,3]<-(var_prior_logbmin-(matrixStats::colSds(log(sbc_bmin))^2))/var_prior_logbmin
 pc[,4]<-(var_prior_beta13-(matrixStats::colSds(sbc_beta13)^2))/var_prior_beta13
 pc[,5]<-(var_prior_I_fished-(matrixStats::colSds(sbc_I_fished)^2))/var_prior_I_fished
 pc[,6]<-(var_prior_sigma_e-(matrixStats::colSds(sbc_sigma_e)^2))/var_prior_sigma_e
 pc[,7]<-(var_prior_sigma_r-(matrixStats::colSds(sbc_sigma_r)^2))/var_prior_sigma_r
 pc[,8]<-(var_prior_sigma_f-(matrixStats::colSds(sbc_sigma_f)^2))/var_prior_sigma_f
 pc[,9]<-(var_prior_beta-(matrixStats::colSds(sbc_beta1)^2))/var_prior_beta
 pc[,10]<-(var_prior_beta-(matrixStats::colSds(sbc_beta2)^2))/var_prior_beta
 pc[,11]<-(var_prior_beta-(matrixStats::colSds(sbc_beta3)^2))/var_prior_beta
 pc[,12]<-(var_prior_beta-(matrixStats::colSds(sbc_beta4)^2))/var_prior_beta
 pc[,13]<-(var_prior_beta-(matrixStats::colSds(sbc_beta5)^2))/var_prior_beta
 pc[,14]<-(var_prior_beta-(matrixStats::colSds(sbc_beta6)^2))/var_prior_beta
 pc[,15]<-(var_prior_beta-(matrixStats::colSds(sbc_beta7)^2))/var_prior_beta
 pc[,16]<-(var_prior_beta-(matrixStats::colSds(sbc_beta8)^2))/var_prior_beta
 pc[,17]<-(var_prior_beta-(matrixStats::colSds(sbc_beta9)^2))/var_prior_beta
 pc[,18]<-(var_prior_beta-(matrixStats::colSds(sbc_beta10)^2))/var_prior_beta
 pc[,19]<-(var_prior_beta-(matrixStats::colSds(sbc_beta11)^2))/var_prior_beta
 pc[,20]<-(var_prior_beta-(matrixStats::colSds(sbc_beta12)^2))/var_prior_beta
 pc[,21]<-(var_prior_sigma_u-(matrixStats::colSds(sbc_sigma_u)^2))/var_prior_sigma_u
 pc[,22]<-(var_prior_MMSY-(matrixStats::colSds(sbc_MMSY)^2))/var_prior_MMSY
 pc[,23]<-(var_prior_BMMSY-(matrixStats::colSds(sbc_BMMSY)^2))/var_prior_BMMSY
 
 colnames(pc)<-c("B0","r","bmin","beta13","I_fished","sigma_e","sigma_r","sigma_f","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","beta9","beta10","beta11","beta12","sigma_u","MMSY","BMMSY")
 
 pc_zscore<-cbind(melt(pc),melt(z_score)$value)
 ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
    geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Posterior contraction (Nsims=50)")+ylab("Z-score (Nsims=50)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)
 ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
    geom_point(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),size=3)+
    geom_text_repel(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Posterior contraction (Nsims=50)")+ylab("Z-score (Nsims=50)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)
 
 windows()
 pc_zscore_fig2<-ggplot()+geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Median posterior contraction (Nsims=50)")+ylab("Median z-score (Nsims=50)")+theme(panel.background = element_rect(fill="white",colour="black"))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")

 pc_zscore_fig<-ggplot()+geom_point(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMeans2(z_score),x=matrixStats::colMeans2(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Mean posterior contraction (Nsims=50)")+ylab("Mean z-score (Nsims=50)")+theme(panel.background = element_rect(fill="white",colour="black"))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")
 
 zscore_fig<-ggplot()+geom_violin(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2,draw_quantiles = 0.5)+theme(axis.text.x=element_blank(),panel.background = element_rect(fill="white",colour="black"))+ylim(c(5,-5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score")+geom_jitter(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2,pch=21)
 abszscore_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(x=X2,y=abs(melt(z_score)$value),fill=X2),alpha=0.2)+theme(axis.text.x=element_text(angle=90),panel.background = element_rect(fill="white",colour="black"))+ylim(c(0,5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score
(absolute)")+geom_jitter(data=pc_zscore,aes(x=X2,y=abs(melt(z_score)$value),fill=X2),alpha=0.2,pch=21)+
    geom_hline(yintercept = c(3,4),lty=3,col="red")
  zscore_fig<-ggarrange( zscore_fig, abszscore_fig,nrow=2,ncol=1,heights=c(1,0.8))
 pc_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme(panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction")+geom_jitter(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2,pch=21)

 pc_fig_zoom<- ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme(axis.text.y=element_blank(), panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction (zoom)")+geom_jitter(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2,pch=21)+coord_cartesian(xlim=c(0, 1))
  
ggarrange(zscore_fig,pc_fig, pc_fig_zoom, nrow=1,ncol=3,labels = c("a","b",""),widths=c(1,0.8,0.5))
  ggarrange(pc_zscore_fig2,pc_zscore_fig, nrow=1,ncol=2,labels = c("a","b"))
 
 #4. Posterior predictive checks...............................................

joined_sim <- rstan::extract(fit_sbc_full)
n_sims <- length(joined_sim $lp__)
y_rep_reserves <- array(NA, c(n_sims, nrow(reserves_complete)))
y_rep_remote <- array(NA, c(n_sims, nrow(remote_complete)))
y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

for (s in 1:n_sims){
   y_rep_reserves[s,] <- rnorm(nrow(reserves_complete), joined_sim$mu[s,], joined_sim$sigma_e[s])
   y_rep_remote[s,] <- rnorm(nrow(remote_complete), joined_sim$mu2[s,], joined_sim$sigma_r[s])
   y_rep_fished[s,] <- rnorm(nrow(fished_data), joined_sim$mu3[s,], joined_sim$sigma_f[s])
}
bayesplot::color_scheme_set(scheme = "gray")
a<- bayesplot::ppc_dens_overlay(res_sim4[,50],y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<- bayesplot::ppc_dens_overlay(rem_sim4[,50],y_rep_remote[1:4000,])+ggtitle("Remote")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))+guides(col=F)
c<- bayesplot::ppc_dens_overlay(fis_sim4[,50],y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))
windows()
ggarrange(a,b,c,nrow=1,ncol=3)

a2<-ppc_stat(res_sim4[,50],y_rep_reserves[1:4000,],stat="mean")+ggtitle("Reserves")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
b2<-ppc_stat(rem_sim4[,50],y_rep_remote[1:4000,],stat="mean")+ggtitle("Remote")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
c2<-ppc_stat(fis_sim4[,50],y_rep_fished[1:4000,],stat="mean")+ggtitle("Fished")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,a2,b2,c2,nrow=2,ncol=3)


#######################################################################################################################################
#save.image(file='ZamborainMasonetal_ReefSustainability_workflow_gravGF.RData')