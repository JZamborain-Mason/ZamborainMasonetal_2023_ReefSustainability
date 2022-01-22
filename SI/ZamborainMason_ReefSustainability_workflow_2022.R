# Title: Reef sustainability_SI2
# Author: Jessica Zamborain Mason
# Description: This code implements the Bayesian workflow for "Sustainable reference points for multispecies coral reef fisheries"
# R version 3.5.3 (2019-03-11)

#the first part of the code is the same as the original analyses (setting up data for model). However, then I implement the Bayesian workflow analyses to show:
 #(i)  Models produce non-biased results
 #(ii) Model used produces identifiable parameters

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
 pMiss <- function(x){sum(is.na(x))/length(x)*100}

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
         sgravtot=standardise(log(Grav_tot+1)),
         stt_market=standardise(log(TT_market_h)))

#create dummy varibles for the different categorical variables
 reefscale_data$pointintercept<- ifelse(reefscale_data$SampMethod=="Point intercept", 1,0)
 reefscale_data$backreef<- ifelse(reefscale_data$ReefHabitat=="Lagoon_Back reef", 1,0)
 reefscale_data$crest<- ifelse(reefscale_data$ReefHabitat=="Crest", 1,0)
 reefscale_data$flat<- ifelse(reefscale_data$ReefHabitat=="Flat", 1,0)
 reefscale_data$distancesampling<- ifelse(reefscale_data$SampMethod=="Distance sampling",1,0)
 reefscale_data$model_component<- ifelse(reefscale_data$section=="benchmark"&reefscale_data$reserves==1,1,ifelse(reefscale_data$section=="benchmark"&reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,2,3))

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
 fished_data$indexj<- as.numeric(fished_data$Larger)
 remote_complete$indexj<- as.numeric(remote_complete$Larger)
 reserves_complete$indexj<- as.numeric(reserves_complete$Larger)

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
                            m_r3=fished_data$status_restricted,
                            R=nlevels(reserves_complete$Larger),
                            pr=reserves_complete$indexj,
                            R2=nlevels(remote_complete$Larger),
                            pr2=remote_complete$indexj,
                            R3=nlevels(fished_data$Larger),
                            pr3=fished_data$indexj,
                            newage=newdata$Closure.age,
                            Bio=Bio, B=B)
 
 
 #null model without generated quantities (to run faster)
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
 real log_r;//intrinsic growth rate
 real log_bmin; //biomass at reserve age 0
 real log_B0; //unfished biomass 
 real I_fished; //intercept for fished reefs
 real<lower=0,upper=1> p;// proportion of biomass exported
}
transformed parameters {
 vector[res] mu;//mean log-biomass  reserves
 vector[rem] mu2;//mean log-biomass  remote
 vector[fis] mu3;//mean log-biomass  fished
 real<lower=0> r;
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass 
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
 log_bmin ~ normal (log(40),1); //weakly informative prior reserve biomass at age 0
 log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
 I_fished~ normal(5,5); // prior for intercept in fished reefs
 p~ uniform (0,1); // uninformative prior for export proportion
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

"
 writeLines(stancode_null,"~/null_pb.stan")
 Fit_null_ss<- stan("~/null_pb.stan", data = stanDat_full_country, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 windows()
 pairs(Fit_null_ss,pars=c("B0","log_r","p","bmin"))
 #posterior contraction for our model parameters ref points
 #posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
 var_prior_B0<-(1^2)
 var_prior_bmin<-(1^2)
 var_prior_logr<-(1^2)
 var_prior_p<-(sd(runif(4000,0,1))^2)
 #posterior contraction for our data
 contraction_null<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_null_ss,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_null_ss,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_null_ss,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_null_ss,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
 colnames( contraction_null)<-c("posterior contraction","parameter")
 print(contraction_null)
 
### BAYESIAN WORKFLOW ###  
#1. prior predictive checks.....................................................
 nsims<-100
 #sample from priors
 B0_sim<-exp(rnorm(nsims,log(120),1))
 bmin_sim<-exp(rnorm(nsims,log(40),1))
 r_sim<-exp(rnorm(nsims,-2,1))
 sigma_e_sim<-rhcauchy(nsims,1)
 sigma_r_sim<-rhcauchy(nsims,1)
 sigma_f_sim<-rhcauchy(nsims,1)
 I_fished_sim<-rnorm(nsims,5,5) 
 p_sim<-runif (nsims,0,1)
#generate simulated data biomass
 res_sim<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 rem_sim<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 fis_sim<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)

 for (n in 1:nsims){
   for (i in 1:nrow(reserves_complete)){
   res_sim[i,n]<-rnorm(1,log((1-p_sim[n])*(B0_sim[n]/(1+((B0_sim[n]-bmin_sim[n])/bmin_sim[n])*exp(-(r_sim[n])*reserves_complete$Closure.age[i])))),sigma_e_sim[n])
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
 sbc_B0<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_bmin<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_p<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_I_fished<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_e<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_sigma_f<-matrix(NA,nrow=4000,ncol=nsims)
 
 for (n in 1:nsims){
    stanDat_sbc<- list(b = res_sim[,n],ag=reserves_complete$Closure.age,
                                res=nrow(reserves_complete),
                                b2 = rem_sim[,n],rem=nrow(remote_complete), 
                                b3 = fis_sim[,n],fis=nrow(fished_data))
    fit_sbc<-stan(file = "Null_schaefer_ponb_review.stan", data = stanDat_sbc, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
    sbc_B0[,n]<-rstan::extract(fit_sbc,pars=c("B0"))$B0
    sbc_bmin[,n]<-rstan::extract(fit_sbc,pars=c("bmin"))$bmin
    sbc_r[,n]<-rstan::extract(fit_sbc,pars=c("r"))$r
    sbc_p[,n]<-rstan::extract(fit_sbc,pars=c("p"))$p
    sbc_I_fished[,n]<-rstan::extract(fit_sbc,pars=c("I_fished"))$I_fished
    sbc_sigma_e[,n]<-rstan::extract(fit_sbc,pars=c("sigma_e"))$sigma_e
    sbc_sigma_r[,n]<-rstan::extract(fit_sbc,pars=c("sigma_r"))$sigma_r
    sbc_sigma_f[,n]<-rstan::extract(fit_sbc,pars=c("sigma_f"))$sigma_f
    
 }

 
a<-ggplot()+geom_histogram(aes(x=melt(log(sbc_B0))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*100,log(120),1)),fill="red",alpha=0.5)+xlab("log (Unfished biomass, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+
   geom_text(aes(6.2,58000),label="Simulated-data
posterior",col="darkgrey")+
   geom_text(aes(6.5,50000),label="Prior",col="darkred")
b<-ggplot()+geom_histogram(aes(x=melt(log(sbc_r))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*100,-2,1)),fill="red",alpha=0.5)+xlab("log (intrinsic growth rate, 1/t)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
c<-ggplot()+geom_histogram(aes(x=melt(log(sbc_bmin))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*100,log(40),1)),fill="red",alpha=0.5)+xlab("log (biomass reserve age 0, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
d<-ggplot()+geom_histogram(aes(x=melt(sbc_p)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=runif(4000*100,0,1)),fill="red",alpha=0.5)+xlab("Proportion of biomass exported")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
e<-ggplot()+geom_histogram(aes(x=melt(sbc_I_fished)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*100,5,5)),fill="red",alpha=0.5)+xlab("log(intercept biomass in fished reefs, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
f<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_e)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5)+xlab("log(sd reserve biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
g<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_r)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5)+xlab("log(sd remote biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
h<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_f)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5)+xlab("log(sd fished biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
a2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_B0)),y=log(B0_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_B0)),y=log(B0_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
unfished biomass,t/km2)")+ylab(" log(unfished biomass used to simulate data,t/km2) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
b2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_r)),y=log(r_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_r)),y=log(r_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
intrinsic growth rate,1/t)")+ylab(" log(intrinsic growth rate used to simulate data,1/t) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
c2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_bmin)),y=log(bmin_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_bmin)),y=log(bmin_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
biomass reserve age 0,t/km2)")+ylab(" log(biomass reserve age 0 used to simulate data,t/km2) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
d2<-ggplot() +geom_point(aes(x=matrixStats::colMedians(sbc_p),y=p_sim),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=matrixStats::colMedians(sbc_p),y=p_sim),method="lm",alpha=0.5,col="black")+
   xlab("Simulated-data posterior median
export proportion")+ylab("Export proportion used to simulate data ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
e2<-ggplot() +geom_point(aes(x=matrixStats::colMedians(sbc_I_fished),y=I_fished_sim),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=matrixStats::colMedians(sbc_I_fished),y=I_fished_sim),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
intercept biomass fished,t/km2)")+ylab("log(intercept used to simulate data, t/km2) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
f2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_sigma_e)),y=log(sigma_e_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_sigma_e)),y=log(sigma_e_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
sd biomass reserves)")+ylab(" log(sd biomass reserves used to simulate data) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
g2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_sigma_r)),y=log(sigma_r_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_sigma_r)),y=log(sigma_r_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
sd biomass remote)")+ylab(" log(sd biomass remote used to simulate data) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
h2<-ggplot() +geom_point(aes(x=log(matrixStats::colMedians(sbc_sigma_f)),y=log(sigma_f_sim)),alpha=0.5)+geom_abline(slope = 1,intercept = 0,lty=2)+geom_smooth(aes(x=log(matrixStats::colMedians(sbc_sigma_f)),y=log(sigma_f_sim)),method="lm",alpha=0.5,col="black")+
   xlab("log(Simulated-data posterior median
sd biomass fished)")+ylab(" log(sd biomass fished used to simulate data) ")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
windows()
sbc_fig<-ggarrange(a,b,c,d,e,f,g,h,a2,b2,c2,d2,e2,f2,g2,h2,nrow=2,ncol=8)
annotate_figure(sbc_fig, top="Simulation-based calibration")

#3.Model sensitivity: z-score vs posterior contraction.........................

#zscore: This measure estimates how close the posterior mean is to the truth relative to the posterior
#standard deviation, in other words how close the entire posterior distribution is to the
#truth value: posterior (posterior mean - prior)/sd posterior
z_score<-matrix(NA,ncol=8,nrow=nsims)
z_score[,1]<-(matrixStats::colMeans2(sbc_B0)-B0_sim)/matrixStats::colSds(sbc_B0)
z_score[,2]<-(matrixStats::colMeans2(sbc_r)-r_sim)/matrixStats::colSds(sbc_r)
z_score[,3]<-(matrixStats::colMeans2(sbc_bmin)-bmin_sim)/matrixStats::colSds(sbc_bmin)
z_score[,4]<-(matrixStats::colMeans2(sbc_p)-p_sim)/matrixStats::colSds(sbc_p)
z_score[,5]<-(matrixStats::colMeans2(sbc_I_fished)-I_fished_sim)/matrixStats::colSds(sbc_I_fished)
z_score[,6]<-(matrixStats::colMeans2(sbc_sigma_e)-sigma_e_sim)/matrixStats::colSds(sbc_sigma_e)
z_score[,7]<-(matrixStats::colMeans2(sbc_sigma_r)-sigma_r_sim)/matrixStats::colSds(sbc_sigma_r)
z_score[,8]<-(matrixStats::colMeans2(sbc_sigma_f)-sigma_f_sim)/matrixStats::colSds(sbc_sigma_f)
colnames(z_score)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f")
#posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
pc<-matrix(NA,ncol=8,nrow=nsims)
var_prior_B0<-(1^2)
var_prior_bmin<-(1^2)
var_prior_logr<-(1^2)
var_prior_I_fished<-(5^2)
var_prior_p<-(sd(runif(4000,0,1))^2)
var_prior_sigma_e<-(sd(rhcauchy(4000,1))^2)
var_prior_sigma_r<-(sd(rhcauchy(4000,1))^2)
var_prior_sigma_f<-(sd(rhcauchy(4000,1))^2)


pc[,1]<-(var_prior_B0-(matrixStats::colSds(log(sbc_B0))^2))/var_prior_B0
pc[,2]<-(var_prior_logr-(matrixStats::colSds(log(sbc_r))^2))/var_prior_logr
pc[,3]<-(var_prior_bmin-(matrixStats::colSds(log(sbc_bmin))^2))/var_prior_bmin
pc[,4]<-(var_prior_p-(matrixStats::colSds(sbc_p)^2))/var_prior_p
pc[,5]<-(var_prior_I_fished-(matrixStats::colSds(sbc_I_fished)^2))/var_prior_I_fished
pc[,6]<-(var_prior_sigma_e-(matrixStats::colSds(sbc_sigma_e)^2))/var_prior_sigma_e
pc[,7]<-(var_prior_sigma_r-(matrixStats::colSds(sbc_sigma_r)^2))/var_prior_sigma_r
pc[,8]<-(var_prior_sigma_f-(matrixStats::colSds(sbc_sigma_f)^2))/var_prior_sigma_f
colnames(pc)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f")

pc_zscore<-cbind(melt(pc),melt(z_score)$value)
ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
   geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3)+
   geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
   xlab("Posterior contraction (Nsims=100)")+ylab("Z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)

windows()
pc_zscore_fig<-ggplot()+geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3,alpha=0.5)+
   geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score2),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
   xlab("Median posterior contraction (Nsims=100)")+ylab("Median z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")
zscore_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2)+theme_classic()+ylim(c(5,-5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
   xlab("")+ylab("Posterior z-score")

pc_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme_classic()+xlim(c(0,1))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
   ylab("")+xlab("Posterior contraction")
ggarrange(zscore_fig,pc_fig)

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
a<- bayesplot::ppc_dens_overlay(res_sim[,100],y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<- bayesplot::ppc_dens_overlay(rem_sim[,100],y_rep_remote[1:4000,])+ggtitle("Remote")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))+guides(col=F)
c<- bayesplot::ppc_dens_overlay(fis_sim[,100],y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))
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
  mu_reserves_sim[,i,j]<-log((1-sbc_p[,j])*(sbc_B0[,j]/(1+((sbc_B0[,j]-sbc_bmin[,j])/sbc_bmin[,j])*exp(-(sbc_r[,j])*reserves_complete$Closure.age[i]))));
 }  
 for(i in 1:nrow(remote_complete))  {
  mu_remote_sim [,i,j]<-log(sbc_B0[,j])  
 }
 for(i in 1:nrow(fished_data))  {
  mu_fished_sim [,i,j]<-sbc_I_fished[,j]  
 }
 for (s in 1:4000){
      y_rep_reserves[s,,j] <- rnorm(nrow(reserves_complete), mu_reserves_sim[s,,j], sbc_sigma_e[s,j])
      y_rep_remote[s,,j] <- rnorm(nrow(remote_complete), mu_remote_sim[s,,j], sbc_sigma_r[s,j])
      y_rep_fished[s,,j] <- rnorm(nrow(fished_data), mu_fished_sim[s,,j], sbc_sigma_f[s,j])
 }
}
windows()
a<-bayesplot::ppc_dens_overlay(res_sim[,75],y_rep_reserves[1:500,,75])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<-bayesplot::ppc_dens_overlay(rem_sim[,75],y_rep_remote[1:500,,75])+ggtitle("remote ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
c<-bayesplot::ppc_dens_overlay(fis_sim[,75],y_rep_fished[1:500,,75])+ggtitle("fished ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,nrow=1,ncol=3)


a2<-ppc_stat(res_sim[,75],y_rep_reserves[,,75],stat="mean")+ggtitle("Reserves")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
b2<-ppc_stat(rem_sim[,75],y_rep_remote[,,75],stat="mean")+ggtitle("Remote")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
c2<-ppc_stat(fis_sim[,75],y_rep_fished[,,75],stat="mean")+ggtitle("Fished")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,a2,b2,c2,nrow=2,ncol=3)

################################################################################
#Complex model 

 stancode_complex="data {
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
   real hc3[fis]; //predictor hard coral (not used in model because majority fished sites dont have that info but did use in generated quantities (av conditions)))
  //random effects
 int<lower=1> R; //number of data regions reserves (groups)
 int<lower=1, upper=R> pr[res]; //region id 
 int<lower=1> R2; //number of data regions remote (groups)
 int<lower=1, upper=R2> pr2[rem]; //region id 
 int<lower=1> R3; //number of data regions fished(groups)
 int<lower=1, upper=R3> pr3[fis]; //region id 
 //to simulate from posterior
 real newage[res];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
 }
 parameters {
   vector[12] beta; //slopes for sampling/environmental variables 
   real<lower=0> sigma_e; //error sd for biomass reserves
   real<lower=0> sigma_r; //error sd for biomass remote
   real<lower=0> sigma_f; //error sd for biomass fished
   vector[R] u_tilde; // random effects for reserve regions 
   vector[R2] u_tilde2; // random effects for remote regions
   vector[R3] u_tilde3; // random effects for fished regions
   real<lower=0> sigma_u; //region sd
   real log_r;//intrinsic growth rate
   real log_bmin; //biomass at reserve age 0
   real log_B0; //unfished biomass 
   real I_fished; //intercept for fished reefs
   real<lower=0,upper=1> p;// proportion of biomass exported
 }
 transformed parameters {
   vector[res] mu;//mean log-biomass  reserves
   vector[rem] mu2;//mean log-biomass  remote
   vector[fis] mu3;//mean log-biomass  fished
   vector[res] k; //site-specific reserve carrying capacity given enevironmental factors
   vector[rem] k2; //site-specific remote carrying capacity given enevironmental factors
   vector[R] u; // random effects for reserve regions 
   vector[R2] u2; // random effects for remote regions 
   vector[R3] u3; // random effects for fished regions
   real<lower=0> r;
   real<lower=0> bmin; //biomass at reserve age 0
   real<lower=0> B0; //unfished biomass 
   bmin=exp(log_bmin);
   B0=exp(log_B0);
   r=exp(log_r); 
   // non centered version for random effects
   for(i in 1:R){
     u[i] = sigma_u .* u_tilde[i];//non-centered
   }
   for(i in 1:R2){
     u2[i] = sigma_u .* u_tilde2[i];//non-centered
   }
   for(i in 1:R3){
     u3[i] = sigma_u .* u_tilde3[i];//non-centered
   }
   //reserve component
   for (i in 1:res){ 
     k[i] = exp(log(B0) +beta[1] * op[i]+beta[2] * sst[i] +beta[3]*hc[i]+beta[4]*at[i]); 
     mu[i] = log((1-p)*exp(log(k[i]/(1+((k[i]-bmin)/bmin)*exp(-(r)*ag[i])))+beta[12]*d[i]+beta[5]*rh_c[i]+beta[6]*rh_f[i]+beta[7]*rh_b[i]+beta[8]*cm_pc[i]+beta[9]*sa[i]+ beta[10]*si[i]+u[pr[i]]));
   }
   //remote component
   for (i in 1:rem){ 
     k2[i]=exp(log(B0)+beta[1] * op2[i] +beta[2] * sst2[i]+beta[3]*hc2[i]+beta[4]*at2[i]);
     mu2[i] = log(k2[i])+beta[12]*d2[i]+beta[5]*rh_c2[i]+beta[7]*rh_b2[i]+beta[11]*cm_ds2[i]+beta[9]*sa2[i]+u2[pr2[i]];
   }
   //fished component
   for (i in 1:fis){ 
     mu3[i] =I_fished+beta[12]*d3[i]+beta[5]*rh_c3[i]+beta[6]*rh_f3[i]+beta[7]*rh_b3[i]+beta[8]*cm_pc3[i]+beta[9]*sa3[i]+beta[11]*cm_ds3[i]+ u3[pr3[i]];
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
   sigma_u ~ cauchy (0,1); //prior sd for group varying intercept
   u_tilde ~ normal(0,1); //prior re reserves
   u_tilde2 ~ normal(0,1); //prior re fished
   u_tilde3 ~ normal(0,1); //prior re fished
   log_r ~ normal (-2,1); //weekly informative prior  biomass intrinsic growth rate
   sigma_e ~ cauchy(0,1); //uninformative prior sd
   sigma_r ~ cauchy(0,1); //uninformative prior sd
   sigma_f ~ cauchy(0,1); //uninformative prior sd
   log_bmin ~ normal (log(40),1); //weakly informative prior reserve biomass at age 0
   log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
   I_fished~ normal(5,5); // prior for intercept in fished reefs
   p~ uniform (0,1); // uninformative prior for export proportion
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

"
 writeLines(stancode_complex,"~/complex_pb.stan")
 Fit_complex_ss_pb<- stan("~/complex_pb.stan",  data = stanDat_full_country, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 windows()
 pairs( Fit_complex_ss_pb,pars=c("p","log_r","B0","bmin"))
 
 #posterior contraction for complex model parameters to our data
 contraction_complex<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_complex_ss_pb,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_complex_ss_pb,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_complex_ss_pb,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_complex_ss_pb,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","p")))
 colnames( contraction_complex)<-c("posterior contraction","parameter")
 print(contraction_complex)
 
 #now we follow the bayesian workflow with simulated data
  ### BAYESIAN WORKFLOW ###  
 #1. prior predictive checks.....................................................
 nsims<-100
 #sample from priors
 B0_sim2<-exp(rnorm(nsims,log(120),1))
 bmin_sim2<-exp(rnorm(nsims,log(40),1))
 r_sim2<-exp(rnorm(nsims,-2,1))
 beta1_sim2<-rnorm(nsims,0,2)
 beta2_sim2<-rnorm(nsims,0,2)
 beta3_sim2<-rnorm(nsims,0,2)
 beta4_sim2<-rnorm(nsims,0,2)
 beta5_sim2<-rnorm(nsims,0,2)
 beta6_sim2<-rnorm(nsims,0,2)
 beta7_sim2<-rnorm(nsims,0,2)
 beta8_sim2<-rnorm(nsims,0,2)
 beta9_sim2<-rnorm(nsims,0,2)
 beta10_sim2<-rnorm(nsims,0,2)
 beta11_sim2<-rnorm(nsims,0,2)
 beta12_sim2<-rnorm(nsims,0,2)
 sigma_u_sim2<-rhcauchy(nsims,1)
 u_tilde_sim2<-matrix(rnorm(nsims*nlevels(reserves_complete$Larger),0,1),ncol=nlevels(reserves_complete$Larger),nrow=nsims)
 u_tilde2_sim2<-matrix(rnorm(nsims*nlevels(remote_complete$Larger),0,1),ncol=nlevels(remote_complete$Larger),nrow=nsims)
 u_tilde3_sim2<-matrix(rnorm(nsims*nlevels(fished_data$Larger),0,1),ncol=nlevels(fished_data$Larger),nrow=nsims)
 sigma_e_sim2<-rhcauchy(nsims,1)
 sigma_r_sim2<-rhcauchy(nsims,1)
 sigma_f_sim2<-rhcauchy(nsims,1)
 I_fished_sim2<-rnorm(nsims,5,5) 
 p_sim2<-runif (nsims,0,1)
 #generate simulated data biomass
 k_res_sim2<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 mu_res_sim2<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 mu_obs_res_sim2<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 k_rem_sim2<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 mu_rem_sim2<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 mu_fis_sim2<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)
 res_sim2<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims)
 rem_sim2<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims)
 fis_sim2<-matrix(NA,nrow=nrow(fished_data),ncol=nsims)
 u_sim2<-matrix(NA,nrow=nlevels(reserves_complete$Larger),ncol=nsims )
 u2_sim2<-matrix(NA,nrow=nlevels(remote_complete$Larger),ncol=nsims )
 u3_sim2<-matrix(NA,nrow=nlevels(fished_data$Larger),ncol=nsims )
 
 
 for (n in 1:nsims){
    for(r1 in 1:nlevels(reserves_complete$Larger)){
       u_sim2[r1,n] =  sigma_u_sim2[n] * u_tilde_sim2[n,r1];
    }
    for(r2 in 1:nlevels(remote_complete$Larger)){
       u2_sim2[r2,n] =  sigma_u_sim2[n] * u_tilde2_sim2[n,r2];
    }
    for(r3 in 1:nlevels(fished_data$Larger)){
       u3_sim2[r3,n] =  sigma_u_sim2[n] * u_tilde3_sim2[n,r3];
    }
    for (i in 1:nrow(reserves_complete)){
       k_res_sim2[i,n]<-exp(log(B0_sim2[n]) +beta1_sim2[n] * reserves_complete$sOcean_prod[i]+beta2_sim2[n]* reserves_complete$sSST[i] +beta3_sim2[n]*reserves_complete$sHardCoral[i]+beta4_sim2[n]*reserves_complete$Atoll[i])
       mu_res_sim2[i,n]<-  log((1-p_sim2[n])*(k_res_sim2[i,n]/(1+((k_res_sim2[i,n]-bmin_sim2[n])/bmin_sim2[n])*exp(-(r_sim2[n])*reserves_complete$Closure.age[i]))))
       mu_obs_res_sim2[i,n]<-mu_res_sim2[i,n]+beta12_sim2[n]*reserves_complete$sDepth[i]+beta5_sim2[n]*reserves_complete$crest[i]+beta6_sim2[n]*reserves_complete$flat[i]+beta7_sim2[n]*reserves_complete$backreef[i]+beta8_sim2[n]*reserves_complete$pointintercept[i]+beta9_sim2[n]*reserves_complete$sSampArea[i]+ beta10_sim2[n]*reserves_complete$sClosure.size[i]+u_sim2[reserves_complete$indexj[i],n]
       res_sim2[i,n]<-rnorm(1,mu_obs_res_sim2[i,n],sigma_e_sim2[n])
    }
    for (i in 1:nrow(remote_complete)){
       k_rem_sim2[i,n]<-exp(log(B0_sim2[n])+beta1_sim2[n] *remote_complete$sOcean_prod[i] +beta2_sim2[n] * remote_complete$sSST[i]+beta3_sim2[n] *remote_complete$sHardCoral[i]+beta4_sim2[n] *remote_complete$Atoll[i]);
       mu_rem_sim2[i,n]<- log(k_rem_sim2[i,n])+beta12_sim2[n]*remote_complete$sDepth[i]+beta5_sim2[n]*remote_complete$crest[i]+beta7_sim2[n]*remote_complete$backreef[i]+beta11_sim2[n]*remote_complete$distancesampling[i]+beta9_sim2[n]*remote_complete$sSampArea[i]+u2_sim2[remote_complete$indexj[i],n]
       rem_sim2[i,n]<-rnorm(1,mu_rem_sim2[i,n],sigma_r_sim2[n])
       
    }   
    for (i in 1:nrow(fished_data)){
       mu_fis_sim2[i,n]<- I_fished_sim2[n]+beta12_sim2[n]*fished_data$sDepth[i]+beta5_sim2[n]*fished_data$crest[i]+beta6_sim2[n]*fished_data$flat[i]+beta7_sim2[n]*fished_data$backreef[i]+beta8_sim2[n]*fished_data$pointintercept[i]+beta11_sim2[n]*fished_data$distancesampling[i]+beta9_sim2[n]*fished_data$sSampArea[i]+u3_sim2[fished_data$indexj[i],n]
       fis_sim2[i,n]<-rnorm(1,mu_fis_sim2[i,n],sigma_f_sim2[n])
    }
 }
 
 res_sim3<-as.data.frame(res_sim2)
 res_sim3$Closure.age<-reserves_complete$Closure.age
 res_sim3 <- melt(res_sim3 ,  id.vars = 'Closure.age')
 
 windows()
 a<-ggplot()+geom_density(aes(x=exp(melt(res_sim2)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in reserves (t/km2)")
 b<-ggplot() + geom_point(data=res_sim3,aes(Closure.age,exp(value),group = variable),alpha=0.1,col="grey")+ geom_smooth(data=res_sim3,method="gam",aes(Closure.age,exp(value),group = variable),alpha=0.1,lty=2,lwd=0.1,col="grey")+
    geom_smooth(aes(x=reserves_complete$Closure.age,y=exp(matrixStats::rowMeans2(res_sim2))),col="black",lwd=2,method="gam")+ylim(c(0,500))+ggtitle("")+guides(col=F)+ylab("Simulated biomass in reserves (t/km2)")+theme_classic()
 c<- ggplot()+geom_density(aes(x=exp(melt(rem_sim2)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+ggtitle("")+xlab("Simulated biomass in remote (t/km2)")
 d<- ggplot()+geom_density(aes(x=exp(melt(fis_sim2)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in fished (t/km2)")
 prior_pred_fig<-ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d")) 
 annotate_figure(prior_pred_fig,top="Prior predictive check")
 
 #2.Computational faithfulness: simulation-based calibration....................
 sbc_B0<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_bmin<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_r<-matrix(NA,nrow=4000,ncol=nsims)
 sbc_p<-matrix(NA,nrow=4000,ncol=nsims)
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
 sbc_u_tilde<-array(NA,dim=c(4000,nsims_nre,nlevels(reserves_complete$Larger)))
 sbc_u2_tilde<-array(NA,dim=c(4000,nsims_nre,nlevels(remote_complete$Larger)))
 sbc_u3_tilde<-array(NA,dim=c(4000,nsims_nre,nlevels(fished_data$Larger)))
 
 
 for (n in 1:nsims){
    print(paste("Starting simulation",n,sep = " "))
   stanDat_sbc<- list(b = res_sim2[,n],    b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                       res=nrow(reserves_complete), d=reserves_complete$sDepth, hc=reserves_complete$sHardCoral, si=reserves_complete$sClosure.size,
                       at=reserves_complete$Atoll,sst=reserves_complete$sSST,
                       op=reserves_complete$sOcean_prod,rh_c=reserves_complete$crest,rh_b=reserves_complete$backreef,rh_f=reserves_complete$flat,
                       cm_pc=reserves_complete$pointintercept, sa=reserves_complete$sSampArea,
                       gr=reserves_complete$sgravtot,
                       b2 = rem_sim2[,n],rem=nrow(remote_complete), d2=remote_complete$sDepth, hc2=remote_complete$sHardCoral, 
                       at2=remote_complete$Atoll,sst2=remote_complete$sSST,
                       op2=remote_complete$sOcean_prod,rh_c2=remote_complete$crest,rh_b2=remote_complete$backreef,
                       sa2=remote_complete$sSampArea,
                       cm_ds2=remote_complete$distancesampling,
                       b3 = fis_sim2[,n],fis=nrow(fished_data), d3=fished_data$sDepth, 
                       at3=fished_data$Atoll,sst3=fished_data$sSST,
                       op3=fished_data$sOcean_prod,rh_c3=fished_data$crest,rh_b3=fished_data$backreef,rh_f3=fished_data$flat,
                       cm_pc3=fished_data$pointintercept, sa3=fished_data$sSampArea,cm_ds3=fished_data$distancesampling,hc3=fished_data$sHardCoral2,
                       m_r3=fished_data$status_restricted,
                       R=nlevels(reserves_complete$Larger),
                       pr=reserves_complete$indexj,
                       R2=nlevels(remote_complete$Larger),
                       pr2=remote_complete$indexj,
                       R3=nlevels(fished_data$Larger),
                       pr3=fished_data$indexj)
    fit_sbc<-stan("~/complex_pb.stan", data = stanDat_sbc, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
    sbc_B0[,n]<-rstan::extract(fit_sbc,pars=c("B0"))$B0
    sbc_bmin[,n]<-rstan::extract(fit_sbc,pars=c("bmin"))$bmin
    sbc_r[,n]<-rstan::extract(fit_sbc,pars=c("r"))$r
    sbc_p[,n]<-rstan::extract(fit_sbc,pars=c("p"))$p
    sbc_I_fished[,n]<-rstan::extract(fit_sbc,pars=c("I_fished"))$I_fished
    sbc_sigma_e[,n]<-rstan::extract(fit_sbc,pars=c("sigma_e"))$sigma_e
    sbc_sigma_r[,n]<-rstan::extract(fit_sbc,pars=c("sigma_r"))$sigma_r
    sbc_sigma_f[,n]<-rstan::extract(fit_sbc,pars=c("sigma_f"))$sigma_f
    sbc_beta1[,n]<- rstan::extract(fit_sbc,pars=c("beta"))$beta[,1]
    sbc_beta2[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,2]
    sbc_beta3[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,3]
    sbc_beta4[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,4]
    sbc_beta5[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,5]
    sbc_beta6[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,6]
    sbc_beta7[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,7]
    sbc_beta8[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,8]
    sbc_beta9[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,9]
    sbc_beta10[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,10]
    sbc_beta11[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,11]
    sbc_beta12[,n]<-rstan::extract(fit_sbc,pars=c("beta"))$beta[,12]
    sbc_sigma_u[,n]<-rstan::extract(fit_sbc,pars=c("sigma_u"))$sigma_u
    sbc_u_tilde[,n,]<-rstan::extract(fit_sbc,pars=c("u_tilde"))$u_tilde
    sbc_u2_tilde[,n,]<-rstan::extract(fit_sbc,pars=c("u_tilde2"))$u_tilde2
    sbc_u3_tilde[,n,]<-rstan::extract(fit_sbc,pars=c("u_tilde3"))$u_tilde3
    
 }
 
 #check whether posteriors resemble priors
 a<-ggplot()+geom_histogram(aes(x=melt(log(sbc_B0))$value),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=rnorm(4000*100,log(120),1)),fill="red",alpha=0.5,bins=10)+xlab("log (Unfished biomass, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+
    geom_text(aes(7,130000),label="Simulated-data
posterior",col="darkgrey")+
    geom_text(aes(6.9,100000),label="Prior",col="darkred")
 b<-ggplot()+geom_histogram(aes(x=melt(log(sbc_r))$value),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=rnorm(4000*100,-2,1)),fill="red",alpha=0.5,bins=10)+xlab("log (intrinsic growth rate, 1/t)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 c<-ggplot()+geom_histogram(aes(x=melt(log(sbc_bmin))$value),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=rnorm(4000*100,log(40),1)),fill="red",alpha=0.5,bins=10)+xlab("log (biomass reserve age 0, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 d<-ggplot()+geom_histogram(aes(x=melt(sbc_p)$value),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=runif(4000*100,0,1)),fill="red",alpha=0.5,bins=10)+xlab("Proportion of biomass exported")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 e<-ggplot()+geom_histogram(aes(x=melt(sbc_I_fished)$value),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=rnorm(4000*100,5,5)),fill="red",alpha=0.5,bins=10)+xlab("log(intercept biomass in fished reefs, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 f<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_e)$value)),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5,bins=10)+xlab("log(sd reserve biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 g<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_r)$value)),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5,bins=10)+xlab("log(sd remote biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 h<-ggplot()+geom_histogram(aes(x=log(melt(sbc_sigma_f)$value)),fill="black",alpha=0.5,bins=10)+geom_histogram(aes(x=log(rhcauchy(4000*100,1))),fill="red",alpha=0.5,bins=10)+xlab("log(sd fished biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 
 ggarrange(a,b,c,d,e,f,g,h,nrow=2,ncol=4)
 
 #they look ok, but check SBC rank for export proportion (Schad et al. 2021)
 sbc_rank_p <- NA
 for (i in 1:nsims) { # Compute SBC rank
    post_p_i <- sbc_p[,i]
    idx_so <- seq(1,length(sbc_p[,i]),8) # thin to remove autocorrelation (if wanted)
    sbc_rank_p[i] <- sum( p_sim2[i] < post_p_i [idx_so] )
    sbc_rank_p[i] <- sum( p_sim2[i] < post_p_i  )
 }
 est_pars <- data.frame(sbc_rank_p)
 summary(sbc_rank_p)
 hist(sbc_rank_p)
 #check if it follows a uniform distribution
 ks.test(sbc_rank_p,runif(length(sbc_rank_p),min(sbc_rank_p),max(sbc_rank_p))) #both come from uniform distributions

 #3.Model sensitivity: z-score vs posterior contraction.........................
 
 #zscore: This measure estimates how close the posterior mean is to the truth relative to the posterior
 #standard deviation, in other words how close the entire posterior distribution is to the
 #truth value: posterior (posterior mean - prior)/sd posterior
 z_score<-matrix(NA,ncol=8,nrow=nsims)
 z_score[,1]<-(matrixStats::colMeans2(sbc_B0)-B0_sim2)/matrixStats::colSds(sbc_B0)
 z_score[,2]<-(matrixStats::colMeans2(sbc_r)-r_sim2)/matrixStats::colSds(sbc_r)
 z_score[,3]<-(matrixStats::colMeans2(sbc_bmin)-bmin_sim2)/matrixStats::colSds(sbc_bmin)
 z_score[,4]<-(matrixStats::colMeans2(sbc_p)-p_sim2)/matrixStats::colSds(sbc_p)
 z_score[,5]<-(matrixStats::colMeans2(sbc_I_fished)-I_fished_sim2)/matrixStats::colSds(sbc_I_fished)
 z_score[,6]<-(matrixStats::colMeans2(sbc_sigma_e)-sigma_e_sim2)/matrixStats::colSds(sbc_sigma_e)
 z_score[,7]<-(matrixStats::colMeans2(sbc_sigma_r)-sigma_r_sim2)/matrixStats::colSds(sbc_sigma_r)
 z_score[,8]<-(matrixStats::colMeans2(sbc_sigma_f)-sigma_f_sim2)/matrixStats::colSds(sbc_sigma_f)
 colnames(z_score)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f")
 #posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
 pc<-matrix(NA,ncol=8,nrow=nsims)
 var_prior_B0<-(1^2)
 var_prior_bmin<-(1^2)
 var_prior_logr<-(1^2)
 var_prior_I_fished<-(5^2)
 var_prior_p<-(sd(runif(4000,0,1))^2)
 var_prior_sigma_e<-(sd(rhcauchy(4000,1))^2)
 var_prior_sigma_r<-(sd(rhcauchy(4000,1))^2)
 var_prior_sigma_f<-(sd(rhcauchy(4000,1))^2)
 
 
 pc[,1]<-(var_prior_B0-(matrixStats::colSds(log(sbc_B0))^2))/var_prior_B0
 pc[,2]<-(var_prior_logr-(matrixStats::colSds(log(sbc_r))^2))/var_prior_logr
 pc[,3]<-(var_prior_bmin-(matrixStats::colSds(log(sbc_bmin))^2))/var_prior_bmin
 pc[,4]<-(var_prior_p-(matrixStats::colSds(sbc_p)^2))/var_prior_p
 pc[,5]<-(var_prior_I_fished-(matrixStats::colSds(sbc_I_fished)^2))/var_prior_I_fished
 pc[,6]<-(var_prior_sigma_e-(matrixStats::colSds(sbc_sigma_e)^2))/var_prior_sigma_e
 pc[,7]<-(var_prior_sigma_r-(matrixStats::colSds(sbc_sigma_r)^2))/var_prior_sigma_r
 pc[,8]<-(var_prior_sigma_f-(matrixStats::colSds(sbc_sigma_f)^2))/var_prior_sigma_f
 colnames(pc)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f")
 
 pc_zscore<-cbind(melt(pc),melt(z_score)$value)
 ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
    geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Posterior contraction (Nsims=100)")+ylab("Z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)
 
 windows()
 pc_zscore_fig<-ggplot()+geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Median posterior contraction (Nsims=100)")+ylab("Median z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")
 zscore_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2)+theme_classic()+ylim(c(5,-5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score")
 
 pc_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme_classic()+xlim(c(0,1))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction")
 ggarrange(zscore_fig,pc_fig, pc_zscore_fig,nrow=1,ncol=3,labels = c("a","b","c"))

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
a<- bayesplot::ppc_dens_overlay(res_sim2[,100],y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<- bayesplot::ppc_dens_overlay(rem_sim2[,100],y_rep_remote[1:4000,])+ggtitle("Remote")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))+guides(col=F)
c<- bayesplot::ppc_dens_overlay(fis_sim2[,100],y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))
windows()
ggarrange(a,b,c,nrow=1,ncol=3)

a2<-ppc_stat(res_sim2[,100],y_rep_reserves[1:4000,],stat="mean")+ggtitle("Reserves")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
b2<-ppc_stat(rem_sim2[,100],y_rep_remote[1:4000,],stat="mean")+ggtitle("Remote")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
c2<-ppc_stat(fis_sim2[,100],y_rep_fished[1:4000,],stat="mean")+ggtitle("Fished")+ labs(x = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b,c,a2,b2,c2,nrow=2,ncol=3)

 ########################################################################


#run full model (without random effects in referece point model component given space-for-time substitution)

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
   real hc3[fis]; //predictor hard coral (not used in model because majority fished sites dont have that info but did use in generated quantities (av conditions)))
 int<lower=1> R3; //number of data regions fished(groups)
 int<lower=1, upper=R3> pr3[fis]; //region id 
 }
 parameters {
   vector[12] beta; //slopes for sampling/environmental variables 
   real<lower=0> sigma_e; //error sd for biomass reserves
   real<lower=0> sigma_r; //error sd for biomass remote
   real<lower=0> sigma_f; //error sd for biomass fished
   vector[R3] u3; // random effects for fished regions
   real<lower=0> sigma_u; //region sd
   real log_r;//intrinsic growth rate
   real log_bmin; //biomass at reserve age 0
   real log_B0; //unfished biomass 
   real I_fished; //intercept for fished reefs
   real<lower=0,upper=1> p;// proportion of biomass exported
 }
 transformed parameters {
   vector[res] mu;//mean log-biomass  reserves
   vector[rem] mu2;//mean log-biomass  remote
   vector[fis] mu3;//mean log-biomass  fished
   vector[res] k; //site-specific reserve carrying capacity given enevironmental factors
   vector[rem] k2; //site-specific remote carrying capacity given enevironmental factors
   real<lower=0> r;
   real<lower=0> bmin; //biomass at reserve age 0
   real<lower=0> B0; //unfished biomass 
   r=exp(log_r); 
   bmin=exp(log_bmin);
   B0=exp(log_B0);
   //reserve component
   for (i in 1:res){ 
     k[i] = exp(log(B0) +beta[1] * op[i]+beta[2] * sst[i] +beta[3]*hc[i]+beta[4]*at[i]); 
     mu[i] = log((1-p)*exp(log(k[i]/(1+((k[i]-bmin)/bmin)*exp(-(r)*ag[i])))+beta[12]*d[i]+beta[5]*rh_c[i]+beta[6]*rh_f[i]+beta[7]*rh_b[i]+beta[8]*cm_pc[i]+beta[9]*sa[i]+ beta[10]*si[i]));
   }
   //remote component
   for (i in 1:rem){ 
     k2[i]=exp(log(B0)+beta[1] * op2[i] +beta[2] * sst2[i]+beta[3]*hc2[i]+beta[4]*at2[i]);
     mu2[i] = log(k2[i])+beta[12]*d2[i]+beta[5]*rh_c2[i]+beta[7]*rh_b2[i]+beta[11]*cm_ds2[i]+beta[9]*sa2[i];
   }
   //fished component
   for (i in 1:fis){ 
     mu3[i] =I_fished+beta[12]*d3[i]+beta[5]*rh_c3[i]+beta[6]*rh_f3[i]+beta[7]*rh_b3[i]+beta[8]*cm_pc3[i]+beta[9]*sa3[i]+beta[11]*cm_ds3[i]+ u3[pr3[i]];
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
   sigma_u ~ cauchy (0,1); //prior sd for group varying intercept
   u3 ~ normal(0,sigma_u); //prior re fished
   log_r ~ normal (-2,1); //weekly informative prior  biomass intrinsic growth rate
   sigma_e ~ cauchy(0,1); //uninformative prior sd
   sigma_r ~ cauchy(0,1); //uninformative prior sd
   sigma_f ~ cauchy(0,1); //uninformative prior sd
   bmin ~ lognormal (log(10),1); //weakly informative prior reserve biomass at age 0
   B0 ~ lognormal(log(120),1);//weakly informative prior for unfished biomass
   I_fished~ normal(5,5); // prior for intercept in fished reefs
   p~ uniform (0,1); // uninformative prior for export proportion
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
 
"
 writeLines(stancode_full,"~/full_pb.stan") 
 
 
 Fit_full_ss_pb_nreref<- stan("~/full_pb.stan", data = stanDat_full_country, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 windows()
 pairs( Fit_full_ss_pb_nreref,pars=c("p","log_r","B0"))
 
 #posterior contraction for full model parameters fitted to our data
 var_prior_B0<-(1^2)
 var_prior_bmin<-(1^2)
 var_prior_logr<-(1^2)
 var_prior_I_fished<-(5^2)
 var_prior_p<-(sd(runif(4000,0,1))^2)
 #posterior contraction for our data
 (var_prior_B0-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_B0"))$log_B0)^2))/var_prior_B0
 (var_prior_p-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("p"))$p)^2))/var_prior_p
 (var_prior_logr-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_r"))$log_r)^2))/var_prior_logr
 (var_prior_bmin-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin
 
 #posterior contraction for our data
 contraction_fullre<-as.data.frame(cbind(round(c((var_prior_B0-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_B0"))$log_B0)^2))/var_prior_B0,(var_prior_p-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("p"))$p)^2))/var_prior_p,(var_prior_logr-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_bmin-(sd(rstan::extract(Fit_full_ss_pb_nreref,pars=c("log_bmin"))$log_bmin)^2))/var_prior_bmin),2),c("B0","p","log_r","bmin")))
 colnames( contraction_fullre)<-c("posterior contraction","parameter")
 print(contraction_fullre)
 
 #rank plots 
 samp_full_open<- as.array(Fit_full_ss_pb_nre)
 mcmc_hist_r_scale(samp_full_open[, , "p"])
 
 
 #now we follow the bayesian workflow with simulated data
 ### BAYESIAN WORKFLOW ###  
 #1. prior predictive checks.....................................................
 nsims_nreref<-100
 #sample from priors
 B0_sim4<-exp(rnorm(nsims_nreref,log(120),1))
 bmin_sim4<-exp(rnorm(nsims_nreref,log(40),1))
 r_sim4<-exp(rnorm(nsims_nreref,-2,1))
 beta1_sim4<-rnorm(nsims_nreref,0,2)
 beta2_sim4<-rnorm(nsims_nreref,0,2)
 beta3_sim4<-rnorm(nsims_nreref,0,2)
 beta4_sim4<-rnorm(nsims_nreref,0,2)
 beta5_sim4<-rnorm(nsims_nreref,0,2)
 beta6_sim4<-rnorm(nsims_nreref,0,2)
 beta7_sim4<-rnorm(nsims_nreref,0,2)
 beta8_sim4<-rnorm(nsims_nreref,0,2)
 beta9_sim4<-rnorm(nsims_nreref,0,2)
 beta10_sim4<-rnorm(nsims_nreref,0,2)
 beta11_sim4<-rnorm(nsims_nreref,0,2)
 beta12_sim4<-rnorm(nsims_nreref,0,2)
 sigma_u_sim4<-rhcauchy(nsims_nreref,1)
 u3_sim4<-matrix(NA,nrow=nlevels(fished_data$Larger),ncol=nsims_nreref)
 for (i in 1: length( sigma_u_sim4)){
    u3_sim4[,i]<-  rnorm(nlevels(fished_data$Larger),0,sigma_u_sim4[i])
 
 }
 sigma_e_sim4<-rhcauchy(nsims_nreref,1)
 sigma_r_sim4<-rhcauchy(nsims_nreref,1)
 sigma_f_sim4<-rhcauchy(nsims_nreref,1)
 I_fished_sim4<-rnorm(nsims_nreref,5,5) 
 p_sim4<-runif (nsims_nreref,0,1)
 #generate simulated data biomass
 k_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims_nreref)
 mu_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims_nreref)
 mu_obs_res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims_nreref)
 k_rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims_nreref)
 mu_rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims_nreref)
 mu_fis_sim4<-matrix(NA,nrow=nrow(fished_data),ncol=nsims_nreref)
 res_sim4<-matrix(NA,nrow=nrow(reserves_complete),ncol=nsims_nreref)
 rem_sim4<-matrix(NA,nrow=nrow(remote_complete),ncol=nsims_nreref)
 fis_sim4<-matrix(NA,nrow=nrow(fished_data),ncol=nsims_nreref)

 for (n in 1:nsims_nreref){
    
    for (i in 1:nrow(reserves_complete)){
       k_res_sim4[i,n]<-exp(log(B0_sim4[n]) +beta1_sim4[n] * reserves_complete$sOcean_prod[i]+beta2_sim4[n]* reserves_complete$sSST[i] +beta3_sim4[n]*reserves_complete$sHardCoral[i]+beta4_sim4[n]*reserves_complete$Atoll[i])
       mu_res_sim4[i,n]<-  log((1-p_sim4[n])*(k_res_sim4[i,n]/(1+((k_res_sim4[i,n]-bmin_sim4[n])/bmin_sim4[n])*exp(-(r_sim4[n])*reserves_complete$Closure.age[i]))))
       mu_obs_res_sim4[i,n]<-mu_res_sim4[i,n]+beta12_sim4[n]*reserves_complete$sDepth[i]+beta5_sim4[n]*reserves_complete$crest[i]+beta6_sim4[n]*reserves_complete$flat[i]+beta7_sim4[n]*reserves_complete$backreef[i]+beta8_sim4[n]*reserves_complete$pointintercept[i]+beta9_sim4[n]*reserves_complete$sSampArea[i]+ beta10_sim4[n]*reserves_complete$sClosure.size[i]
       res_sim4[i,n]<-rnorm(1,mu_obs_res_sim4[i,n],sigma_e_sim4[n])
    }
    for (i in 1:nrow(remote_complete)){
       k_rem_sim4[i,n]<-exp(log(B0_sim4[n])+beta1_sim4[n] *remote_complete$sOcean_prod[i] +beta2_sim4[n] * remote_complete$sSST[i]+beta3_sim4[n] *remote_complete$sHardCoral[i]+beta4_sim4[n] *remote_complete$Atoll[i]);
       mu_rem_sim4[i,n]<- log(k_rem_sim4[i,n])+beta12_sim4[n]*remote_complete$sDepth[i]+beta5_sim4[n]*remote_complete$crest[i]+beta7_sim4[n]*remote_complete$backreef[i]+beta11_sim4[n]*remote_complete$distancesampling[i]+beta9_sim4[n]*remote_complete$sSampArea[i]
       rem_sim4[i,n]<-rnorm(1,mu_rem_sim4[i,n],sigma_r_sim4[n])
       
    }   
    for (i in 1:nrow(fished_data)){
       mu_fis_sim4[i,n]<- I_fished_sim4[n]+beta12_sim4[n]*fished_data$sDepth[i]+beta5_sim4[n]*fished_data$crest[i]+beta6_sim4[n]*fished_data$flat[i]+beta7_sim4[n]*fished_data$backreef[i]+beta8_sim4[n]*fished_data$pointintercept[i]+beta11_sim4[n]*fished_data$distancesampling[i]+beta9_sim4[n]*fished_data$sSampArea[i]+u3_sim4[fished_data$indexj[i],n]
       fis_sim4[i,n]<-rnorm(1,mu_fis_sim4[i,n],sigma_f_sim4[n])
    }
 }
 
 res_sim5<-as.data.frame(res_sim4)
 res_sim5$Closure.age<-reserves_complete$Closure.age
 res_sim6 <- melt(res_sim5 ,  id.vars = 'Closure.age')
 
 windows()
 a<-ggplot()+geom_density(aes(x=exp(melt(res_sim6)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in reserves (t/km2)")
 b<-ggplot() + geom_point(data=res_sim6,aes(Closure.age,exp(value),group = variable),alpha=0.1,col="grey")+ geom_smooth(data=res_sim6,method="gam",aes(Closure.age,exp(value),group = variable),alpha=0.1,lty=2,lwd=0.1,col="grey")+
    geom_smooth(aes(x=reserves_complete$Closure.age,y=exp(matrixStats::rowMeans2(res_sim4))),col="black",lwd=2,method="gam")+ylim(c(0,500))+ggtitle("")+guides(col=F)+ylab("Simulated biomass in reserves (t/km2)")+theme_classic()
 c<- ggplot()+geom_density(aes(x=exp(melt(rem_sim4)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+ggtitle("")+xlab("Simulated biomass in remote (t/km2)")
 d<- ggplot()+geom_density(aes(x=exp(melt(fis_sim4)$value)),fill="black",alpha=0.5,col="black")+xlim(c(0,500))+theme_classic()+xlab("Simulated biomass in fished (t/km2)")
 prior_pred_fig<-ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d")) 
 annotate_figure(prior_pred_fig,top="Prior predictive check")
 
 #2.Computational faithfullness: simulation-based calibration....................
 sbc_nreref_B0<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_bmin<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_r<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_p<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_I_fished<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_sigma_e<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_sigma_r<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_sigma_f<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta1<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta2<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta3<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta4<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta5<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta6<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta7<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta8<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta9<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta10<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta11<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_beta12<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_sigma_u<-matrix(NA,nrow=4000,ncol=nsims_nreref)
 sbc_nreref_u3<-array(NA,dim=c(4000,nsims_nreref,nlevels(fished_data$Larger)))
 
 
 for (n in 1:nsims_nreref){
    print(paste("Starting simulation",n,sep = " "))
    
    stanDat_sbc_nreref<- list(b = res_sim4[,n],    b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
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
                           m_r3=fished_data$status_restricted,
                           R2=nlevels(remote_complete$Larger),
                           pr2=remote_complete$indexj,
                           R3=nlevels(fished_data$Larger),
                           pr3=fished_data$indexj)
    fit_sbc_nreref<-stan("~/full_pb.stan", data = stanDat_sbc_nreref,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9) )
    sbc_nreref_B0[,n]<-rstan::extract(fit_sbc_nreref,pars=c("B0"))$B0
    sbc_nreref_bmin[,n]<-rstan::extract(fit_sbc_nreref,pars=c("bmin"))$bmin
    sbc_nreref_r[,n]<-rstan::extract(fit_sbc_nreref,pars=c("r"))$r
    sbc_nreref_p[,n]<-rstan::extract(fit_sbc_nreref,pars=c("p"))$p
    sbc_nreref_I_fished[,n]<-rstan::extract(fit_sbc_nreref,pars=c("I_fished"))$I_fished
    sbc_nreref_sigma_e[,n]<-rstan::extract(fit_sbc_nreref,pars=c("sigma_e"))$sigma_e
    sbc_nreref_sigma_r[,n]<-rstan::extract(fit_sbc_nreref,pars=c("sigma_r"))$sigma_r
    sbc_nreref_sigma_f[,n]<-rstan::extract(fit_sbc_nreref,pars=c("sigma_f"))$sigma_f
    sbc_nreref_beta1[,n]<- rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,1]
    sbc_nreref_beta2[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,2]
    sbc_nreref_beta3[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,3]
    sbc_nreref_beta4[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,4]
    sbc_nreref_beta5[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,5]
    sbc_nreref_beta6[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,6]
    sbc_nreref_beta7[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,7]
    sbc_nreref_beta8[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,8]
    sbc_nreref_beta9[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,9]
    sbc_nreref_beta10[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,10]
    sbc_nreref_beta11[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,11]
    sbc_nreref_beta12[,n]<-rstan::extract(fit_sbc_nreref,pars=c("beta"))$beta[,12]
    sbc_nreref_sigma_u[,n]<-rstan::extract(fit_sbc_nreref,pars=c("sigma_u"))$sigma_u
    sbc_nreref_u3[,n,]<-rstan::extract(fit_sbc_nreref,pars=c("u3"))$u3
    
 }
 
 
 windows()
 a<-ggplot()+geom_histogram(aes(x=melt(log(sbc_nreref_B0))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(120),1)),fill="red",alpha=0.5)+xlab("log (Unfished biomass, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+
    geom_text(aes(6.2,30000),label="Simulated-data
posterior",col="darkgrey")+
    geom_text(aes(6.5,28000),label="Prior",col="darkred")
 b<-ggplot()+geom_histogram(aes(x=melt(log(sbc_nreref_r))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,-2,1)),fill="red",alpha=0.5)+xlab("log (intrinsic growth rate, 1/t)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 c<-ggplot()+geom_histogram(aes(x=melt(log(sbc_nreref_bmin))$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,log(40),1)),fill="red",alpha=0.5)+xlab("log (biomass reserve age 0, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 d<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_p)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=runif(4000*50,0,1)),fill="red",alpha=0.5)+xlab("Proportion of biomass exported")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 e<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_I_fished)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,5,5)),fill="red",alpha=0.5)+xlab("log(intercept biomass in fished reefs, t/km2)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 f<-ggplot()+geom_histogram(aes(x=log(melt(sbc_nreref_sigma_e)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd reserve biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 g<-ggplot()+geom_histogram(aes(x=log(melt(sbc_nreref_sigma_r)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd remote biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 h<-ggplot()+geom_histogram(aes(x=log(melt(sbc_nreref_sigma_f)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd fished biomass)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 i<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta1)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size ocean productivity")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 j<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta2)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size sst")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 k<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta3)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size coral cover")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 l<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta4)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size atoll")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 m<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta5)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size crest")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 n<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta6)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size flat")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 o<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta7)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size backreef")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 p<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta8)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size point count meth.")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 q<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta9)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size sampling area")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 r<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta10)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size reserve size")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 s<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta11)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size distance sampling meth.")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 t<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_beta12)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=rnorm(4000*50,0,2)),fill="red",alpha=0.5)+xlab("effect size depth")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
 u<-ggplot()+geom_histogram(aes(x=log(melt(sbc_nreref_sigma_u)$value)),fill="black",alpha=0.5)+geom_histogram(aes(x=log(rhcauchy(4000*50,1))),fill="red",alpha=0.5)+xlab("log(sd random effects)")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))
ov_prioru<-matrix(NA,nrow=4000*50,ncol=nlevels(fished_data$Larger))
ov_sigmau<-rhcauchy(4000*50,1)
for (i in 1: length(ov_sigmau)){
    ov_prioru[i,]<-  rnorm(nlevels(fished_data$Larger),0,ov_sigmau[i])
 }
length(melt(sbc_nreref_u3)$value)
length(melt(ov_prioru)$value)
x<-ggplot()+geom_histogram(aes(x=melt(sbc_nreref_u3)$value),fill="black",alpha=0.5)+geom_histogram(aes(x=melt(ov_prioru)$value),fill="red",alpha=0.5)+xlab("Random effects fished")+theme(text = element_text(size=8),panel.background = element_rect(fill="white",colour="black"))+xlim(c(-50,50))
 
 ggarrange(a,b,c,d,e,f,g,h,nrow=2,ncol=4)
 ggarrange(i,j,k,l,m,n,o,p,q,r,s,t,nrow=2,ncol=6)
x
 
#rank plot for p 
 sbc_rank_p <- NA
 for (i in 1:nsims_nreref) { # Compute SBC rank
    post_p_i <- sbc_nreref_p[,i]
    idx_so <- seq(1,length(sbc_nreref_p[,i]),8) # thin to remove autocorrelation (if wanted)
    sbc_rank_p[i] <- sum( p_sim4[i] < post_p_i [idx_so] )
    sbc_rank_p[i] <- sum( p_sim4[i] < post_p_i  )
 }
 est_pars <- data.frame(sbc_rank_p)
 summary(sbc_rank_p)
 hist(sbc_rank_p)
 #test distribution is uniform
 ks.test(sbc_rank_p,runif(length(sbc_rank_p),min(sbc_rank_p),max(sbc_rank_p)))
 
 #3.Model sensitivity: z-score vs posterior contraction.........................
 
 #zscore: This measure estimates how close the posterior mean is to the truth relative to the posterior
 #standard deviation, in other words how close the entire posterior distribution is to the
 #truth value: posterior (posterior mean - prior)/sd posterior
 z_score<-matrix(NA,ncol=21,nrow=nsims_nreref)
 z_score[,1]<-(matrixStats::colMeans2(sbc_nreref_B0)-B0_sim4)/matrixStats::colSds(sbc_nreref_B0)
 z_score[,2]<-(matrixStats::colMeans2(sbc_nreref_r)-r_sim4)/matrixStats::colSds(sbc_nreref_r)
 z_score[,3]<-(matrixStats::colMeans2(sbc_nreref_bmin)-bmin_sim4)/matrixStats::colSds(sbc_nreref_bmin)
 z_score[,4]<-(matrixStats::colMeans2(sbc_nreref_p)-p_sim4)/matrixStats::colSds(sbc_nreref_p)
 z_score[,5]<-(matrixStats::colMeans2(sbc_nreref_I_fished)-I_fished_sim4)/matrixStats::colSds(sbc_nreref_I_fished)
 z_score[,6]<-(matrixStats::colMeans2(sbc_nreref_sigma_e)-sigma_e_sim4)/matrixStats::colSds(sbc_nreref_sigma_e)
 z_score[,7]<-(matrixStats::colMeans2(sbc_nreref_sigma_r)-sigma_r_sim4)/matrixStats::colSds(sbc_nreref_sigma_r)
 z_score[,8]<-(matrixStats::colMeans2(sbc_nreref_sigma_f)-sigma_f_sim4)/matrixStats::colSds(sbc_nreref_sigma_f)
 z_score[,9]<-(matrixStats::colMeans2(sbc_nreref_beta1)-beta1_sim4)/matrixStats::colSds(sbc_nreref_beta1)
 z_score[,10]<-(matrixStats::colMeans2(sbc_nreref_beta2)-beta2_sim4)/matrixStats::colSds(sbc_nreref_beta2)
 z_score[,11]<-(matrixStats::colMeans2(sbc_nreref_beta3)-beta3_sim4)/matrixStats::colSds(sbc_nreref_beta3)
 z_score[,12]<-(matrixStats::colMeans2(sbc_nreref_beta4)-beta4_sim4)/matrixStats::colSds(sbc_nreref_beta4)
 z_score[,13]<-(matrixStats::colMeans2(sbc_nreref_beta5)-beta5_sim4)/matrixStats::colSds(sbc_nreref_beta5)
 z_score[,14]<-(matrixStats::colMeans2(sbc_nreref_beta6)-beta6_sim4)/matrixStats::colSds(sbc_nreref_beta6)
 z_score[,15]<-(matrixStats::colMeans2(sbc_nreref_beta7)-beta7_sim4)/matrixStats::colSds(sbc_nreref_beta7)
 z_score[,16]<-(matrixStats::colMeans2(sbc_nreref_beta8)-beta8_sim4)/matrixStats::colSds(sbc_nreref_beta8)
 z_score[,17]<-(matrixStats::colMeans2(sbc_nreref_beta9)-beta9_sim4)/matrixStats::colSds(sbc_nreref_beta9)
 z_score[,18]<-(matrixStats::colMeans2(sbc_nreref_beta10)-beta10_sim4)/matrixStats::colSds(sbc_nreref_beta10)
 z_score[,19]<-(matrixStats::colMeans2(sbc_nreref_beta11)-beta11_sim4)/matrixStats::colSds(sbc_nreref_beta11)
 z_score[,20]<-(matrixStats::colMeans2(sbc_nreref_beta12)-beta12_sim4)/matrixStats::colSds(sbc_nreref_beta12)
 z_score[,21]<-(matrixStats::colMeans2(sbc_nreref_sigma_u)-sigma_u_sim4)/matrixStats::colSds(sbc_nreref_sigma_u)
 
 colnames(z_score)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","beta9","beta10","beta11","beta12","sigma_u")
 #posterior contraction: estimates how much prior uncertainty is reduced in the posterior estimation
 pc<-matrix(NA,ncol=21,nrow=nsims_nreref)
 var_prior_B0<-(1^2)
 var_prior_bmin<-(1^2)
 var_prior_logr<-(1^2)
 var_prior_I_fished<-(5^2)
 var_prior_p<-(sd(runif(4000,0,1))^2)
 var_prior_sigma_e<-(sd(rhcauchy(4000,1))^2)
 var_prior_sigma_r<-(sd(rhcauchy(4000,1))^2)
 var_prior_sigma_f<-(sd(rhcauchy(4000,1))^2)
 var_prior_beta<-(2^2)
 var_prior_sigma_u<-(sd(rhcauchy(4000,1))^2)
 
 pc[,1]<-(var_prior_B0-(matrixStats::colSds(log(sbc_nreref_B0))^2))/var_prior_B0
 pc[,2]<-(var_prior_logr-(matrixStats::colSds(log(sbc_nreref_r))^2))/var_prior_logr
 pc[,3]<-(var_prior_bmin-(matrixStats::colSds(log(sbc_nreref_bmin))^2))/var_prior_bmin
 pc[,4]<-(var_prior_p-(matrixStats::colSds(sbc_nreref_p)^2))/var_prior_p
 pc[,5]<-(var_prior_I_fished-(matrixStats::colSds(sbc_nreref_I_fished)^2))/var_prior_I_fished
 pc[,6]<-(var_prior_sigma_e-(matrixStats::colSds(sbc_nreref_sigma_e)^2))/var_prior_sigma_e
 pc[,7]<-(var_prior_sigma_r-(matrixStats::colSds(sbc_nreref_sigma_r)^2))/var_prior_sigma_r
 pc[,8]<-(var_prior_sigma_f-(matrixStats::colSds(sbc_nreref_sigma_f)^2))/var_prior_sigma_f
 pc[,9]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta1)^2))/var_prior_beta
 pc[,10]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta2)^2))/var_prior_beta
 pc[,11]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta3)^2))/var_prior_beta
 pc[,12]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta4)^2))/var_prior_beta
 pc[,13]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta5)^2))/var_prior_beta
 pc[,14]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta6)^2))/var_prior_beta
 pc[,15]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta7)^2))/var_prior_beta
 pc[,16]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta8)^2))/var_prior_beta
 pc[,17]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta9)^2))/var_prior_beta
 pc[,18]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta10)^2))/var_prior_beta
 pc[,19]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta11)^2))/var_prior_beta
 pc[,20]<-(var_prior_beta-(matrixStats::colSds(sbc_nreref_beta12)^2))/var_prior_beta
 pc[,21]<-(var_prior_sigma_u-(matrixStats::colSds(sbc_nreref_sigma_u)^2))/var_prior_sigma_u
 
 colnames(pc)<-c("B0","r","bmin","p","I_fished","sigma_e","sigma_r","sigma_f","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","beta9","beta10","beta11","beta12","sigma_u")
 
 pc_zscore<-cbind(melt(pc),melt(z_score)$value)
 ggplot()+geom_point(data=pc_zscore,aes(x=value,y=melt(z_score)$value,col=X2),alpha=0.2)+
    geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="darkgrey",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Posterior contraction (Nsims=100)")+ylab("Z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2)+geom_vline(xintercept = 0.5,lty=2)
 
 windows()
 pc_zscore_fig<-ggplot()+geom_point(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),size=3,alpha=0.5)+
    geom_text_repel(aes(y=matrixStats::colMedians(z_score),x=matrixStats::colMedians(pc)),label=colnames(z_score),col="black",family="Helvetica")+ylim(-5,5)+xlim(0,1)+
    xlab("Median posterior contraction (Nsims=100)")+ylab("Median z-score (Nsims=100)")+theme_classic()+geom_hline(yintercept = 0,lty=2,col="darkgrey")+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")
 zscore_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(x=X2,y=melt(z_score)$value,fill=X2),alpha=0.2)+theme(axis.text.x = element_text(angle = 90, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+ylim(c(5,-5))+geom_hline(yintercept = 0,lty=2,col="darkgrey")+guides(fill=F)+
    xlab("")+ylab("Posterior z-score")
 
 pc_fig<-ggplot()+geom_boxplot(data=pc_zscore,aes(y=X2,x=value,fill=X2),alpha=0.2)+theme_classic()+xlim(c(0,1))+geom_vline(xintercept = 0.5,lty=2,col="darkgrey")+guides(fill=F)+
    ylab("")+xlab("Posterior contraction")
 ggarrange(zscore_fig,pc_fig)
 
 #4. Posterior predictive checks...............................................

joined_sim <- rstan::extract(fit_sbc_nreref)
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
#save.image(file='ZamborainMasonetal_ReefSustainability_workflow3.RData')#null model
#save.image(file='ZamborainMasonetal_ReefSustainability_workflow3_full4.RData')#null model
 load(file='ZamborainMasonetal_ReefSustainability_workflow3_full4.RData') 
 #load(file='ZamborainMasonetal_ReefSustainability_workflow.RData')