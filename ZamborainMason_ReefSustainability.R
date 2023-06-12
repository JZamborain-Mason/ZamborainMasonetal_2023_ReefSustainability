# Title: Reef sustainability
# Author: Jessica Zamborain Mason
# Description: This code implements the analyses for "Sustainable reference points for multispecies coral reef fisheries"
# R version 4.2.1 (2022-06-23)
#Note that, given the nature of our models, results may vary slightly from run to run

#clear R
rm(list=ls())


#set working directory (where the data and code of this repository is stored)
#setwd("")
#setwd("C:/Users/jzamb/Dropbox/PhD/chapter 1/ALL ANALYSES/FINAL CODE AND DATA/FINAL AFTER REVISIONS (GITHUB)/CODE AFTER REVISIONS 2/CODE AFTER REVISIONS FROM SCIENCE/REVISITIONS AFTER NATURE ROUND 2")
setwd("C:/Users/jez297/Dropbox/PhD/chapter 1/Manuscript/NATURE COMMUNICATIONS/REVIEW/ZamborainMason_Data and Code_REVISED")


#load required libraries
library(ggrepel)  #Kamil Slowikowski (2018). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.8.0. https://CRAN.R-project.org/package=ggrepel
library(ggplot2) #H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
library(ggpubr) #  Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.2.https://CRAN.R-project.org/package=ggpubr
library(plyr) #Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
library(rstan) #Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. http://mc-stan.org/.
library(lme4) #  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
library(rstan)#Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. http://mc-stan.org/.
library(dplyr)#Hadley Wickham, Romain Fran?ois, Lionel Henry and Kirill M?ller (2019). dplyr: A Grammar of Data Manipulation. R package version 0.8.0.1. https://CRAN.R-project.org/package=dplyr
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
library(loo)#Vehtari A, Gabry J, Yao Y, Gelman A (2019). "loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models." R package version 2.1.0, <URL: https://CRAN.R-project.org/package=loo>.
library(geoR)#Paulo J. Ribeiro Jr and Peter J. Diggle (2018). geoR: Analysis of Geostatistical Data. R package version 1.7-5.2.1. https://CRAN.R-project.org/package=geoR
library(DHARMa)#Florian Hartig (2019). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.2.4. https://CRAN.R-project.org/package=DHARMa
library(brms)# Paul-Christian B?rkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of StatisticalSoftware, 80(1), 1-28. doi:10.18637/jss.v080.i01
library(ggeffects)#L?decke D (2018). "ggeffects: Tidy Data Frames of Marginal Effects from Regression Models." _Journal of Open SourceSoftware_, *3*(26), 772. doi: 10.21105/joss.00772 (URL: http://doi.org/10.21105/joss.00772).
library (loo)#Vehtari A, Gabry J, Magnusson M, Yao Y, Bürkner P, Paananen T, Gelman A (2022). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.5.1, <https://mc-stan.org/loo/>.

#to run chains in parallel
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
#backup_options <- options()
#options(backup_options)

#upload data (data is available in the paper's supplementary information)
reefscale_data<- read.csv("reefscale_data_revised.csv", head=T)
jurisdictionscale_data<- read.csv("jurisdictionscale_data_revised.csv", head=T)
jurisdictionscale_data2<-jurisdictionscale_data


#functions used in this script.........................................................................................................

#mean center covariates: standardise (Following Gelman and Hill 2007)
standardise <- function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#Pearson's correlation and histogram for pairs plot
panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*1.5)
}

#histogram for pairs
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

#variance inflation factors for hierarchical models
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

#Function to check missingness
pMiss <- function(x){sum(is.na(x))/length(x)*100}

#functions from Vehtari et al  2019 to monitor convergence
#obtain and use these documents through https://github.com/avehtari/rhat_ess
#if you use please cite the paper: Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, Paul-Christian B?rkner (2019): Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
source('monitornew_Vehtarietal2019.R')
source('monitorplot_Vehtarietal2019.R')

#number of sites sampled per jurisidiction
sites<- ddply(reefscale_data,.(Larger),summarize, reefsites=length(UniqueSite)) #bonaire is separated from neteherlands antilles but for reef area and catch it is the same
reefscale_data$Larger<- as.factor(ifelse(reefscale_data$Larger=="Bonaire","Netherlands Antilles", as.character(reefscale_data$Larger)))
sites<- ddply(reefscale_data,.(Larger),summarize, reefsites=length(UniqueSite)) 

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
         sgravtot=log(Grav_tot+1),#not standardized because we do not want potential exports to influence our estimated parameters
         stt_market=standardise(log(TT_market_h)))

#example average conditions
mean(sqrt(reefscale_data$Depth))^2
mean(sqrt(reefscale_data$HardCoral),na.rm=T)^2
exp(mean(log(reefscale_data$SampArea)))

#map of our sites
newmap <- getMap(resolution = "high")
#jitter points 
reefscale_data$Site_Lat2<- reefscale_data$Site_Lat+runif(length(reefscale_data$Site_Lat), min=0, max=3)
reefscale_data$Site_Long2<- reefscale_data$Site_Long+runif(length(reefscale_data$Site_Long), min=0, max=4)
reefscale_data$Site_Lat2<- ifelse(reefscale_data$Site_Lat2>23.5, reefscale_data$Site_Lat,reefscale_data$Site_Lat2)

#map of defined protection
#first make sure that MGMT_type as defined in MacNeil et al. 2015, matches our remote classification: re-clasify those he had defined as remote but that were <20h from human settlements
unique(reefscale_data$MGMT_type[reefscale_data$dataset=="M"])
unique(reefscale_data$MGMT_type[reefscale_data$dataset=="M"&reefscale_data$Management=="Remote" & !(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0)])
reefscale_data$definedprotection<- ifelse(reefscale_data$section=="benchmark"& reefscale_data$reserves==1,"HC Marine reserves", ifelse(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,"Remote",ifelse(reefscale_data$Management=="Restricted"|(reefscale_data$Management=="Remote" & !(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0)&reefscale_data$MGMT_type=="Access restricted"),"Restricted","Openly fished")))
ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=reefscale_data,aes(x=Site_Long2, y=Site_Lat2, fill = as.factor(reefscale_data$definedprotection)),colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="HC Marine reserves",],aes(x=Site_Long2, y=Site_Lat2), fill ="darkorchid1",colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="Remote",],aes(x=Site_Long2, y=Site_Lat2), fill ="green",colour="black", pch=21,size=3, alpha=0.7)+
  scale_fill_manual (name="Protection",values=c( "Openly fished"="red","Restricted"="turquoise1","HC Marine reserves"="darkorchid1", "Remote"="green"))+geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.line =element_blank())

#create dummy varibles for the different categorical variables
reefscale_data$pointintercept<- ifelse(reefscale_data$SampMethod=="Point intercept", 1,0)
reefscale_data$backreef<- ifelse(reefscale_data$ReefHabitat=="Lagoon_Back reef", 1,0)
reefscale_data$crest<- ifelse(reefscale_data$ReefHabitat=="Crest", 1,0)
reefscale_data$flat<- ifelse(reefscale_data$ReefHabitat=="Flat", 1,0)
reefscale_data$distancesampling<- ifelse(reefscale_data$SampMethod=="Distance sampling",1,0)
reefscale_data$model_component<- ifelse(reefscale_data$section=="benchmark"&reefscale_data$reserves==1,1,ifelse(reefscale_data$section=="benchmark"&reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,2,3))

#MacNeil et al. 2015 data defined as remote that does not satisfy our remote is classified as restricted (overlaps with some data providor classifications) 
reefscale_data$status_restricted<- ifelse(reefscale_data$model_component==3 &reefscale_data$definedprotection=="Restricted",1,0)
reefscale_data$status_fished<- ifelse(reefscale_data$model_component==3 &reefscale_data$definedprotection=="Fished",1,0)
reefscale_data$status_management<- ifelse(reefscale_data$status_restricted==1,"Restricted",ifelse(reefscale_data$status_fished==1,"Fished", reefscale_data$definedprotection))

#we do not have reserve size/coral cover for most sites so we  create variable to add 0 to those components (coral cover for the post-model estimation: average conditions)
reefscale_data$sHardCoral2<- ifelse(is.na(reefscale_data$sHardCoral),0,reefscale_data$sHardCoral)
reefscale_data$sClosure.size2<- ifelse(is.na(reefscale_data$sClosure.size),0,reefscale_data$sClosure.size)
reefscale_data$Closure.age2<- ifelse(is.na(reefscale_data$Closure.age),0,reefscale_data$Closure.age)

#correlation among covariates (check multicolinearity)
windows()
pairs(~ sDepth+ ReefHabitat+
        sSampArea+
        SampMethod+sOcean_prod+Atoll+sHardCoral2+sClosure.size2+sSST+Closure.age2+sgravtot,data=reefscale_data,lower.panel=panel.cor)

VIF.table.method<- as.data.frame(vif.mer(lmer(log(FamBiomass_tkm2)~
                                                sDepth+
                                                ReefHabitat+
                                                sSampArea+sSST+
                                                SampMethod+sOcean_prod+Atoll+sHardCoral2+sClosure.size2+Closure.age2+sgravtot+
                                                (1|Larger),data=reefscale_data)))
colnames(VIF.table.method)<- "VIF"
print(VIF.table.method)

#Overall no colinearity(whole dataset).
#none of those variables have multicolinearity with entire dataset
#write.csv(VIF.table.method, "VIF_refpointmodel.csv")

#separate reserves, remote and fished data and add index for random effects
fished_data<- reefscale_data[reefscale_data$model_component==3,]
fished_data<- droplevels(fished_data)
reserves_complete<- reefscale_data[reefscale_data$model_component==1,]
reserves_complete<- droplevels(reserves_complete)
remote_complete<- reefscale_data[reefscale_data$model_component==2,]
remote_complete<- droplevels(remote_complete)


#however, because different submodels inform different data subsets, I check potential colinearity
#between factors only included in one submodel
windows()
pairs(~ sDepth+ ReefHabitat+
        sSampArea+
        SampMethod+sOcean_prod+Atoll+sHardCoral+sSST+sClosure.size+Closure.age+sgravtot,data=reserves_complete,lower.panel=panel.cor)


pairs(~ sDepth+
        ReefHabitat+sSST+
        SampMethod+sOcean_prod+Atoll+sHardCoral+sSampArea,data=remote_complete,lower.panel=panel.cor)

pairs(~ sDepth+ ReefHabitat+
        sSampArea+
        SampMethod+sOcean_prod+Atoll+sHardCoral2+sSST+sgravtot,data=fished_data,lower.panel=panel.cor)


ggplot(rbind(reserves_complete,remote_complete),aes(x=as.factor(Atoll),y=sOcean_prod))+geom_boxplot()

#index for random effects 
fished_data$indexj<- as.numeric(factor(fished_data$Larger))
remote_complete$indexj<- as.numeric(factor(remote_complete$Larger))
reserves_complete$indexj<- as.numeric(factor(reserves_complete$Larger))

#data benchmark used in original model
refpointsdata<- rbind(reserves_complete,remote_complete)
#visualize reserve and remote data
a<- ggplot(reserves_complete, aes(x=Closure.age, y=FamBiomass_tkm2))+geom_point(aes(shape=reserves_complete$Region),col="darkgrey",alpha=0.5)+guides(shape=F)+theme_classic()+ labs(y = expression ("Biomass ("~t/km^2*")"))
b<- ggplot(remote_complete, aes(x=Locality, y=log(FamBiomass_tkm2)))+geom_boxplot(fill="grey", col="black",alpha=0.5)+geom_jitter(width = 0.2, height = 0,col="darkgrey",alpha=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill="white",colour="black"))+ labs(y = expression ("Biomass (log("~t/km^2*"))"))
ggarrange(a,b, nrow=1, ncol=2)

#create new data to plot model predictions
#As we have standardized all continuous variables, we create newdata  based only on reserve age (i.e., at average conditions)
newdata<- with(reserves_complete, data.frame(Closure.age=seq(min(Closure.age), max(Closure.age), len=length(reserves_complete$Closure.age)))) 
Bio<- seq(0,350,1) #biomass data for surplus production curve
B<- length(Bio)


#plot of sites
mapsites<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(-35, 35), ratio = 1.3)+
  geom_point(data=reefscale_data,aes(x=Site_Long2, y=Site_Lat2, fill = as.factor(reefscale_data$definedprotection)),colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="HC Marine reserves",],aes(x=Site_Long2, y=Site_Lat2), fill ="yellow",colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="Remote",],aes(x=Site_Long2, y=Site_Lat2), fill ="green",colour="black", pch=21,size=3, alpha=0.7)+
  scale_fill_manual (name="Protection",values=c( "Openly fished"="red","Restricted"="turquoise1","HC Marine reserves"="yellow", "Remote"="green"))+geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.line =element_blank(),panel.background = element_rect(fill="white",color="black"))

#check no major differences in environmental factors among management categories
a<-ggplot(reefscale_data,aes(fill=Management))+geom_jitter(aes(x=Management,y=sHardCoral),alpha=0.1,pch=21)+geom_violin(aes(x=Management,y=sHardCoral),alpha=0.3,draw_quantiles = c( 0.5))+ylab("Std. hard coral cover")+xlab("")+guides(fill=F)+theme_classic()+
  scale_fill_manual (name="",values=c( "Fished"="red","Restricted"="turquoise1","UnfishedHigh"="yellow", "Remote"="green"))
b<-ggplot(reefscale_data,aes(fill=Management))+geom_jitter(aes(x=Management,y=sOcean_prod),alpha=0.1,pch=21)+geom_violin(aes(x=Management,y=sOcean_prod),alpha=0.3,draw_quantiles = c( 0.5))+ylab("Std. ocean productivity")+xlab("")+guides(fill=F)+theme_classic()+
  scale_fill_manual (name="",values=c( "Fished"="red","Restricted"="turquoise1","UnfishedHigh"="yellow", "Remote"="green"))
c<-ggplot(reefscale_data,aes(fill=Management))+geom_jitter(aes(x=Management,y=sSST),alpha=0.1,pch=21)+geom_violin(aes(x=Management,y=sSST),alpha=0.3,draw_quantiles = c( 0.5))+ylab("Std. SST")+xlab("")+guides(fill=F)+theme_classic()+
  scale_fill_manual (name="",values=c( "Fished"="red","Restricted"="turquoise1","UnfishedHigh"="yellow", "Remote"="green"))
d<-ggplot(reefscale_data,aes(fill=Management))+geom_jitter(aes(x=Management,y=Atoll),alpha=0.5,pch=21)+xlab("")+guides(fill=F)+theme_classic()+
  scale_fill_manual (name="",values=c( "Fished"="red","Restricted"="turquoise1","UnfishedHigh"="yellow", "Remote"="green"))


#check that distribution of species richness in high complaince marine reserves is within those of fished sites
windows()
ggplot(reefscale_data[!is.na(reefscale_data$total_sp),],aes(fill=Management))+geom_jitter(aes(x=Management,y=log(total_sp)),alpha=0.5,pch=21)+geom_violin(aes(x=Management,y=log(total_sp)),alpha=0.3,draw_quantiles = c( 0.5))+ylab("log(species richness)")+xlab("")+guides(fill=F)+theme_classic()+
  scale_fill_manual (name="",values=c( "Fished"="red","Restricted"="turquoise1","UnfishedHigh"="yellow", "Remote"="green"))

windows()
envfig<-ggarrange(a,b,c,d, nrow=1,ncol=4,widths = c(1,1,1,1),labels=c("b","c","d","e"))
ggarrange(mapsites,envfig,nrow=2,ncol=1,labels=c("a",""))

#######################################################################################################################################
## Run reference point and status model  ..............................................................................................

#data
stanDat_full_country<- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                            res=nrow(reserves_complete), d=reserves_complete$sDepth, hc=reserves_complete$sHardCoral, si=reserves_complete$sClosure.size,
                            at=reserves_complete$Atoll,sst=reserves_complete$sSST,
                            op=reserves_complete$sOcean_prod,rh_c=reserves_complete$crest,rh_b=reserves_complete$backreef,rh_f=reserves_complete$flat,
                            cm_pc=reserves_complete$pointintercept, sa=reserves_complete$sSampArea, gr=reserves_complete$sgravtot,
                            b2 = log(remote_complete$FamBiomass_tkm2),rem=nrow(remote_complete), d2=remote_complete$sDepth, hc2=remote_complete$sHardCoral, 
                            at2=remote_complete$Atoll,sst2=remote_complete$sSST,
                            op2=remote_complete$sOcean_prod,rh_c2=remote_complete$crest,rh_b2=remote_complete$backreef,
                            sa2=remote_complete$sSampArea,
                            cm_ds2=remote_complete$distancesampling,
                            b3 = log(fished_data$FamBiomass_tkm2),fis=nrow(fished_data), d3=fished_data$sDepth, 
                            at3=fished_data$Atoll,sst3=fished_data$sSST,gr3=fished_data$sgravtot,
                            op3=fished_data$sOcean_prod,rh_c3=fished_data$crest,rh_b3=fished_data$backreef,rh_f3=fished_data$flat,
                            cm_pc3=fished_data$pointintercept, sa3=fished_data$sSampArea,cm_ds3=fished_data$distancesampling,hc3=fished_data$sHardCoral2,
                            m_r3=fished_data$status_restricted,
                            R3=nlevels(factor(fished_data$Larger)),
                            pr3=fished_data$indexj,
                            newage=newdata$Closure.age,
                            Bio=Bio, B=B)

#null model
Fit_null_gf<- stan(file = "Null_gomperztfox_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
windows()
pairs(Fit_null_gf,pars=c("r","B0","bmin"))

#full model
Fit_full_gf<- stan(file = "Full_gompertzfox_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
#save(Fit_full_gf,file="Fit_full_inreview.RData")
windows()
pairs(Fit_full_gf,pars=c("log_r","B0","log_bmin"))


#model selection through loo: favours full model
full_loglik<- extract_log_lik(Fit_full_gf, merge_chains = F)
r_eff_full <- relative_eff(exp(full_loglik)) 
loo_full <- loo(full_loglik,r_eff =r_eff_full)

null_loglik<- extract_log_lik(Fit_null_gf, merge_chains = F)
r_eff_null <- relative_eff(exp(null_loglik)) 
loo_null <- loo(null_loglik,r_eff =r_eff_null)
comp <- loo_compare(loo_null, loo_full)
print(comp, simplify=F)


#estimate posterior contraction for model
var_prior_logB0<-(1^2)
var_prior_logbmin<-(1^2)
var_prior_logr<-(1^2)
var_prior_beta<-(2^2)
var_priorI_fished<-(5^5)
var_prior_sigma<-(sd(rcauchy(4000,0,1))^2)
var_prior_MMSY<-(sd(exp(rnorm(4000,-2,1))*exp(rnorm(4000,log(120),1))/2.718281828)^2)
var_prior_BMMSY<-(sd(exp(rnorm(4000,log(120),1))/2.718281828)^2)
contraction_full<-as.data.frame(cbind(round(c((var_prior_logB0-(sd(rstan::extract(Fit_full_gf,pars=c("log_B0"))$log_B0)^2))/var_prior_logB0,(var_prior_logr-(sd(rstan::extract(Fit_full_gf,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_logbmin-(sd(rstan::extract(Fit_full_gf,pars=c("log_bmin"))$log_bmin)^2))/var_prior_logbmin,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,1])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,2])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,3])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,4])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,5])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,6])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,7])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,8])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,9])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,10])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,11])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,12])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf,pars=c("beta"))$beta[,13])^2))/var_prior_beta,(var_priorI_fished-(sd(rstan::extract(Fit_full_gf,pars=c("I_fished"))$I_fished)^2))/var_priorI_fished,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf,pars=c("sigma_e"))$sigma_e)^2))/var_prior_sigma,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf,pars=c("sigma_r"))$sigma_r)^2))/var_prior_sigma,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf,pars=c("sigma_f"))$sigma_f)^2))/var_prior_sigma,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf,pars=c("sigma_u"))$sigma_u)^2))/var_prior_sigma,(var_prior_MMSY-(sd(rstan::extract(Fit_full_gf,pars=c("MMSY"))$MMSY)^2))/var_prior_MMSY,(var_prior_BMMSY-(sd(rstan::extract(Fit_full_gf,pars=c("BMMSY"))$BMMSY)^2))/var_prior_BMMSY),2),c("log_B0","log_r","log_bmin","b_op","b_sst","b_hc","b_at","b_rhc","b_rhf","b_rhbr","b_cmpc","b_sa","b_rs","b_cmds","b_d","b_grav","I_fished", "sd_res","sd_rem","sd_fis","sd_u","MMSY","BMMSY")))
colnames( contraction_full)<-c("posterior contraction","parameter")
print(contraction_full)


#posterior draws of best-ranked model
list_of_draws_full <- as.data.frame(Fit_full_gf) #rows are iterations and columns are parameters
#write.csv(list_of_draws_full,"list_of_draws_full.csv")

#now we look at variograms on random effects to see if there is strong spatial structure even when we include environmental covariates
fished_Larger<- ddply(fished_data,.(Larger),summarize, Site_Lat=Site_Lat[1], Site_Long=Site_Long[1]) 
#random effects
randomef<- c(matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("u3"))$u3))
vardata <- data.frame(x = c(fished_Larger$Site_Lat), y = c(fished_Larger$Site_Long), resid = randomef) 
geodat <- geoR::as.geodata(vardata)
geodat_variog <- geoR::variog(geodat, breaks = seq(0, 310, 5),option = "bin",
                        estimator.type =  "classic")
fit_geodat_variog<- geoR::variofit(geodat_variog )
#variogram random effects: no spatial structure
windows()
plot(geodat_variog,pts.range = c(1,3)); lines(fit_geodat_variog, col='red')

#model diagnostics
mon_full_all_open<- monitor(Fit_full_gf)
print(mon_full_all_open)

#rank plots 
samp_full_open<- as.array(Fit_full_gf)
mcmc_hist_r_scale(samp_full_open[, , "MMSY"])
mcmc_hist_r_scale(samp_full_open[, , "BMMSY"])
modeldiagplot1<- stan_trace(Fit_full_gf, pars=c("MMSY","BMMSY"))
modeldiagplot2<- ggplot(mon_full_all_open, aes(x=Rhat))+geom_histogram()+xlab("Rhat statistic")+theme_classic()+geom_vline(xintercept = 1.01,lty=2)+geom_text(aes(x=1.0101,y=2000),label="1.01")
modeldiagplot3<- ggplot(mon_full_all_open, aes(x=Tail_ESS))+geom_histogram()+xlab("Tail effective sample sizes")+theme_classic()+geom_vline(xintercept =400,lty=2)+geom_text(aes(x=410,y=2000),label="400")
diagnostics<- ggarrange(modeldiagplot1,modeldiagplot2, modeldiagplot3,nrow=1,ncol=3, labels=c("a","b","c"), widths = c(1.5,1,1))


#posterior vs prior
BO_prior<- exp(rnorm(4000, log(120),1))
bmin_prior<- exp(rnorm(4000, log(10),1))
r_prior<- exp(rnorm(4000, -2,1))
beta_prior<- rnorm(4000, 0,2)
MMSY_prior<- BO_prior*r_prior/2.718281828
BMMSY_prior<- BO_prior/2.718281828

c<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$B0),fill="black",col="black")+geom_histogram(aes(x=BO_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Unfished Biomass")
d<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$bmin),fill="black",col="black")+geom_histogram(aes(x=bmin_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Biomass reserve age 0")+ylab("")
e<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$r ),fill="black",col="black")+geom_histogram(aes(x=r_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("community growth rate")+ylab("")
f<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$`beta[1]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[2]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[3]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[4]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[5]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[6]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[7]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[8]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[9]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[10]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[11]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[12]`),fill="black",col="black")+
  geom_histogram(aes(x=list_of_draws_full$`beta[13]`),fill="black",col="black")+
  geom_histogram(aes(x=beta_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Effect sizes")+ylab("")
g<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$MMSY),fill="black",col="black")+geom_histogram(aes(x=MMSY_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("MMSY")
h<- ggplot(NULL)+geom_histogram(aes(x=list_of_draws_full$BMMSY),fill="black",col="black")+geom_histogram(aes(x=BMMSY_prior), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("BMMSY")+ylab("")
postprior<- ggarrange(c,d,e,f,g,h,nrow=1,ncol=6, labels=c("d","e","f","g","h","i"))

#model fit
Fit_full_gf_summary <- summary(Fit_full_gf,probs=c(0.1,0.25,0.5,0.75, 0.9))
output_Fit_full_gf<- as.data.frame(Fit_full_gf_summary$summary)

#examining model fit
pred_full<- c(matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("mu2"))$mu2),matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("mu3"))$mu3))
resid_full<- c(log(reserves_complete$FamBiomass_tkm2),log(remote_complete$FamBiomass_tkm2), log(fished_data$FamBiomass_tkm2))-pred_full
a_full_open<- ggplot(data=NULL,aes(x=pred_full,y=resid_full))+geom_point()+theme_classic()+ggtitle("")+xlab("Fitted ")+ylab("Residuals ")
b_full_open<- ggplot(NULL, aes(x = resid_full)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("")+ylab("Count")+xlab("Residuals ")

#posterior predictive checks
joined_sim <- rstan::extract(Fit_full_gf)
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
a<- bayesplot::ppc_dens_overlay(log(reserves_complete$FamBiomass_tkm2),y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b<- bayesplot::ppc_dens_overlay(log(remote_complete$FamBiomass_tkm2),y_rep_remote[1:4000,])+ggtitle("Remote")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))+guides(col=F)
c<- bayesplot::ppc_dens_overlay(log(fished_data$FamBiomass_tkm2),y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~t/km^2*"))"))

resid_fig<- ggarrange(a_full_open,b_full_open,nrow=1, ncol=2,labels=c("j","k"))
ppcheckfig<- ggarrange(a,b,c,nrow=1,ncol=3, widths=c(1,1,1.2),labels=c("l","m","n"))
windows()
ggarrange(diagnostics,postprior,resid_fig,ppcheckfig,nrow=4,ncol=1)

## Reference point posteriors .........................................................................................................
#effect sizes from covariates
betas<- output_Fit_full_gf[1:13,c("50%","10%","90%")]
betas$variable<- c("Ocean productivity", "SST", "Hard coral","Atoll","Crest", "Flat", "Backreef/lagoon","Point count","Sampling area","Reserve size", "Distance sampling","Depth","Gravity")
betas$order<- c(1,4,2,3,6,7,8,9,10,11,12,5,13)
betas[order(betas$order),]
betas$variable <- factor(betas$variable, levels = betas$variable[order(betas$order)])
betasfig2<- ggplot(betas,aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`))+
  geom_pointrange(size=0.7,col="navyblue")+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("Effect size")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 90, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))

#model posteriors
#obtain median 50% and 90% credible intervals for parameters of interest
r_full_open<- c(median(list_of_draws_full$r), mean(list_of_draws_full$r), quantile(list_of_draws_full$r,  probs=c(0.1,0.25,0.75, 0.9)))
B0_full_open<- c(median(list_of_draws_full$B0), mean(list_of_draws_full$B0), quantile(list_of_draws_full$B0,  probs=c(0.1,0.25,0.75, 0.9)))
bmin_full_open<- c(median(list_of_draws_full$bmin), mean(list_of_draws_full$bmin), quantile(list_of_draws_full$bmin,  probs=c(0.1,0.25,0.75, 0.9)))
bmmsy_full_open<- c(median(list_of_draws_full$BMMSY), mean(list_of_draws_full$BMMSY), quantile(list_of_draws_full$BMMSY,  probs=c(0.1,0.25,0.75, 0.9)))
mmsy_full_open<- c(median(list_of_draws_full$MMSY), mean(list_of_draws_full$MMSY), quantile(list_of_draws_full$MMSY,  probs=c(0.1,0.25,0.75, 0.9)))

Fit_full_gf_open_pars<- cbind(r_full_open,B0_full_open,bmin_full_open,bmmsy_full_open,mmsy_full_open)
rownames(Fit_full_gf_open_pars)=c("median","mean","10%","25%","75%","90%")
print(Fit_full_gf_open_pars)
#write.csv(Fit_full_gf_open_pars, "pars_bestfirmodel_ss_cen.csv", row.names=F)

#plot median reserve trajectory with 90 credible intervals
newpred_full_open<- rstan::extract(Fit_full_gf,pars=c("predreserveB"))[[1]]
newpred_full_open<- broom.mixed::tidyMCMC(coda::as.mcmc(newpred_full_open),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
newpred_full_open$Closure.age=newdata$Closure.age
full_tra<- ggplot(NULL)+
  geom_point(data=reserves_complete, aes(y= FamBiomass_tkm2, x=Closure.age),col="darkgrey",fill="darkgrey", alpha=0.35)+
  geom_line(aes(y=newpred_full_open$estimate, x=newpred_full_open$Closure.age), col="navyblue",lwd=2)+geom_ribbon(aes(x=newpred_full_open$Closure.age,ymin=newpred_full_open$conf.low, ymax=newpred_full_open$conf.high), alpha=0.3, fill="blue")+theme_classic()+xlab("MPA age (years)")+ylab("Biomass (t/km2)")+ labs(y = expression ("Biomass ("~t/km^2*")"))
full_tra<- ggplot(NULL)+
  geom_point( aes(y= matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("site_Bmarg"))[[1]])[1:nrow(reserves_complete)], x=reserves_complete$Closure.age),col="darkgrey",fill="darkgrey", alpha=0.35)+
  geom_line(aes(y=newpred_full_open$estimate, x=newpred_full_open$Closure.age), col="navyblue",lwd=2)+geom_ribbon(aes(x=newpred_full_open$Closure.age,ymin=newpred_full_open$conf.low, ymax=newpred_full_open$conf.high), alpha=0.3, fill="blue")+theme_classic()+xlab("MPA age (years)")+ylab("Biomass (t/km2)")+ labs(y = expression ("Biomass ("~t/km^2*")"))

#plot posterior unfished biomass with 90% credible intervals on remote boxplot
full_remote<- ggplot(NULL)+geom_boxplot(data=remote_complete, aes(x=Locality, y=log(FamBiomass_tkm2)),fill="darkgrey", alpha=0.35)+geom_jitter(data=remote_complete, aes(x=Locality, y=log(FamBiomass_tkm2)),width = 0.2, height = 0,col="darkgrey",alpha=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+geom_abline(intercept=log(Fit_full_gf_open_pars[1,2]), slope=0, col="navyblue",lwd=2)+
  geom_rect(aes(xmin=0, xmax=Inf,ymin=log(Fit_full_gf_open_pars[3,2]),ymax=log(Fit_full_gf_open_pars[6,2])),fill="blue",alpha=0.3)+ labs(y = expression ("Biomass (log("~t/km^2*"))"))
full_remote<- ggplot(NULL)+geom_boxplot( aes(x=remote_complete$Locality, y=log(matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("site_Bmarg"))[[1]])[(nrow(reserves_complete)+1):(nrow(reserves_complete)+nrow(remote_complete))])),fill="darkgrey", alpha=0.35)+geom_jitter(aes(x=remote_complete$Locality, y=log(matrixStats::colMedians(rstan::extract(Fit_full_gf,pars=c("site_Bmarg"))[[1]])[(nrow(reserves_complete)+1):(nrow(reserves_complete)+nrow(remote_complete))])),width = 0.2, height = 0,col="darkgrey",alpha=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+geom_abline(intercept=log(Fit_full_gf_open_pars[1,2]), slope=0, col="navyblue",lwd=2)+
  geom_rect(aes(xmin=0, xmax=Inf,ymin=log(Fit_full_gf_open_pars[3,2]),ymax=log(Fit_full_gf_open_pars[6,2])),fill="blue",alpha=0.3)+ labs(y = expression ("Biomass (log("~t/km^2*"))"))


#plot surplus production curve along a gradient of biomass
popgrowth_full_open<- rstan::extract(Fit_full_gf,pars=c("sustyield"))[[1]]
popgrowth_full_open[is.nan(popgrowth_full_open)] <- 0
popgrowth_full_open<- broom.mixed::tidyMCMC(coda::as.mcmc(popgrowth_full_open),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
#assign almost 0 surplus for biomass values that have no surplus
popgrowth_full_open$median=ifelse(popgrowth_full_open$estimate<0,0.00000001,popgrowth_full_open$estimate)
popgrowth_full_open$Biomass=Bio
#places that should have no surplus production assign very small (to not get inf)
popgrowth_full_open$`5%`<- ifelse(popgrowth_full_open$conf.low<0,0.000001,popgrowth_full_open$conf.low)
popgrowth_full_open$`95%`<- ifelse(popgrowth_full_open$conf.high<0,0.000001,popgrowth_full_open$conf.high)

#for plotting purposes
popgrowth_full_open2<- popgrowth_full_open
popgrowth_full_open2$median<- ifelse(popgrowth_full_open2$median==0.00000001,NA,popgrowth_full_open2$median)
popgrowth_full_open2$`5%`<- ifelse(popgrowth_full_open2$`5%`==0.000001,NA,popgrowth_full_open2$`5%`)
popgrowth_full_open2$`95%`<- ifelse(popgrowth_full_open2$`95%`==0.000001,NA,popgrowth_full_open2$`95%`)

#shading regions
shade1<- rbind(c(0,0), subset(popgrowth_full_open2, Biomass >= 
                                0 & Biomass <= popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)]), c(popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)], 0))
shade3<- rbind(c(0,Inf), subset(popgrowth_full_open2, Biomass >= 
                                  0 & Biomass <=popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)]), c(popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)], Inf))
#surplus curve
full_grow1<- ggplot(NULL)+geom_polygon(data=shade1,aes(Biomass, median), fill="cyan3",alpha=0.7)+
  geom_rect(aes(xmin=popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)], xmax=Inf,ymin=0,ymax=max(popgrowth_full_open$median)),fill="navyblue",alpha=0.6)+geom_polygon(data=shade3,aes(Biomass, median), fill="red3",alpha=0.7)+
  geom_rect(aes(xmin=popgrowth_full_open$Biomass[popgrowth_full_open$median==max(popgrowth_full_open$median)], xmax=Inf,ymin=max(popgrowth_full_open$median)),ymax=Inf, fill="goldenrod1",alpha=0.7)+
  geom_line(aes(y=popgrowth_full_open2$median[!is.na(popgrowth_full_open2$median)], x=popgrowth_full_open2$Biomass[!is.na(popgrowth_full_open2$median)]), col="black", lwd=2)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression (atop("Potential sustainable", paste("yield ("~t/km^2/y*")"))))+scale_x_continuous(expand = c(0, 0),limits=c(0,180)) + scale_y_continuous(expand = c(0, 0),limits=c(0,20))

full_grow<- ggplot(NULL)+geom_line(aes(y=popgrowth_full_open2$median[!is.na(popgrowth_full_open2$median)], x=popgrowth_full_open2$Biomass[!is.na(popgrowth_full_open2$median)]), col="navyblue", lwd=2)+
  geom_ribbon(aes(x=popgrowth_full_open$Biomass,ymin=popgrowth_full_open$`5%`[!is.na(popgrowth_full_open$`5%`)],ymax=popgrowth_full_open$`95%`[!is.na(popgrowth_full_open$`95%`)]), fill="blue", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))+ylim(c(0,20))+xlim(c(0,200))
modelfit1<- ggarrange(full_tra,full_remote, full_grow, nrow=1, ncol=3, labels=c("a","b","c"))

#unfished biomass, intrinsic growth rate and biomass at reserve age 0 posterior distributions
B0_full_open<- ggplot(list_of_draws_full)+geom_histogram(aes(x=B0, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("Unfished Biomass(t/km2)")+ylab("")+ labs(x = expression ("Unfished Biomass ("~t/km^2*")"))
r_full_open<- ggplot(list_of_draws_full)+geom_histogram(aes(x=log_r, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("")+ylab("Posterior density")+ labs(x = expression ("log(community growth rate (1/y))"))
bmin_full_open<- ggplot(list_of_draws_full)+geom_histogram(aes(x=bmin, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("Biomass age 0 (t/km2)")+ylab("")+ labs(x = expression ("Biomass age 0 ("~t/km^2*")"))
MMSY_full_open<- ggplot(list_of_draws_full)+geom_histogram(aes(x=MMSY, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("MMSY (t/km2/y)")+ylab("")+ labs(x = expression ("MMSY("~t/km^2/y*")"))
BMMSY_full_open<- ggplot(list_of_draws_full)+geom_histogram(aes(x=BMMSY, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("BMMSY (t/km2)")+ylab("")+ labs(x = expression ("BMMSY ("~t/km^2*")"))
logistic_par_post<- ggarrange(r_full_open,B0_full_open,bmin_full_open,MMSY_full_open,BMMSY_full_open,nrow=1, ncol=5, labels=c("d","e","f","g","h"))

#suplementary figure of posteriors
modelfit<- ggarrange(modelfit1, logistic_par_post,nrow=2, ncol=1, heights=c(1.5,1))
windows()
ggarrange(modelfit, betasfig2,nrow=2,ncol=1,heights = c(2,1),labels=c("","i"))

#Simulate data to show the estmated effect of environmental covariates on reference points (assuming all other covaruates at average values)
#coral cover 
coraltrends<- expand.grid(sHardCoral=seq(min(reefscale_data$sHardCoral,na.rm=T),max(reefscale_data$sHardCoral,na.rm=T),length=101),
                          Atoll=c(0,1))
X <- model.matrix(~sHardCoral+Atoll,data=coraltrends)
coefs<- cbind(log(list_of_draws_full$B0),list_of_draws_full$`beta[3]`,list_of_draws_full$`beta[4]`)
fit_B0<-exp(coefs%*% t(X))
B0estimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_B0),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
B0estimates$term<-NULL
colnames(B0estimates)<-c("median_B0","std.error_B0","conf.low_B0","conf.high_B0")
fit_BMMSY<-fit_B0/2.718281828
BMMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_BMMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
BMMSYestimates$term<-NULL
colnames(BMMSYestimates)<-c("median_BMMSY","std.error_BMMSY","conf.low_BMMSY","conf.high_BMMSY")
fit_MMSY<-fit_B0*list_of_draws_full$r/2.718281828
MMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_MMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
MMSYestimates$term<-NULL
colnames(MMSYestimates)<-c("median_MMSY","std.error_MMSY","conf.low_MMSY","conf.high_MMSY")
coraltrends<- coraltrends%>% cbind(B0estimates,BMMSYestimates,MMSYestimates)
#ocean productivity
prodtrends<- expand.grid(sOcean_prod=seq(min(reefscale_data$sOcean_prod,na.rm=T),max(reefscale_data$sOcean_prod,na.rm=T),length=101),
                         Atoll=c(0,1))
X <- model.matrix(~sOcean_prod+Atoll,data=prodtrends)
coefs<- cbind(log(list_of_draws_full$B0),list_of_draws_full$`beta[1]`,list_of_draws_full$`beta[4]`)
fit_B0<-exp(coefs%*% t(X))
B0estimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_B0),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
B0estimates$term<-NULL
colnames(B0estimates)<-c("median_B0","std.error_B0","conf.low_B0","conf.high_B0")
fit_BMMSY<-fit_B0/2.718281828
BMMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_BMMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
BMMSYestimates$term<-NULL
colnames(BMMSYestimates)<-c("median_BMMSY","std.error_BMMSY","conf.low_BMMSY","conf.high_BMMSY")
fit_MMSY<-fit_B0*list_of_draws_full$r/2.718281828
MMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_MMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
MMSYestimates$term<-NULL
colnames(MMSYestimates)<-c("median_MMSY","std.error_MMSY","conf.low_MMSY","conf.high_MMSY")
prodtrends<- prodtrends %>% cbind(B0estimates,BMMSYestimates,MMSYestimates)
#SST trends
ssttrends<- expand.grid(sSST=seq(min(reefscale_data$sSST,na.rm=T),max(reefscale_data$sSST,na.rm=T),length=101),
                        Atoll=c(0,1))
X <- model.matrix(~sSST+Atoll,data=ssttrends)
coefs<- cbind(log(list_of_draws_full$B0),list_of_draws_full$`beta[2]`,list_of_draws_full$`beta[4]`)
fit_B0<-exp(coefs%*% t(X))
B0estimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_B0),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
B0estimates$term<-NULL
colnames(B0estimates)<-c("median_B0","std.error_B0","conf.low_B0","conf.high_B0")
fit_BMMSY<-fit_B0/2.718281828
BMMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_BMMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
BMMSYestimates$term<-NULL
colnames(BMMSYestimates)<-c("median_BMMSY","std.error_BMMSY","conf.low_BMMSY","conf.high_BMMSY")
fit_MMSY<-fit_B0*list_of_draws_full$r/2.718281828
MMSYestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_MMSY),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
MMSYestimates$term<-NULL
colnames(MMSYestimates)<-c("median_MMSY","std.error_MMSY","conf.low_MMSY","conf.high_MMSY")
ssttrends<- ssttrends %>% cbind(B0estimates,BMMSYestimates,MMSYestimates)
#atoll trends
MyData2<- expand.grid(Atoll=c(0,1))
X_2<- model.matrix(~Atoll,data=MyData2)
coefs2<- cbind(log(list_of_draws_full$B0),list_of_draws_full$`beta[4]`)
fit2<-coefs2%*% t(X_2)
fit_B02<-exp(fit2)
colnames(fit_B02)=c("Non-atoll","Atoll")
atolltrend<- melt(fit_B02)
fit_BMMSY2<-fit_B02/2.718281828
colnames(fit_BMMSY2)<- c("Non-atoll","Atoll")
atolltrend2<- melt(fit_BMMSY2)
fit_MMSY2<-fit_B02*list_of_draws_full$r/2.718281828
colnames(fit_MMSY2)<- c("Non-atoll","Atoll")
atolltrend3<- melt(fit_MMSY2)
colnames(atolltrend3)

#combine to create SI fig
a<- ggplot(coraltrends)+
  geom_line(aes(x=sHardCoral,y=as.numeric(median_B0),col=as.factor(coraltrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sHardCoral,ymax=conf.high_B0,ymin=conf.low_B0,fill=as.factor(coraltrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+labs(y = expression ("Unfished biomass ("~t/km^2*")"))+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+guides(fill=F,col=F)
a2<- ggplot(coraltrends)+
  geom_line(aes(x=sHardCoral,y=median_BMMSY,col=as.factor(coraltrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sHardCoral,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(coraltrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+labs(y = expression ("BMMSY ("~t/km^2*")"))+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+guides(fill=F,col=F)
a3<- ggplot(coraltrends)+
  geom_line(aes(x=sHardCoral,y=median_MMSY,col=as.factor(coraltrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sHardCoral,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(coraltrends$Atoll)),alpha=0.2)+theme_classic()+xlab("Stnd. coral cover")+labs(y = expression ("MMSY("~t/km^2/y*")"))+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+guides(fill=F,col=F)
b<- ggplot(prodtrends)+
  geom_line(aes(x=sOcean_prod,y=median_B0,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sOcean_prod,ymax=conf.high_B0,ymin=conf.low_B0,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
b2<- ggplot(prodtrends)+
  geom_line(aes(x=sOcean_prod,y=median_BMMSY,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sOcean_prod,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
b3<- ggplot(prodtrends)+
  geom_line(aes(x=sOcean_prod,y=median_MMSY,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sOcean_prod,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("Stnd. ocean productivity")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
c<- ggplot(ssttrends)+
  geom_line(aes(x=sSST,y=median_B0,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sSST,ymax=conf.high_B0,ymin=conf.low_B0,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
c2<- ggplot(ssttrends)+
  geom_line(aes(x=sSST,y=median_BMMSY,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sSST,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
c3<- ggplot(ssttrends)+
  geom_line(aes(x=sSST,y=median_MMSY,col=as.factor(prodtrends$Atoll)),size=1)+
  geom_ribbon(aes(x=sSST,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("Stnd. sea surface temperature")+
  scale_fill_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+
  scale_color_manual(values = c("1" = "mediumspringgreen",  "0"="navyblue"))+ylab("")+guides(fill=F,col=F)
d<- ggplot(atolltrend)+
  geom_violin(aes(x=X2,y=value,fill=X2,col=X2),draw_quantiles = c( 0.5),alpha=0.2)+guides(fill=F,col=F)+
  scale_colour_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+
  scale_fill_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+xlab("")+ylab("")+theme_classic()
d2<- ggplot(atolltrend2)+
  geom_violin(aes(x=X2,y=value,fill=X2,col=X2),draw_quantiles = c( 0.5),alpha=0.2)+guides(fill=F,col=F)+
  scale_colour_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+
  scale_fill_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+xlab("")+ylab("")+theme_classic()
d3<- ggplot(atolltrend3)+
  geom_violin(aes(x=X2,y=value,fill=X2,col=X2),draw_quantiles = c( 0.5),alpha=0.2)+guides(fill=F,col=F)+
  scale_colour_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+ylim(c(0,100))+
  scale_fill_manual(values = c("Non-atoll" = "navyblue",  "Atoll"="mediumspringgreen"))+xlab("(other environmental factors at mean)")+ylab("")+theme_classic()

B0fig<- ggarrange(a,b,c,d,nrow=1,ncol=4,labels=c("a","b","c","d"))
BMMSYfig<- ggarrange(a2,b2,c2,d2,nrow=1,ncol=4,labels=c("e","f","g","h"))
MMSYfig<- ggarrange(a3,b3,c3,d3,nrow=1,ncol=4,labels=c("i","j","k","l"))
trendsfig<- ggarrange(B0fig,BMMSYfig,MMSYfig,ncol=1,nrow=3)

#Atoll surplus
B0_atoll<- exp(log(list_of_draws_full$B0)+list_of_draws_full$`beta[4]`)
median(B0_atoll)
popgrowth_full_atoll<- matrix(NA, nrow=length(list_of_draws_full$r),ncol=length(Bio))

for (i in 1:length(Bio)){
  popgrowth_full_atoll[,i]<- log(B0_atoll)*list_of_draws_full$r*Bio[i] * (1-(log(Bio[i])/log(B0_atoll)))
}
popgrowth_full_atoll[is.nan(popgrowth_full_atoll)] <- 0
popgrowth_full_atoll<- broom.mixed::tidyMCMC(coda::as.mcmc(popgrowth_full_atoll),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
popgrowth_full_atoll$median<- ifelse(popgrowth_full_atoll$estimate<0,0.00000001,popgrowth_full_atoll$estimate)
popgrowth_full_atoll$Biomass<- Bio
popgrowth_full_atoll$`5%`<- ifelse(popgrowth_full_atoll$conf.low<0,0.000001,popgrowth_full_atoll$conf.low)
popgrowth_full_atoll$`95%`<- ifelse(popgrowth_full_atoll$conf.high<0,0.000001,popgrowth_full_atoll$conf.high)
#for plotting purposes
popgrowth_full_atoll2<- popgrowth_full_atoll
popgrowth_full_atoll2$median<- ifelse(popgrowth_full_atoll2$median==0.00000001,NA,popgrowth_full_atoll2$median)
popgrowth_full_atoll2$`5%`<- ifelse(popgrowth_full_atoll2$`5%`==0.000001,NA,popgrowth_full_atoll2$`5%`)
popgrowth_full_atoll2$`95%`<- ifelse(popgrowth_full_atoll2$`95%`==0.000001,NA,popgrowth_full_atoll2$`95%`)

#plot surplus production curve along a gradient of biomass
popgrowth_full_open<- rstan::extract(Fit_full_gf,pars=c("sustyield"))[[1]]
popgrowth_full_open[is.nan(popgrowth_full_open)] <- 0
popgrowth_full_open<- broom.mixed::tidyMCMC(coda::as.mcmc(popgrowth_full_open),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
#assign almost 0 surplus for biomass values that have no surplus
popgrowth_full_open$median=ifelse(popgrowth_full_open$estimate<0,0.00000001,popgrowth_full_open$estimate)
popgrowth_full_open$Biomass=Bio
#places that should have no surplus production assign very small (to not get inf)
popgrowth_full_open$`5%`<- ifelse(popgrowth_full_open$conf.low<0,0.000001,popgrowth_full_open$conf.low)
popgrowth_full_open$`95%`<- ifelse(popgrowth_full_open$conf.high<0,0.000001,popgrowth_full_open$conf.high)
#for plotting purposes
popgrowth_full_open2<- popgrowth_full_open
popgrowth_full_open2$median<- ifelse(popgrowth_full_open2$median==0.00000001,NA,popgrowth_full_open2$median)
popgrowth_full_open2$`5%`<- ifelse(popgrowth_full_open2$`5%`==0.000001,NA,popgrowth_full_open2$`5%`)
popgrowth_full_open2$`95%`<- ifelse(popgrowth_full_open2$`95%`==0.000001,NA,popgrowth_full_open2$`95%`)

atfig<- ggplot(NULL)+geom_line(aes(y=popgrowth_full_open2$median[!is.na(popgrowth_full_open2$median)], x=popgrowth_full_open2$Biomass[!is.na(popgrowth_full_open2$median)]), col="navyblue", lwd=2)+
  geom_ribbon(aes(x=popgrowth_full_open$Biomass,ymin=popgrowth_full_open$`5%`[!is.na(popgrowth_full_open$`5%`)],ymax=popgrowth_full_open$`95%`[!is.na(popgrowth_full_open$`95%`)]), fill="navyblue", alpha=0.3)+
  geom_line(aes(y=popgrowth_full_atoll2$median[!is.na(popgrowth_full_atoll2$median)], x=popgrowth_full_atoll2$Biomass[!is.na(popgrowth_full_atoll2$median)]), col="aquamarine3", lwd=2)+
  geom_ribbon(aes(x=popgrowth_full_atoll$Biomass,ymin=popgrowth_full_atoll$`5%`[!is.na(popgrowth_full_atoll$`5%`)],ymax=popgrowth_full_atoll$`95%`[!is.na(popgrowth_full_atoll$`95%`)]), fill="mediumspringgreen", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))+ylim(c(0,35))+xlim(c(0,300))+
  geom_text(aes(x=140,y=8),label="Atoll", col="mediumspringgreen",size=4.5)+
  geom_text(aes(x=70,y=4),label="Non-atoll", col="navyblue",size=4.5)

#supplementary figure
windows()  
ggarrange(trendsfig,atfig,ncol=2,nrow=1,widths=c(2.5,1),labels=c("","j"))

#check estimated bmin is within (and to the lower end) of fished reefs
windows()
ggplot()+geom_density(aes(x=list_of_draws_full$log_bmin),fill="darkred",alpha=0.5)+
  geom_density(aes(x=log(fished_data$FamBiomass_tkm2[fished_data$definedprotection=="Openly fished"])),fill="grey",alpha=0.5)+theme_classic()+xlab ("log-Biomass (t/km2))")+
  geom_text(aes(x=4,y=0.75),label="Bmin",col="darkred")+geom_text(aes(x=6,y=0.4),label="Openly fished sites",col="grey")+ylab("Density")

#######################################################################################################################################
#Reference points per jurisdiction and status of fished reefs .........................................................................

#extract model posterior
everything<- rstan::extract(Fit_full_gf)

#status for each reef site (note this includes unfished sites as well, but we are doing status of fished reefs)
sites_status<- everything$site_status
dim(sites_status)
ggplot()+geom_point(aes(x=sites_status[,1],y=sites_status[,20]))

#calculate median and 90% uncertainity intervals for each reef site
ordereddata<- rbind(reserves_complete,remote_complete,fished_data)
ordereddata<- droplevels(ordereddata)
colnames(sites_status)<- paste(ordereddata$UniqueSite,"status",sep="_")
ordereddata$site_status<- matrixStats::colMedians(sites_status)
ordereddata$sd_logsitesstatus<- matrixStats::colSds(log(sites_status))
ordereddata$mean_logsitesstatus<- matrixStats::colMeans2(log(sites_status))
ordereddata$site_status_op<- matrixStats::colQuantiles(sites_status,prob=0.9)
ordereddata$site_status_prec<- matrixStats::colQuantiles(sites_status,prob=0.1)

#marginalized biomass to account for sampling
sites_Bmarg<- everything$site_Bmarg
ordereddata$site_Bmarg<- matrixStats::colMedians(sites_Bmarg)

#marginalized biomass for sampling and environmental variables (for ecosystem metric section)
sites_Bmarg_all<- everything$site_Bmarg_all
ordereddata$site_Bmarg_all<- matrixStats::colMedians(sites_Bmarg_all)

#B0
sites_B0<- everything$site_B0
ordereddata$site_B0<- matrixStats::colMedians(sites_B0)
summary(ordereddata$site_B0)

#BMMSY
sites_BMMSY<- everything$site_BMMSY
sites_BMMSY2<- melt(sites_BMMSY)
ordereddata$site_BMMSY<- matrixStats::colMedians(sites_BMMSY)
ordereddata$site_BMMSY_op<- matrixStats::colQuantiles(sites_BMMSY,probs = 0.9)
ordereddata$site_BMMSY_prec<- matrixStats::colQuantiles(sites_BMMSY,probs = 0.1)

bmmsy_entiredist3<- ggplot()+geom_density(aes(x=sites_BMMSY2$value),fill="grey",alpha=0.5)+xlab("BMMSY (t/km2)")+ labs(x = expression ("B"["MMSY"]* " ("~t/km^2*")"))+xlim(c(0,350))+geom_rug(data=ordereddata, aes(x=site_BMMSY),col="gray41")+theme(axis.title.y = element_blank(),panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept =median(list_of_draws_full$BMMSY),lty=2 )+geom_text(aes(x=median(list_of_draws_full$BMMSY)+45,y=0.0150),label=paste("~", round(median(list_of_draws_full$BMMSY),2)),col="darkgrey")
summary(ordereddata$site_BMMSY)

#MMSY
sites_MMSY<- everything$site_MMSY
sites_MMSY2<- melt(sites_MMSY)
ordereddata$site_MMSY<- matrixStats::colMedians(sites_MMSY)
ordereddata$site_MMSY_op<- matrixStats::colQuantiles(sites_MMSY,probs = 0.9)
ordereddata$site_MMSY_prec<- matrixStats::colQuantiles(sites_MMSY,probs = 0.1)
mmsy_entiredist3<- ggplot()+geom_density(aes(x=sites_MMSY2$value),fill="grey",alpha=0.5)+xlab("MMSY (t/km2/y)")+ labs(x = expression ("MMSY ("~t/km^2/y*")"))+xlim(c(0,40))+geom_rug(data=ordereddata, aes(x=site_MMSY),col="gray41")+ylab("")+theme(panel.background = element_rect(fill="white",colour="black"))+geom_vline(xintercept =median(list_of_draws_full$MMSY),lty=2 )+geom_text(aes(x=median(list_of_draws_full$MMSY)+4,y=0.13),label=paste("~", round(median(list_of_draws_full$MMSY),2)),col="darkgrey")
summary(ordereddata$site_MMSY)
fig1a<- ggarrange(bmmsy_entiredist3,mmsy_entiredist3,nrow=2,ncol=1,labels=c("a","b"))

#PGMY
sites_PGMY<- 0.8*sites_MMSY
ordereddata$site_PGMY<- matrixStats::colMedians(sites_PGMY)

#Biomass at PGMY
#ESTIMATE THE BIOMASS AT PGMY (solving the equation where site_surplus is the sites PGMY (we will get two b but we keep the high one):site_surplus[,i]=everything$r*sites_Bmarg[,i]*(1-(sites_Bmarg[,i]/sites_B0[,i])))
Biomass_increase<-seq(0,round(max(ordereddata$site_B0)),0.1)
site_sustyield<-matrix(NA,nrow(ordereddata),length(Biomass_increase))
for(i in 1:length(Biomass_increase)){
site_sustyield[,i]<-log(ordereddata$site_B0)*median(everything$r)*Biomass_increase[i] * (1-(log(Biomass_increase[i])/log(ordereddata$site_B0)))  
}
site_sustyield[site_sustyield<0]<-0
site_sustyield[is.na(site_sustyield)]<-0
colnames(site_sustyield)<-Biomass_increase
#round to match the rounded pgmy
site_sustyield<-round(site_sustyield,1)
#get biomass value that matches those PGMY
for(i in 1:nrow(ordereddata)){
  ordereddata$site_BPGMY_low[i]<-as.numeric(first(names(site_sustyield[i,which(site_sustyield[i,]==round(ordereddata$site_PGMY[i],1))])))
  ordereddata$site_BPGMY_high[i]<-as.numeric(last(names(site_sustyield[i,which(site_sustyield[i,]==round(ordereddata$site_PGMY[i],1))])))
}
#check they are not the same
ggplot(ordereddata,aes(x=site_BPGMY_low,y=site_BPGMY_high))+geom_point()+geom_abline()


#combine for figure 1
#backtransform
coraltrends$HardCoral<-(((coraltrends$sHardCoral)*2*sd(sqrt(reefscale_data$HardCoral),na.rm=T))+mean(sqrt(reefscale_data$HardCoral),na.rm=T))^2
prodtrends$Ocean_prod<-(((prodtrends$sOcean_prod)*2*sd(log(reefscale_data$Prod_mgCm2d_finer),na.rm=T))+mean(log(reefscale_data$Prod_mgCm2d_finer),na.rm=T))
ssttrends$SST<-(((ssttrends$sSST)*2*sd(reefscale_data$SST_dC ,na.rm=T))+mean(reefscale_data$SST_dC ,na.rm=T))

a2<- ggplot()+
  #geom_point(data=ordereddata,aes(y=site_BMMSY,x=HardCoral,col=as.factor(Atoll)),alpha=0.1)+
  geom_line(data=coraltrends,aes(x=HardCoral,y=median_BMMSY,col=as.factor(coraltrends$Atoll)),size=2)+
  geom_ribbon(data=coraltrends,aes(x=HardCoral,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(coraltrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+labs(y = expression ("BMMSY ("~t/km^2*")"))+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)+
  geom_text(aes(60,30),label="Non-atoll",col="magenta4",size=4)+
  geom_text(aes(25,200),label="Atoll",col="steelblue2",size=4)+
  theme(axis.text.x = element_text(color="white"),panel.background = element_rect(fill="white",colour="black"))+xlab("")
a3<- ggplot()+
  #geom_point(data=ordereddata,aes(y=site_MMSY,x=HardCoral,col=as.factor(Atoll)),alpha=0.2)+
  geom_line(data=coraltrends,aes(x=HardCoral,y=median_MMSY,col=as.factor(coraltrends$Atoll)),size=2)+
  geom_ribbon(data=coraltrends,aes(x=HardCoral,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(coraltrends$Atoll)),alpha=0.3)+theme(panel.background = element_rect(fill="white",colour="black"))+xlab("% Hard Coral")+labs(y = expression ("MMSY("~t/km^2/y*")"))+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)

b2<- ggplot()+
  #geom_point(data=ordereddata,aes(y=site_BMMSY,x=log(Prod_mgCm2d_finer),col=as.factor(Atoll)),alpha=0.1)+
  geom_line(data=prodtrends,aes(x=Ocean_prod,y=median_BMMSY,col=as.factor(prodtrends$Atoll)),size=2)+
  geom_ribbon(data=prodtrends,aes(x=Ocean_prod,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+ylab("")+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)+
  theme(axis.text.x = element_text(color="white"),panel.background = element_rect(fill="white",colour="black"))+xlab("")
b3<- ggplot()+
 # geom_point(data=ordereddata,aes(y=site_MMSY,x=log(Prod_mgCm2d_finer),col=as.factor(Atoll)),alpha=0.2)+
  geom_line(data=prodtrends,aes(x=Ocean_prod,y=median_MMSY,col=as.factor(prodtrends$Atoll)),size=2)+
  geom_ribbon(data=prodtrends,aes(x=Ocean_prod,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(prodtrends$Atoll)),alpha=0.3)+theme(panel.background = element_rect(fill="white",colour="black"))+ylab("")+xlab("Productivity (log-mg/cm2/d)")+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)

c2<- ggplot()+
  #geom_point(data=ordereddata,aes(y=site_BMMSY,x=SST_dC,col=as.factor(Atoll)),alpha=0.1)+
  geom_line(data=ssttrends,aes(x=SST,y=median_BMMSY,col=as.factor(ssttrends$Atoll)),size=2)+
  geom_ribbon(data=ssttrends,aes(x=SST,ymax=conf.high_BMMSY,ymin=conf.low_BMMSY,fill=as.factor(ssttrends$Atoll)),alpha=0.2)+theme_classic()+xlab("")+ylab("")+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)+
  theme(axis.text.x = element_text(color="white"),panel.background = element_rect(fill="white",colour="black"))+xlab("")
c3<- ggplot()+
  #geom_point(data=ordereddata,aes(y=site_MMSY,x=SST_dC,col=as.factor(Atoll)),alpha=0.2)+
  geom_line(data=ssttrends,aes(x=SST,y=median_MMSY,col=as.factor(ssttrends$Atoll)),size=2)+
  geom_ribbon(data=ssttrends,aes(x=SST,ymax=conf.high_MMSY,ymin=conf.low_MMSY,fill=as.factor(ssttrends$Atoll)),alpha=0.3)+theme(panel.background = element_rect(fill="white",colour="black"))+xlab("SST (C)")+ylab("")+
  scale_fill_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+
  scale_color_manual(values = c("1" = "steelblue2",  "0"="magenta4"))+guides(fill=F,col=F)
windows()
fig1b<- ggarrange(a2,b2,c2,a3,b3,c3,nrow=2,ncol=3, heights = c(1,1),labels=c("c","d","e","f","g","h"))
fig1b

#sites actual  surplus based on their estimated standing stock  biomass
sites_surplus<-log(sites_B0)*list_of_draws_full$r*sites_Bmarg * (1-(log(sites_Bmarg)/log(sites_B0)))
sites_surplus[sites_surplus<0] <- 0

ordereddata$site_surplus<- matrixStats::colMedians(sites_surplus)
ordereddata$site_surplus_op<- matrixStats::colQuantiles(sites_surplus,prob=0.9)
ordereddata$site_surplus_prec<- matrixStats::colQuantiles(sites_surplus,prob=0.1)

#sites surplus relative to MMSY 
ordereddata$site_surplusrelMMSY<- ordereddata$site_surplus/ordereddata$site_MMSY
ordereddata$site_surplusrelMMSY_op<- ordereddata$site_surplus_prec/ordereddata$site_MMSY_prec
ordereddata$site_surplusrelMMSY_prec<-ordereddata$site_surplus_op/ordereddata$site_MMSY_op


#ADDED AFTER REVISIONS: site's individual catch status based on estimated cell mean spatial catch per km2 (reported and reconstructed)
ordereddata$site_catchstatus_reported<-ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_surplus
ordereddata$site_catchstatus_reported_op<- ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_surplus_op
ordereddata$site_catchstatus_reported_prec<- ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_surplus_prec

ordereddata$site_catchstatus_reportedMMSY<-ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_MMSY
ordereddata$site_catchstatus_reportedMMSY_op<- ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_MMSY_op
ordereddata$site_catchstatus_reportedMMSY_prec<- ordereddata$mean_tonnes_reported_km2cell/ordereddata$site_MMSY_prec

ordereddata$site_catchstatus<-ordereddata$mean_tonnes_km2cell/ordereddata$site_surplus
ordereddata$site_catchstatus_op<- ordereddata$mean_tonnes_km2cell/ordereddata$site_surplus_op
ordereddata$site_catchstatus_prec<- ordereddata$mean_tonnes_km2cell/ordereddata$site_surplus_prec

ordereddata$site_catchstatusMMSY<-ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY
ordereddata$site_catchstatusMMSY_op<- ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY_op
ordereddata$site_catchstatusMMSY_prec<- ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY_prec


#note that NAS are returned when surplus and catch is 0, inf is returned when surplus is 0
#for that shake, we assign 0 when the dividion returns NA and assume 0.0001 surplus when surplus is 0

ordereddata$site_catchstatus<-ifelse(is.na(ordereddata$site_catchstatus),0,
                                     ifelse(ordereddata$site_catchstatus==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus))
ordereddata$site_catchstatus_op<-ifelse(is.na(ordereddata$site_catchstatus_op),0,
                                     ifelse(ordereddata$site_catchstatus_op==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus_op))

ordereddata$site_catchstatus_prec<-ifelse(is.na(ordereddata$site_catchstatus_prec),0,
                                     ifelse(ordereddata$site_catchstatus_prec==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus_prec))

ordereddata$site_catchstatus_reported<-ifelse(is.na(ordereddata$site_catchstatus_reported),0,
                                     ifelse(ordereddata$site_catchstatus_reported==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus_reported))
ordereddata$site_catchstatus_reported_op<-ifelse(is.na(ordereddata$site_catchstatus_reported_op),0,
                                        ifelse(ordereddata$site_catchstatus_reported_op==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus_reported_op))

ordereddata$site_catchstatus_reported_prec<-ifelse(is.na(ordereddata$site_catchstatus_reported_prec),0,
                                          ifelse(ordereddata$site_catchstatus_reported_prec==Inf,ordereddata$mean_tonnes_km2cell/0.0001,ordereddata$site_catchstatus_reported_prec))


#probability of being beloW BMMSY for each site
#for each reef, calculate the probability of it being belowBMMSY based on their distributions
for (i in 1:length(ordereddata$site_status)){
  ordereddata$prob_belowBMMSY[i]=length(sites_status[,i][sites_status[,i]<1])/length(sites_status[,i])
}

#status fig:sites open to extraction (excluding those in high compliance reserves and remote)
extracteddata<-ordereddata[ordereddata$model_component==3,]
extracteddata<-droplevels(extracteddata)
ord <- with(extracteddata, reorder(definedprotection, site_status, median, order = TRUE, na.rm=T))
extracteddata<- within(extracteddata, 
                       definedprotection <- factor(definedprotection, 
                                                   levels=levels(ord)))



#categorize
extracteddata$belowBMMSY<-ifelse(extracteddata$site_status<1,"belowBMMSY","Not belowBMMSY")
extracteddata$aboveMMSY<-ifelse(extracteddata$site_catchstatusMMSY>1,"aboveMMSY","Not aboveBMMSY")
extracteddata$overfishing_site<-ifelse(extracteddata$site_catchstatus>1,"abovesurplus","Not abovesurplus")

#figure
#add cut-off (locations predicted to be above their median unfished biomass, make them be at their unfished biomass)
#likewise, add cut-off for places ctahing above MMSY
extracteddata$site_status_cutoff<-ifelse(extracteddata$site_status>2,2,extracteddata$site_status)
extracteddata$site_catchstatusMMSY_cutoff<-ifelse(extracteddata$site_catchstatusMMSY>2,2,extracteddata$site_catchstatusMMSY)

#color also if collapsed
extracteddata$site_status_color<-ifelse(extracteddata$site_status>1,"navyblue",ifelse(extracteddata$site_status<0.2,"black","red"))
extracteddata$site_catchstatusMMSY_color<-ifelse(extracteddata$site_status>1,"red","navyblue")

#perc of sites below BMMSY
(length(extracteddata$site_status[extracteddata$site_status<1])/length(extracteddata$site_status))*100
(length(extracteddata$site_status[extracteddata$site_status_op<1])/length(extracteddata$site_status))*100
(length(extracteddata$site_status[extracteddata$site_status_prec<1])/length(extracteddata$site_status))*100

#percent of sites above MMSY
(length(extracteddata$site_catchstatusMMSY[extracteddata$site_catchstatusMMSY>1])/length(extracteddata$site_catchstatusMMSY))*100
(length(extracteddata$site_catchstatusMMSY[extracteddata$site_catchstatusMMSY_op>1])/length(extracteddata$site_catchstatusMMSY))*100
(length(extracteddata$site_catchstatusMMSY[extracteddata$site_catchstatusMMSY_prec>1])/length(extracteddata$site_catchstatusMMSY))*100

#perc of sites below 0.1of B0
(length(extracteddata$site_status[extracteddata$site_status<0.2])/length(extracteddata$site_status))*100
(length(extracteddata$site_status[extracteddata$site_status_op<0.2])/length(extracteddata$site_status))*100
(length(extracteddata$site_status[extracteddata$site_status_prec<0.2])/length(extracteddata$site_status))*100

#percent of sites overfishing
(length(extracteddata$site_catchstatus[extracteddata$site_catchstatus>1])/length(extracteddata$site_catchstatus))*100
(length(extracteddata$site_catchstatus[extracteddata$site_catchstatus_op>1])/length(extracteddata$site_catchstatus))*100
(length(extracteddata$site_catchstatus[extracteddata$site_catchstatus_prec>1])/length(extracteddata$site_catchstatus))*100

#percent of sites below BPGMSY LOW
extracteddata$site_status_BPGMY_low<-extracteddata$site_Bmarg/extracteddata$site_BPGMY_low
(length(extracteddata$site_status_BPGMY_low[extracteddata$site_status_BPGMY_low<1])/length(extracteddata$site_status_BPGMY_low))*100

#percentage of fished reefs at PGMSY (surplus is within 80% of MMSY)
(nrow(extracteddata[!extracteddata$site_surplusrelMMSY<0.8,])/nrow(extracteddata))*100

#sustainability status 
extracteddata$fisherystatus<- as.factor(ifelse(extracteddata$site_catchstatus<1 & extracteddata$site_status>1, "In good condition", ifelse(extracteddata$site_catchstatus>1 & extracteddata$site_status<1, "Unsustainable", ifelse(extracteddata$site_catchstatus>1 & extracteddata$site_status>1, "Warning", "Recovering"))))
extracteddata$fisherystatus_prec<- as.factor(ifelse(extracteddata$site_catchstatus_prec<1 & extracteddata$site_status_prec>1, "In good condition", ifelse(extracteddata$site_catchstatus_prec>1 & extracteddata$site_status_prec<1, "Unsustainable", ifelse(extracteddata$site_catchstatus_prec>1 & extracteddata$site_status_prec>1, "Warning", "Recovering"))))
extracteddata$fisherystatus_op<- as.factor(ifelse(extracteddata$site_catchstatus_op<1 & extracteddata$site_status_op>1, "In good condition", ifelse(extracteddata$site_catchstatus_op>1 & extracteddata$site_status_op<1, "Unsustainable", ifelse(extracteddata$site_catchstatus_op>1 & extracteddata$site_status_op>1, "Warning", "Recovering"))))

#(incorporating quasi-sustainable) as sustainable
extracteddata$fisherystatus<- ifelse(extracteddata$fisherystatus=="Warning" &extracteddata$site_catchstatusMMSY<1, "In good condition",as.character(extracteddata$fisherystatus))
extracteddata$fisherystatus_prec<- ifelse(extracteddata$fisherystatus_prec=="Warning" &extracteddata$site_catchstatusMMSY_prec<1, "In good condition",as.character(extracteddata$fisherystatus_prec))
extracteddata$fisherystatus_op<- ifelse(extracteddata$fisherystatus_op=="Warning" &extracteddata$site_catchstatusMMSY_op<1, "In good condition",as.character(extracteddata$fisherystatus_op))

summary(as.factor(extracteddata$fisherystatus))
paste(c(round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Unsustainable"])/length(extracteddata$fisherystatus))*100),
  round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Warning"])/length(extracteddata$fisherystatus))*100),
    round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Recovering"])/length(extracteddata$fisherystatus))*100),
      round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="In good condition"])/length(extracteddata$fisherystatus))*100)), "%", sep=" ")

paste(c(round((length(extracteddata$fisherystatus_prec[extracteddata$fisherystatus_prec=="Unsustainable"])/length(extracteddata$fisherystatus_prec))*100),
        round((length(extracteddata$fisherystatus_prec[extracteddata$fisherystatus_prec=="Warning"])/length(extracteddata$fisherystatus_prec))*100),
        round((length(extracteddata$fisherystatus_prec[extracteddata$fisherystatus_prec=="Recovering"])/length(extracteddata$fisherystatus_prec))*100),
        round((length(extracteddata$fisherystatus_prec[extracteddata$fisherystatus_prec=="In good condition"])/length(extracteddata$fisherystatus_prec))*100)), "%", sep=" ")

paste(c(round((length(extracteddata$fisherystatus_op[extracteddata$fisherystatus_op=="Unsustainable"])/length(extracteddata$fisherystatus_op))*100),
        round((length(extracteddata$fisherystatus_op[extracteddata$fisherystatus_op=="Warning"])/length(extracteddata$fisherystatus_op))*100),
        round((length(extracteddata$fisherystatus_op[extracteddata$fisherystatus_op=="Recovering"])/length(extracteddata$fisherystatus_op))*100),
        round((length(extracteddata$fisherystatus_op[extracteddata$fisherystatus_op=="In good condition"])/length(extracteddata$fisherystatus_op))*100)), "%", sep=" ")

# stats using reported data (e.g., conservation concern)
length(extracteddata$site_status[extracteddata$site_catchstatus_reportedMMSY>1])/length(extracteddata$site_status)
length(extracteddata$site_status[extracteddata$site_catchstatus_reported>1])/length(extracteddata$site_status)
length(extracteddata$site_status[extracteddata$site_status<1|extracteddata$site_catchstatus_reportedMMSY>1])/length(extracteddata$site_status)

length(extracteddata$site_status[extracteddata$site_catchstatusMMSY>1])/length(extracteddata$site_status)
length(extracteddata$site_status[extracteddata$site_catchstatus>1])/length(extracteddata$site_status)
length(extracteddata$site_status[extracteddata$site_status<1|extracteddata$site_catchstatusMMSY>1])/length(extracteddata$site_status)

c<-ggplot() + 
  geom_bar(data=extracteddata, aes(fill=factor(fisherystatus, levels=c("In good condition","Recovering","Warning","Unsustainable")), x=""),position="fill")+
  scale_fill_manual (name="",values=c( "Recovering"="cyan3",
                                                   "In good condition"="navyblue",
                                                   "Unsustainable"="red", 
                                                   "Warning"="goldenrod1"),
                   labels = c("Recovering"="Potentially
recovering",
                              "In good condition"="In good condition",
                              "Unsustainable"="Unsustainable", 
                              "Warning"="Warning"))+theme(axis.title.x = element_blank(),panel.background = element_rect(fill="white",colour="black"))+scale_y_continuous(breaks=c(0.00,0.25,0.50,0.75,1.00), limits=c(0,1),labels=c("0","25","50","75","100"))+labs(y = expression ("Fishery status (% of sites)"))+guides(fill=F)+
  geom_text(aes(x="",y=c(0.17,0.45,0.67,0.89)),label=paste(c(round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Unsustainable"])/length(extracteddata$fisherystatus))*100),
                                                      round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Warning"])/length(extracteddata$fisherystatus))*100),
                                                      round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="Recovering"])/length(extracteddata$fisherystatus))*100),
                                                      round((length(extracteddata$fisherystatus[extracteddata$fisherystatus=="In good condition"])/length(extracteddata$fisherystatus))*100)), "%", sep=" "),col="white")
  

#figure of status
a<-ggplot()+geom_jitter(data=extracteddata,aes(y=site_status_cutoff,x="",fill=as.factor(extracteddata$site_status_cutoff<1),color=as.factor(extracteddata$site_status_cutoff<1)),pch=21,alpha=0.2)+geom_violin(data=extracteddata,aes(y=site_status_cutoff,x="" ),scale="width",fill="gray",col="darkgrey",draw_quantiles = c( 0.5),alpha=0.5,lwd=1.2)+ 
  geom_text(aes(x="",y=-0.1),col="darkred",lwd=5,label=paste(round((length(extracteddata$site_status[extracteddata$site_status<1])/length(extracteddata$site_status))*100),"%"))+
  
  geom_hline(aes(yintercept=1),lty=2,lwd=1.1)+xlab("")+
  scale_fill_manual(values = c("FALSE" = "navyblue",
                               "TRUE" = "red"),
                    labels = c("FALSE" ="Not belowBMMSY", 
                               "TRUE" ="belowBMMSY"))+
  scale_color_manual(values = c("FALSE" = "navyblue",
                                "TRUE" = "red"),
                     labels = c("FALSE" ="Not belowBMMSY", 
                                "TRUE" ="belowBMMSY"))+guides(col=F,fill=F)+theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(color="black"),axis.title = element_text(family="Helvetica",size=12),panel.background = element_rect(fill="white",colour="black"))+scale_y_continuous(breaks=c(0,1,2), limits=c(-0.1,2.1),labels=c("0","1",">2"))+labs(y = expression ("Biomass status (B/ B"["MMSY"]* ")"))


b<-ggplot()+geom_jitter(data=extracteddata,aes(y=site_catchstatusMMSY_cutoff,x="",fill=as.factor(extracteddata$site_catchstatusMMSY_cutoff>1),color=as.factor(extracteddata$site_catchstatusMMSY_cutoff>1)),pch=21,alpha=0.2)+geom_violin(data=extracteddata,aes(y=site_catchstatusMMSY_cutoff,x="" ),scale="width",fill="gray",col="darkgrey",draw_quantiles = c( 0.5),alpha=0.5,lwd=1.2)+ 
  geom_text(aes(x="",y=2.15),col="darkred",lwd=5,label=paste(round((length(extracteddata$site_catchstatusMMSY[extracteddata$site_catchstatusMMSY>1])/length(extracteddata$site_catchstatusMMSY))*100),"%"))+
  
  geom_hline(aes(yintercept=1),lty=2,lwd=1.1)+xlab("")+
  scale_fill_manual(values = c("FALSE" = "navyblue",
                               "TRUE" = "red"),
                    labels = c("FALSE" ="Not belowBMMSY", 
                               "TRUE" ="belowBMMSY"))+
  scale_color_manual(values = c("FALSE" = "navyblue",
                                "TRUE" = "red"),
                     labels = c("FALSE" ="belowMMSY", 
                                "TRUE" ="aboveMMSY"))+guides(col=F,fill=F)+theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x = element_text(color="black"),axis.title = element_text(family="Helvetica",size=12),panel.background = element_rect(fill="white",colour="black"))+scale_y_continuous(breaks=c(0,1,2), limits=c(-0.1,2.15),labels=c("0","1",">2"))+labs(y = expression ("Fishing status (C/MMSY)"))

fig1c<-ggarrange(a,b,nrow=2,ncol=1,labels=c("i","j"),heights=c(1,1))
fig1c<-ggarrange(fig1c,c,nrow=1,ncol=2,labels=c("","k"),widths=c(1.5,1))


windows()
fig1<-ggarrange(fig1a,fig1b,fig1c,nrow=1,ncol=3,widths=c(1,2,1.1))
annotate_figure(fig1,left="Posterior density")


#Estimates by jurisdiction
#jurisidction number of sites
country_sites<- ddply(ordereddata,.(Larger),summarize, overall_count=length(site_MMSY))
country_fishedsites<- ddply(fished_data,.(Larger),summarize, fished_count=length(FamBiomass_tkm2))
country_sites<- merge(country_sites, country_fishedsites,by="Larger",all.x=T)
country_fsites<- country_sites[!is.na(country_sites$fished_count),]


#B0/BMMSY/MMSY distribution per jurisdiction: calculating the mean of samples(ITERATIONS) for sites within jurisdictions
sites_B0_juris<-sites_B0
colnames(sites_B0_juris)<-ordereddata$Larger
sites_B0_juris<-melt(sites_B0_juris)
sites_B0_juris$iter_country<-paste(sites_B0_juris$iterations,sites_B0_juris$Var.2,sep="::")
sites_B0_juris$model_component<-rep(ordereddata$model_component,each=4000)

juris_B0distributions<-ddply(sites_B0_juris,.(iter_country),summarise,country=Var.2[1],iter=iterations[1],meanlogB0=mean(log(value)),n_sites=length(Var.2))
juris_sites<-ddply(sites_B0_juris,.(Var.2),summarise,n_sites=length(Var.2)/4000)

fac2 <- with( juris_B0distributions, reorder(country, meanlogB0, median, order = TRUE, na.rm=T))
juris_B0distributions<- within(juris_B0distributions, 
                               country<- factor(country, 
                                                levels=levels(fac2)))

#jurisdiction B0,BMMSY and MMSY
juris_B0distributions$B0<-exp(juris_B0distributions$meanlogB0)
summary(juris_B0distributions$B0)
juris_B0distributions$meanBMMSY<-juris_B0distributions$B0/2.718281828
juris_B0distributions_h<-reshape2::dcast(juris_B0distributions, iter ~ country, value.var="B0")
juris_MMSYdistributions<-melt((list_of_draws_full$r*juris_B0distributions_h[,-1])/2.718281828)
juris_MMSYdistributions$iter<-rep(seq(1:4000),length(unique(juris_MMSYdistributions$variable)))
juris_MMSYdistributions$iter_country<-paste(juris_MMSYdistributions$iter,juris_MMSYdistributions$variable,sep="::")
juris_MMSYdistributions$meanMMSY<-juris_MMSYdistributions$value

#Jurisdiction SURPLUS
Bio2<-seq(0,600,1)
popgrowth_full_juris<- array(NA, c(length(Bio2),dim(juris_B0distributions_h[,-1])))
for (i in 1:length(Bio2)){
  popgrowth_full_juris[i,,]<- as.matrix(log(juris_B0distributions_h[,-1])*list_of_draws_full$r*Bio2[i]  * (1-(log(Bio2[i] )/log(juris_B0distributions_h[,-1]))))  
} 
popgrowth_full_juris[1,,]
popgrowth_full_juris[popgrowth_full_juris<0]<-0
popgrowth_full_juris[is.nan(popgrowth_full_juris)]<-0

#calculate median and quantile surplus per jurisdiction for each biomass value
mediansurplus_perjuris<-as.data.frame(apply(popgrowth_full_juris, c(1,3), median))
colnames(mediansurplus_perjuris)<-colnames(juris_B0distributions_h[,-1])
precsurplus_perjuris<-as.data.frame(apply(popgrowth_full_juris, c(1,3), quantile,prob=0.1))
colnames(precsurplus_perjuris)<-colnames(juris_B0distributions_h[,-1])
opsurplus_perjuris<-as.data.frame(apply(popgrowth_full_juris, c(1,3), quantile,prob=0.9))
colnames(opsurplus_perjuris)<-colnames(juris_B0distributions_h[,-1])

mediansurplus_perjuris2<- melt(mediansurplus_perjuris)
precsurplus_perjuris2<- melt(precsurplus_perjuris)
opsurplus_perjuris2<- melt(opsurplus_perjuris)
mediansurplus_perjuris2$Biomass<- rep(Bio2,length.out=length(Bio2)*ncol(mediansurplus_perjuris))
precsurplus_perjuris2$Biomass<- rep(Bio2,length.out=length(Bio2)*ncol(precsurplus_perjuris))
opsurplus_perjuris2$Biomass<- rep(Bio2,length.out=length(Bio2)*ncol(opsurplus_perjuris))

#distributions per jurisdiction
BMMSY_perjuris_fig<-ggplot()+geom_density_ridges(data=juris_B0distributions,aes(x=meanBMMSY,y=country),alpha=0.5,rel_min_height = 0.01)+ylab("")+xlab("Jurisdiction-level BMMSY (t/km2)")+geom_vline(xintercept = median(list_of_draws_full$BMMSY),lty=2)+xlim(c(0,300))+theme(panel.background = element_rect(fill="white",color="black"))+ labs(x =  expression ("Jurisdiction-level "*B["MMSY "]*"("~t/km^2*")"))
MMSY_perjuris_fig<-ggplot()+geom_density_ridges(data=juris_MMSYdistributions,aes(x=value,y=variable),alpha=0.5,rel_min_height = 0.01)+ylab("")+xlab("Jurisdiction-level MMSY (t/km2/y)")+geom_vline(xintercept = median(list_of_draws_full$MMSY),lty=2)+xlim(c(0,40))+theme(axis.text.y =element_blank(),panel.background = element_rect(fill="white",color="black"))+ labs(x =  expression ("Jurisdiction-level MMSY"*"("~t/km^2/y*")"))

#country-specific surplus 
countryspcfig<-ggplot(NULL)+geom_ribbon(aes(x=precsurplus_perjuris2$Biomass,ymin=precsurplus_perjuris2$value,ymax=opsurplus_perjuris2$value, group=precsurplus_perjuris2$variable), fill="darkgrey", alpha=0.1)+
  geom_line(data=mediansurplus_perjuris2,aes(y=value, x=Biomass,group=variable), col="grey25",alpha=0.7,lty=3)+ylim(c(0,41))+
  geom_line(data=popgrowth_full_open2 ,aes(y=median, x=Biomass), col="black",lwd=1.5,lty=1)+xlim(c(0,600))+
   theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))
windows()
ggarrange(BMMSY_perjuris_fig,MMSY_perjuris_fig,countryspcfig,nrow=1,ncol=3,labels = c("a","b","c"),widths=c(1.4,1,1))

#get biomass values of extracted reefs per jurisdiction. Assuming protected areas are at B0 values
sites_B_juris<-sites_Bmarg
colnames(sites_B_juris)<-ordereddata$Larger
sites_B_juris<-melt(sites_B_juris)
colnames(sites_B_juris)<-c("iterations","jurisdiction","Bmarg")
sites_B_juris$iter_country<-paste(sites_B_juris$iterations,sites_B_juris$jurisdiction,sep="::")
sites_B_juris$model_component<-rep(ordereddata$model_component,each=4000)
sites_B_juris$B0<-sites_B0_juris$value
sites_B_juris$definedprotection<-rep(ordereddata$definedprotection ,each=4000)

#only extracted reefs
sites_B_juris_extracted<-sites_B_juris[sites_B_juris$model_component==3,]
#get mean of sites B0 and r for each sample per jurisdiction
juris_Bdistributions<-ddply(sites_B_juris_extracted,.(iter_country),summarise,jurisdiction=jurisdiction[1],iterations=iterations[1],meanBmarg=exp(mean(log(Bmarg))))
length(unique(juris_Bdistributions$iter_country))
length(unique(juris_Bdistributions$jurisdiction))
length(unique(juris_MMSYdistributions$iter_country))
length(unique(juris_MMSYdistributions$variable))
dim(juris_Bdistributions)

#add proportion of waters protected
juris_Bdistributions<- merge(juris_Bdistributions,jurisdictionscale_data2[,c("Larger","mpa_perc")],by.x="jurisdiction",by.y="Larger",all.x=T)

#calculate weighted biomass
length(unique(juris_B0distributions$country))

#merge with country MMSY and BMMSY distributions
juris_Bdistributions<-merge(juris_Bdistributions,juris_B0distributions[c("iter_country","meanBMMSY","B0")],by="iter_country",all.x=T)
juris_Bdistributions<-merge(juris_Bdistributions,juris_MMSYdistributions[,c("iter_country","meanMMSY")],by="iter_country",all.x=T)

#calculate status (using only extracted reefs , and also doing a weighted average assuming the reported waters protected are at B0)
juris_Bdistributions$weightedB_B0<-juris_Bdistributions$meanBmarg*(1-(juris_Bdistributions$mpa_perc/100))+juris_Bdistributions$B0*(juris_Bdistributions$mpa_perc/100)
juris_Bdistributions$Bstatus_onlyfishedreefs<-juris_Bdistributions$meanBmarg/juris_Bdistributions$meanBMMSY
juris_Bdistributions$Bstatus_B0<-juris_Bdistributions$weightedB_B0/juris_Bdistributions$meanBMMSY
summary(juris_Bdistributions$Bstatus_B0)
summary(juris_Bdistributions$Bstatus_onlyfishedreefs)

#Recall this is examining the status of fished reefs, PRIA AND Chagos did have reefs open to fishing at the time of sampling
ordereddata$Site[ordereddata$Larger=="British Indian Ocean Territory"& (ordereddata$Management=="Fished"|ordereddata$Management=="Restricted")]
ordereddata$Site[ordereddata$Larger=="PRIA"& (ordereddata$Management=="Fished"|ordereddata$Management=="Restricted")]

#figure of biomass status
fac2 <- with( juris_Bdistributions, reorder(jurisdiction, Bstatus_onlyfishedreefs, median, order = TRUE, na.rm=T))
juris_Bdistributions_ord<- within(juris_Bdistributions, 
                               jurisdiction<- factor(jurisdiction, 
                                                levels=levels(fac2)))

prob_belowBMMSY_fig<-ggplot(NULL)+geom_density_ridges_gradient(data=juris_Bdistributions_ord,aes(x=log2(Bstatus_onlyfishedreefs), y=jurisdiction, fill=stat(x)),scale=2,rel_min_height = 0.01)+geom_vline(xintercept = log2(1), lty=2,lwd=1.3)+
  theme(axis.text.y = element_text(size=6))+guides(fill=F)+
  scale_fill_gradient2("log2(B/BMMSY)",high="navyblue",mid="white",low="darkred",midpoint = log2(1))+xlab("log2(B/BMMSY)")+ylab("")+theme(axis.text.y = element_text(size=7),panel.background = element_rect(fill="white",color="black"))+xlim(c(-5,5))

prob_belowBMMSY_figweighted<-ggplot(NULL)+geom_density_ridges_gradient(data=juris_Bdistributions_ord,aes(x=log2(Bstatus_B0), y=jurisdiction, fill=stat(x)),scale=2,rel_min_height = 0.01)+geom_vline(xintercept = log2(1), lty=2,lwd=1.3)+
  theme(axis.text.y = element_text(size=6))+guides(fill=F)+
  scale_fill_gradient2("log2(B/BMMSY)",high="navyblue",mid="white",low="darkred",midpoint = log2(1))+xlab("log2(Bweighted/BMMSY)")+ylab("")+theme(axis.text.y = element_blank(),panel.background = element_rect(fill="white",color="black"))+xlim(c(-5,5))
windows()
prob_belowBMMSY_fig_all<-ggarrange(prob_belowBMMSY_fig,prob_belowBMMSY_figweighted, widths = c(1.8,1),labels=c("b","c"))

#country ref points/status
country_refpoints<- ddply(juris_Bdistributions,.(jurisdiction),summarize, larger_MMSY=median(meanMMSY,na.rm=T),larger_B0=median(B0,na.rm=T),larger_BMMSY=median(meanBMMSY,na.rm=T),larger_MMSY_low=quantile(meanMMSY, probs=c(0.1),na.rm=T),larger_MMSY_high=quantile(meanMMSY, probs=c(0.9),na.rm=T),larger_BMMSY_low=quantile(meanBMMSY, probs=c(0.1),na.rm=T),larger_BMMSY_high=quantile(meanBMMSY, probs=c(0.9),na.rm=T),larger_B0_low=quantile(B0, probs=c(0.1),na.rm=T),larger_B0_high=quantile(B0, probs=c(0.9),na.rm=T), larger_Bstatus=median(Bstatus_onlyfishedreefs),larger_Bstatus_B0=median(Bstatus_B0,na.rm=T), larger_Bstatus_low=quantile(Bstatus_onlyfishedreefs, probs=c(0.1),na.rm=T),larger_Bstatus_B0_low=quantile(Bstatus_B0, probs=c(0.1),na.rm=T), larger_Bstatus_high=quantile(Bstatus_onlyfishedreefs, probs=c(0.9),na.rm=T),larger_Bstatus_B0_high=quantile(Bstatus_B0, probs=c(0.9),na.rm=T), sd_logBstatus=sd(log(Bstatus_onlyfishedreefs),na.rm=T), sd_logBstatus_B0=sd(log(Bstatus_B0),na.rm=T), mean_logBstatus=mean(log(Bstatus_onlyfishedreefs),na.rm=T), mean_logBstatus_B0=mean(log(Bstatus_B0),na.rm=T),larger_weightedB=median(weightedB_B0,na.rm=T),larger_weightedB_low=quantile(weightedB_B0,prob=0.1,na.rm=T),larger_weightedB_high=quantile(weightedB_B0,prob=0.9,na.rm=T), larger_Bmarg=median(meanBmarg,na.rm=T),larger_Bmarg_low=quantile(meanBmarg,prob=0.1,na.rm=T),larger_Bmarg_high=quantile(meanBmarg,prob=0.9,na.rm=T))
country_refpoints$colorBstatus<-ifelse(country_refpoints$larger_Bstatus<1 &country_refpoints$larger_Bstatus_B0<1,"always_below",ifelse(country_refpoints$larger_Bstatus>1 &country_refpoints$larger_Bstatus_B0>1,"always_above","not_consistent"))
ggplot(country_refpoints,aes(x=log(larger_Bstatus),y=log(larger_Bstatus_B0)))+
  geom_pointrange(aes(ymax=log(larger_Bstatus_B0_high),ymin=log(larger_Bstatus_B0_low)),col="grey")+
  geom_errorbarh(aes(xmax=log(larger_Bstatus_high),xmin=log(larger_Bstatus_low)),col="grey")+geom_point(aes(col=colorBstatus))+geom_text_repel(aes(label=country_refpoints$jurisdiction),size=3)+geom_vline(xintercept = 0,lty=2)+geom_hline(yintercept = 0,lty=2)+theme_classic()+ylab("log(Bweighted/BMMSY)")+xlab("log(B/BMMSY)")
country_refpoints$jurisdiction[country_refpoints$colorBstatus=="not_consistent"]
#Jurisdiction Fishing status
#calculate catch per unit area of reef
jurisdictionscale_data2$catch_tkm2<- jurisdictionscale_data2$spatial_tonnes_km2 

# merge and country analyses
jurisdictionscale_data2<- merge(jurisdictionscale_data2,country_refpoints,by.x="Larger", by.y="jurisdiction",all.x=T )


#for jurisidctions without data we use the combined distribution of the jurisdictions
#add the median and quantiles for jurisdictions that dont have biomas data
jurisdictionscale_data2$MMSY_median<- ifelse(is.na(jurisdictionscale_data2$larger_MMSY),median(juris_Bdistributions$meanMMSY),jurisdictionscale_data2$larger_MMSY)
jurisdictionscale_data2$MMSY_low<- ifelse(is.na(jurisdictionscale_data2$larger_MMSY_low),quantile(juris_Bdistributions$meanMMSY,probs = c(0.1)),jurisdictionscale_data2$larger_MMSY_low)
jurisdictionscale_data2$MMSY_high<- ifelse(is.na(jurisdictionscale_data2$larger_MMSY_high),quantile(juris_Bdistributions$meanMMSY,probs = c(0.9)),jurisdictionscale_data2$larger_MMSY_high)
jurisdictionscale_data2$BMMSY_median<- ifelse(is.na(jurisdictionscale_data2$larger_BMMSY),median(juris_Bdistributions$meanBMMSY),jurisdictionscale_data2$larger_BMMSY)
jurisdictionscale_data2$BMMSY_low<- ifelse(is.na(jurisdictionscale_data2$larger_BMMSY_low),quantile(juris_Bdistributions$meanBMMSY,probs = c(0.1)),jurisdictionscale_data2$larger_BMMSY_low)
jurisdictionscale_data2$BMMSY_high<- ifelse(is.na(jurisdictionscale_data2$larger_BMMSY_high),quantile(juris_Bdistributions$meanBMMSY,probs = c(0.9)),jurisdictionscale_data2$larger_BMMSY_high)
jurisdictionscale_data2$B0_median<- ifelse(is.na(jurisdictionscale_data2$larger_B0),median(juris_Bdistributions$B0),jurisdictionscale_data2$larger_B0)
jurisdictionscale_data2$B0_low<- ifelse(is.na(jurisdictionscale_data2$larger_B0_low),quantile(juris_Bdistributions$B0,probs = c(0.1)),jurisdictionscale_data2$larger_B0_low)
jurisdictionscale_data2$B0_high<- ifelse(is.na(jurisdictionscale_data2$larger_B0_high),quantile(juris_Bdistributions$B0,probs = c(0.9)),jurisdictionscale_data2$larger_B0_high)


#plot combined distributions
combinedBMMSY<- ggplot()+geom_density(aes(x=juris_Bdistributions$meanBMMSY),fill="grey",alpha=0.5)+xlab("BMMSY (t/km2)")+ labs(x = expression ("B"["MMSY"]* " ("~t/km^2*")"))+xlim(c(0,300))+geom_rug(data=jurisdictionscale_data2, aes(x=BMMSY_median),col="gray41")+ylab(" ")+theme(panel.background = element_rect(fill="white",colour="black"))
#MMSY
combinedMMSY<-ggplot()+geom_density(aes(x=juris_Bdistributions$meanMMSY),fill="grey",alpha=0.5)+xlab("MMSY (t/km2/y)")+ labs(x = expression ("MMSY ("~t/km^2/y*")"))+xlim(c(0,35))+geom_rug(data=jurisdictionscale_data2, aes(x=MMSY_median),col="gray41")+ylab("")+theme(panel.background = element_rect(fill="white",colour="black"))

#fishing status for jurisdiction without biomass data
jurisdictions_nob<- jurisdictionscale_data2[is.na(jurisdictionscale_data2$larger_MMSY),]
jurisdictions_nob<- droplevels(jurisdictions_nob)
country_proboverfishing<- matrix(NA,nrow=length(juris_Bdistributions$meanMMSY),ncol=length(jurisdictions_nob$Larger))
for (i in 1:length(jurisdictions_nob$catch_tkm2)){
  country_proboverfishing[,i]<- jurisdictions_nob$catch_tkm2[i]/juris_Bdistributions$meanMMSY
}
colnames(country_proboverfishing)<-jurisdictions_nob$area_name
country_proboverfishing2<-melt(country_proboverfishing)
colnames(country_proboverfishing2)<-c("iterations","jurisdiction","fishingstatus")

#probability of overfishing for countries with biomass data 
juris_Bdistributions<-merge(juris_Bdistributions,jurisdictionscale_data2[,c("Larger","catch_tkm2")],by.x="jurisdiction",by.y="Larger",all.x=T)
juris_Bdistributions$fishingstatus<-juris_Bdistributions$catch_tkm2/juris_Bdistributions$meanMMSY
colnames(juris_Bdistributions)
#merge them
country_proboverfishing2<-rbind(country_proboverfishing2,juris_Bdistributions[,c("iterations","jurisdiction","fishingstatus")])
#figure of fishing status
fac2 <- with( country_proboverfishing2, reorder(jurisdiction, fishingstatus, median, order = TRUE, na.rm=T))
country_proboverfishing2_ord<- within(country_proboverfishing2, 
                                  jurisdiction<- factor(jurisdiction, 
                                                        levels=levels(fac2)))
summary(country_proboverfishing2_ord)
prob_aboveMMSY_fig<-ggplot(NULL)+geom_density_ridges_gradient(data=country_proboverfishing2_ord[!is.na(country_proboverfishing2_ord$fishingstatus),],aes(x=sqrt(sqrt(fishingstatus)), y=jurisdiction, fill=stat(x)),scale=2,rel_min_height = 0.01)+geom_vline(xintercept = sqrt(sqrt(1)), lty=2,lwd=1.3)+
  theme(axis.text.y = element_text(size=6))+guides(fill=F)+
  scale_fill_gradient2("fourth-root (C/MMSY)",low="navyblue",mid="white",high="darkred",midpoint = sqrt(sqrt(1)))+xlab("fourth-root (C/MMSY)")+ylab("")+theme(axis.text.y = element_text(size=7),panel.background = element_rect(fill="white",color="black"))

#site-sepcific distribution status for catch for sites with biomass data and the specific conditions MMSY for others
overfishing_country<- ddply(country_proboverfishing2_ord,.(jurisdiction),summarize, medianFstatus=median(fishingstatus,na.rm=T),Fstatus_low=quantile(fishingstatus,na.rm=T,prob=0.1),Fstatus_high=quantile(fishingstatus,na.rm=T,prob=0.9))
jurisdictionscale_data2$Larger2<-ifelse(is.na(jurisdictionscale_data2$mean_logBstatus),jurisdictionscale_data2$area_name,jurisdictionscale_data2$Larger)
jurisdictionscale_data2<-merge(jurisdictionscale_data2,overfishing_country,by.x="Larger2",by.y="jurisdiction",all.x=T)


#jurisdictions fishing status: 
jurisdictionscale_data2$Larger[is.na(jurisdictionscale_data2$medianFstatus)]
jurisdictionscale_data2$overfishing<- ifelse(is.na(jurisdictionscale_data2$medianFstatus),NA,ifelse(jurisdictionscale_data2$medianFstatus>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$overfishing)&jurisdictionscale_data2$overfishing=="Overfishing"])/length(jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$overfishing)])
jurisdictionscale_data2$overfishing_prec<- ifelse(is.na(jurisdictionscale_data2$Fstatus_low),NA,ifelse(jurisdictionscale_data2$Fstatus_low>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data2$overfishing_prec[!is.na(jurisdictionscale_data2$overfishing_prec)&jurisdictionscale_data2$overfishing_prec=="Overfishing"])/length(jurisdictionscale_data2$overfishing_prec[!is.na(jurisdictionscale_data2$overfishing_prec)])
jurisdictionscale_data2$overfishing_op<- ifelse(is.na(jurisdictionscale_data2$Fstatus_high),NA,ifelse(jurisdictionscale_data2$Fstatus_high>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data2$overfishing_op[!is.na(jurisdictionscale_data2$overfishing_op)&jurisdictionscale_data2$overfishing_op=="Overfishing"])/length(jurisdictionscale_data2$overfishing_op[!is.na(jurisdictionscale_data2$overfishing_op)])

#jurisdictions biomass status: 
jurisdictionscale_data2$belowBMMSY<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus),NA,ifelse(jurisdictionscale_data2$larger_Bstatus>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY[!is.na(jurisdictionscale_data2$belowBMMSY)&jurisdictionscale_data2$belowBMMSY=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY[!is.na(jurisdictionscale_data2$belowBMMSY)])
jurisdictionscale_data2$belowBMMSY_prec<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_low),NA,ifelse(jurisdictionscale_data2$larger_Bstatus_low>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY_prec[!is.na(jurisdictionscale_data2$belowBMMSY_prec)&jurisdictionscale_data2$belowBMMSY_prec=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY_prec[!is.na(jurisdictionscale_data2$belowBMMSY_prec)])
jurisdictionscale_data2$belowBMMSY_op<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_high),NA,ifelse(jurisdictionscale_data2$larger_Bstatus_high>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY_op[!is.na(jurisdictionscale_data2$belowBMMSY_op)&jurisdictionscale_data2$belowBMMSY_op=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY_op[!is.na(jurisdictionscale_data2$belowBMMSY_op)])

#check if we assumed a country's reserves were at B0
jurisdictionscale_data2$belowBMMSY_B0<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0),NA,ifelse(jurisdictionscale_data2$larger_Bstatus_B0>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY_B0[!is.na(jurisdictionscale_data2$belowBMMSY_B0)&jurisdictionscale_data2$belowBMMSY_B0=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY_B0[!is.na(jurisdictionscale_data2$belowBMMSY_B0)])
jurisdictionscale_data2$belowBMMSY_B0_prec<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0_low),NA,ifelse(jurisdictionscale_data2$larger_Bstatus_B0_low>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY_B0_prec[!is.na(jurisdictionscale_data2$belowBMMSY_B0_prec)&jurisdictionscale_data2$belowBMMSY_B0_prec=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY_B0_prec[!is.na(jurisdictionscale_data2$belowBMMSY_B0_prec)])
jurisdictionscale_data2$belowBMMSY_B0_op<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0_high),NA,ifelse(jurisdictionscale_data2$larger_Bstatus_B0_high>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data2$belowBMMSY_B0_op[!is.na(jurisdictionscale_data2$belowBMMSY_B0_op)&jurisdictionscale_data2$belowBMMSY_B0_op=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY_B0_op[!is.na(jurisdictionscale_data2$belowBMMSY_B0_op)])

#jurisdictions COLLAPSED status: 
jurisdictionscale_data2$collapsed_B<- jurisdictionscale_data2$larger_B0*0.10
jurisdictionscale_data2$collapsed_B_prec<- jurisdictionscale_data2$larger_B0_low*0.10
jurisdictionscale_data2$collapsed_B_op<- jurisdictionscale_data2$larger_B0_high*0.10

jurisdictionscale_data2$collapsed<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0),NA, ifelse(jurisdictionscale_data2$larger_Bstatus_B0*jurisdictionscale_data2$larger_BMMSY<jurisdictionscale_data2$collapsed_B,"Collapsed","Not collapsed"))
length(jurisdictionscale_data2$collapsed[!is.na(jurisdictionscale_data2$collapsed)&jurisdictionscale_data2$collapsed=="Collapsed"])/length(jurisdictionscale_data2$collapsed[!is.na(jurisdictionscale_data2$collapsed)&jurisdictionscale_data2$belowBMMSY_B0=="belowBMMSY"])
jurisdictionscale_data2$Larger[!is.na(jurisdictionscale_data2$collapsed)&jurisdictionscale_data2$collapsed=="Collapsed"]
jurisdictionscale_data2$collapsed_prec<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0_low),NA, ifelse(jurisdictionscale_data2$larger_Bstatus_B0_low*jurisdictionscale_data2$larger_BMMSY_low<jurisdictionscale_data2$collapsed_B_prec,"Collapsed","Not collapsed"))
length(jurisdictionscale_data2$collapsed_prec[!is.na(jurisdictionscale_data2$collapsed_prec)&jurisdictionscale_data2$collapsed_prec=="Collapsed"])/length(jurisdictionscale_data2$collapsed_prec[!is.na(jurisdictionscale_data2$collapsed_prec)&jurisdictionscale_data2$belowBMMSY_B0_prec=="belowBMMSY"])
jurisdictionscale_data2$collapsed_op<- ifelse(is.na(jurisdictionscale_data2$larger_Bstatus_B0),NA, ifelse(jurisdictionscale_data2$larger_Bstatus_B0_high*jurisdictionscale_data2$larger_BMMSY_high<jurisdictionscale_data2$collapsed_B_op,"Collapsed","Not collapsed"))
length(jurisdictionscale_data2$collapsed_op[!is.na(jurisdictionscale_data2$collapsed_op)&jurisdictionscale_data2$collapsed_op=="Collapsed"])/length(jurisdictionscale_data2$collapsed_op[!is.na(jurisdictionscale_data2$collapsed_op)&jurisdictionscale_data2$belowBMMSY_B0_op=="belowBMMSY"])
jurisdictionscale_data2$Larger[!is.na(jurisdictionscale_data2$collapsed_op)&jurisdictionscale_data2$collapsed_op=="Collapsed"]
jurisdictionscale_data2$Larger[!is.na(jurisdictionscale_data2$collapsed_prec)&jurisdictionscale_data2$collapsed_prec=="Collapsed"]


#Jurisdiction Sustainability status 

#merge with jurisdiction specific data
colnames(mediansurplus_perjuris2)<- c("Larger","surplus","Biomass")
colnames(precsurplus_perjuris2)<- c("Larger","surplus_low","Biomass")
colnames(opsurplus_perjuris2)<- c("Larger","surplus_high","Biomass")
popgrowth_full_country2<- merge(merge(mediansurplus_perjuris2,precsurplus_perjuris2,by=c("Larger","Biomass")),opsurplus_perjuris2,by=c("Larger","Biomass"))

#make the biomass overlap the biomass of the surplus production curve by rounding it to the nearest integer
jurisdictionscale_data2$roundedmedianB<- round(jurisdictionscale_data2$larger_weightedB)
jurisdictionscale_data2$Biomass<- jurisdictionscale_data2$roundedmedianB

#merge with most common avergae surplus and with a countries own surplus production curve
popgrowth_full_open<- rstan::extract(Fit_full_gf,pars=c("sustyield"))[[1]]
popgrowth_full_open[is.nan(popgrowth_full_open)]<-0
popgrowth_full_open<- broom.mixed::tidyMCMC(coda::as.mcmc(popgrowth_full_open),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
#assign almost 0 surplus for biomass values that have no surplus
popgrowth_full_open$median=ifelse(popgrowth_full_open$estimate<0,0.00000001,popgrowth_full_open$estimate)
popgrowth_full_open$Biomass=seq(0,350,1)
#places that should have no surplus production assign very small (to not get inf)
popgrowth_full_open$`05%`<- ifelse(popgrowth_full_open$conf.low<0,0.000001,popgrowth_full_open$conf.low)
popgrowth_full_open$`95%`<- ifelse(popgrowth_full_open$conf.high<0,0.000001,popgrowth_full_open$conf.high)

jurisdictionscale_data2<- merge(jurisdictionscale_data2,popgrowth_full_open,by="Biomass",all.x=T)
jurisdictionscale_data2<- merge(jurisdictionscale_data2,popgrowth_full_country2,by=c("Larger","Biomass"),all.x=T)

#give almost 0 yield to those biomass values that had no surplus production based on our intercept model/ and a countries own surplus
jurisdictionscale_data2$avyield<- ifelse(is.na(jurisdictionscale_data2$median),NA,jurisdictionscale_data2$median)
jurisdictionscale_data2$highyield<- ifelse(is.na(jurisdictionscale_data2$`95%`),NA, ifelse(jurisdictionscale_data2$`95%`==0,0.0001,jurisdictionscale_data2$`95%`))
jurisdictionscale_data2$lowyield<- ifelse(is.na(jurisdictionscale_data2$`05%`),NA, ifelse(jurisdictionscale_data2$`05%`==0,0.0001,jurisdictionscale_data2$`05%`))

jurisdictionscale_data2$surplus<- ifelse(is.na(jurisdictionscale_data2$surplus),NA, ifelse(jurisdictionscale_data2$surplus<0,0.0001,jurisdictionscale_data2$surplus))
jurisdictionscale_data2$surplus_low<- ifelse(is.na(jurisdictionscale_data2$surplus_low),NA, ifelse(jurisdictionscale_data2$surplus_low<0,0.0001,jurisdictionscale_data2$surplus_low))
jurisdictionscale_data2$surplus_high<- ifelse(is.na(jurisdictionscale_data2$surplus_high),NA, ifelse(jurisdictionscale_data2$surplus_high<0,0.0001,jurisdictionscale_data2$surplus_high))

#calculate fishing status based on the surplus production curve based on average conditions and using the country surplus
jurisdictionscale_data2$fishingstatus_curve<- jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$surplus
jurisdictionscale_data2$overfishing_curve<- ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve>1,1,0))
jurisdictionscale_data2$fishingstatus_curve_prec<- jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$surplus_low
jurisdictionscale_data2$overfishing_curve_prec<- ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_prec),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1,1,0))
jurisdictionscale_data2$fishingstatus_curve_op<- jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$surplus_high
jurisdictionscale_data2$overfishing_curve_op<- ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_op),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1,1,0))
length(jurisdictionscale_data2$overfishing_curve[!is.na(jurisdictionscale_data2$overfishing_curve)& jurisdictionscale_data2$overfishing_curve==1])/length(jurisdictionscale_data2$overfishing_curve[!is.na(jurisdictionscale_data2$overfishing_curve)])
length(jurisdictionscale_data2$overfishing_curve_prec[!is.na(jurisdictionscale_data2$overfishing_curve_prec)& jurisdictionscale_data2$overfishing_curve_prec==1])/length(jurisdictionscale_data2$overfishing_curve_prec[!is.na(jurisdictionscale_data2$overfishing_curve_prec)])
length(jurisdictionscale_data2$overfishing_curve_op[!is.na(jurisdictionscale_data2$overfishing_curve_op)& jurisdictionscale_data2$overfishing_curve_op==1])/length(jurisdictionscale_data2$overfishing_curve_op[!is.na(jurisdictionscale_data2$overfishing_curve_op)])


#biomass status the same
jurisdictionscale_data2$biomassstatus_curve<- jurisdictionscale_data2$larger_Bstatus_B0
jurisdictionscale_data2$biomassstatus_curve_prec<- jurisdictionscale_data2$larger_Bstatus_B0_low
jurisdictionscale_data2$biomassstatus_curve_op<- jurisdictionscale_data2$larger_Bstatus_B0_high

#sustainability status 
jurisdictionscale_data2$fisherystatus<- as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve<1 & jurisdictionscale_data2$biomassstatus_curve>1, "In good condition", ifelse(jurisdictionscale_data2$fishingstatus_curve>1 & jurisdictionscale_data2$biomassstatus_curve<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve>1 & jurisdictionscale_data2$biomassstatus_curve>1, "Warning", "Recovering")))))
jurisdictionscale_data2$fisherystatus_prec<- as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_prec),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve_prec<1 & jurisdictionscale_data2$biomassstatus_curve_prec>1, "In good condition", ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1 & jurisdictionscale_data2$biomassstatus_curve_prec<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1 & jurisdictionscale_data2$biomassstatus_curve_prec>1, "Warning", "Recovering")))))
jurisdictionscale_data2$fisherystatus_op<- as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_op),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve_op<1 & jurisdictionscale_data2$biomassstatus_curve_op>1, "In good condition", ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1 & jurisdictionscale_data2$biomassstatus_curve_op<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1 & jurisdictionscale_data2$biomassstatus_curve_op>1, "Warning", "Recovering")))))

#(incorporating quasi-sustainable) as sustainable
jurisdictionscale_data2$fisherystatus<- ifelse(is.na(jurisdictionscale_data2$fisherystatus),NA,ifelse(jurisdictionscale_data2$fisherystatus=="Warning" &jurisdictionscale_data2$medianFstatus<1, "In good condition",as.character(jurisdictionscale_data2$fisherystatus)))
jurisdictionscale_data2$fisherystatus_prec<- ifelse(is.na(jurisdictionscale_data2$fisherystatus_prec),NA,ifelse(jurisdictionscale_data2$fisherystatus_prec=="Warning" &jurisdictionscale_data2$Fstatus_low<1, "In good condition",as.character(jurisdictionscale_data2$fisherystatus_prec)))
jurisdictionscale_data2$fisherystatus_op<- ifelse(is.na(jurisdictionscale_data2$fisherystatus_op),NA,ifelse(jurisdictionscale_data2$fisherystatus_op=="Warning" &jurisdictionscale_data2$Fstatus_high<1, "In good condition",as.character(jurisdictionscale_data2$fisherystatus_op)))

#separating only those that have both biomass and catch
jurisdictionscale_data2<-merge(jurisdictionscale_data2,country_sites,by="Larger",all.x=T)

jurisdictionscale_data3<- jurisdictionscale_data2[!is.na(jurisdictionscale_data2$fisherystatus),]
jurisdictionscale_data3<- droplevels(jurisdictionscale_data3)

#how many in their PGMY range
jurisdictionscale_data3$rel_catchpotential<-jurisdictionscale_data3$surplus/jurisdictionscale_data3$larger_MMSY
jurisdictionscale_data3$rel_catchpotential_prec<-jurisdictionscale_data3$surplus_low/jurisdictionscale_data3$larger_MMSY_low
jurisdictionscale_data3$rel_catchpotential_op<-jurisdictionscale_data3$surplus_high/jurisdictionscale_data3$larger_MMSY_high
length(jurisdictionscale_data3$rel_catchpotential[jurisdictionscale_data3$rel_catchpotential>0.8])/length(jurisdictionscale_data3$rel_catchpotential)
length(jurisdictionscale_data3$rel_catchpotential_prec[jurisdictionscale_data3$rel_catchpotential_prec>0.8])/length(jurisdictionscale_data3$rel_catchpotential_prec)
length(jurisdictionscale_data3$rel_catchpotential_op[jurisdictionscale_data3$rel_catchpotential_op>0.8])/length(jurisdictionscale_data3$rel_catchpotential_op)
length(jurisdictionscale_data3$rel_catchpotential[jurisdictionscale_data3$rel_catchpotential<0.8 &jurisdictionscale_data3$biomassstatus_curve<1])/length(jurisdictionscale_data3$rel_catchpotential)
length(jurisdictionscale_data3$rel_catchpotential_prec[jurisdictionscale_data3$rel_catchpotential_prec<0.8 &jurisdictionscale_data3$biomassstatus_curve_prec<1])/length(jurisdictionscale_data3$rel_catchpotential_prec)
length(jurisdictionscale_data3$rel_catchpotential_op[jurisdictionscale_data3$rel_catchpotential_op<0.8 &jurisdictionscale_data3$biomassstatus_curve_op<1])/length(jurisdictionscale_data3$rel_catchpotential_op)


#percentages in each category 
#in good condition
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="In good condition"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="In good condition"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="In good condition"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="In good condition"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="In good condition"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="In good condition"]

#unsustainable
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="Unsustainable"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="Unsustainable"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="Unsustainable"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="Unsustainable"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="Unsustainable"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="Unsustainable"]

#warning
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="Warning"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="Warning"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="Warning"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="Warning"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="Warning"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="Warning"]

#potentially Recovering
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="Recovering"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="Recovering"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="Recovering"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="Recovering"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="Recovering"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="Recovering"]

#have past one or both MMSY benchmarks
length(jurisdictionscale_data3$fisherystatus[!jurisdictionscale_data3$fisherystatus=="In good condition"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
length(jurisdictionscale_data3$fisherystatus_prec[!jurisdictionscale_data3$fisherystatus_prec=="In good condition"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
length(jurisdictionscale_data3$fisherystatus_op[!jurisdictionscale_data3$fisherystatus_op=="In good condition"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])

#create a data frame to store all the values
status_results<- as.data.frame(c(length(jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$catch_tkm2)]),length(jurisdictionscale_data2$belowBMMSY[!is.na(jurisdictionscale_data2$belowBMMSY)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                                 length(jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$catch_tkm2)&jurisdictionscale_data2$overfishing=="Overfishing"])/length(jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$catch_tkm2)]),
                                 length(jurisdictionscale_data2$collapsed[!is.na(jurisdictionscale_data2$collapsed)&jurisdictionscale_data2$collapsed=="Collapsed"])/length(jurisdictionscale_data2$collapsed[!is.na(jurisdictionscale_data2$collapsed)&jurisdictionscale_data2$belowBMMSY=="belowBMMSY"]),
                                 length(jurisdictionscale_data2$belowBMMSY[!is.na(jurisdictionscale_data2$belowBMMSY)&jurisdictionscale_data2$belowBMMSY=="belowBMMSY"])/length(jurisdictionscale_data2$belowBMMSY[!is.na(jurisdictionscale_data2$belowBMMSY)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="In good condition"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Unsustainable"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Warning"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Recovering"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                                 length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&!jurisdictionscale_data2$fisherystatus=="In good condition"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)])))
row.names(status_results)<- c("n_fishingstatus","n_biomassstatus","n_fisherystatus","overfishing","collapsed","belowBMMSY","In good condition","unsustainable","warning","Recovering","conservation_concern")
#write.csv(status_results, "main_status_results_bestfit_ss.csv")
colnames(status_results)<- "country_percentages"


#map of status
newmap <- getMap(resolution = "high")
#margin(t, r, l, b)

Fig.2d<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(-35, 35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data2[!is.na(jurisdictionscale_data2$overfishing),],aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data2$overfishing[!is.na(jurisdictionscale_data2$overfishing)]),colour="black", pch=21,size=3, alpha=0.9)+
  scale_fill_manual(name = "Fishing status
 (C/MMSY)",
                    values = c("Not overfishing" = "navyblue",
                               "Overfishing" = "red"),
                    labels = c("Not overfishing" ="Not overfishing", 
                               "Overfishing" ="Overfishing"), guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),panel.border = element_rect(colour = "darkgrey", fill=NA),panel.background = element_rect(fill="white",color="darkgrey")) 

Fig.2e<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group), fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(-35, 35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data2[!is.na(jurisdictionscale_data2$belowBMMSY_B0),], aes(x=EEZ_long, y=EEZ_lat, fill =jurisdictionscale_data2$belowBMMSY_B0[!is.na(jurisdictionscale_data2$belowBMMSY_B0)], size=jurisdictionscale_data2$overall_count[!is.na(jurisdictionscale_data2$belowBMMSY)]),col="black", pch=21)+
  scale_fill_manual(name = "Biomass status
 (B/BMMSY)",
                    values = c("Not belowBMMSY" = "navyblue",
                               "belowBMMSY" = "red"),
                    labels = c("Not belowBMMSY" ="Not belowBMMSY", 
                               "belowBMMSY" ="belowBMMSY"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),panel.border = element_rect(colour = "darkgrey", fill=NA),panel.background = element_rect(fill="white",color="darkgrey")) 

Fig.2f<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),  fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(-35, 35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data3, aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data3$fisherystatus, size=jurisdictionscale_data3$overall_count),col="black", pch=21,  alpha=0.9)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  scale_fill_manual (name="Fishery status",values=c( "Recovering"="cyan3",
                                                     "In good condition"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="goldenrod1"),
                     labels = c("Recovering"="Potentially
recovering",
                                "In good condition"="In good condition",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  xlab("")+
  ylab("")+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),panel.border = element_rect(colour = "darkgrey", fill=NA),panel.background = element_rect(fill="white",color="darkgrey")) 

fig2b=ggarrange(Fig.2d,Fig.2e,Fig.2f, ncol=1, nrow=3,heights=c(1,1,1),labels=c("d","e","f"))
windows()
countryspcfig<-ggplot(NULL)+geom_ribbon(aes(x=precsurplus_perjuris2$Biomass,ymin=precsurplus_perjuris2$surplus_low,ymax=opsurplus_perjuris2$surplus_high, group=precsurplus_perjuris2$Larger), fill="darkgrey", alpha=0.1)+
  geom_line(data=mediansurplus_perjuris2,aes(y=surplus, x=Biomass,group=Larger), col="grey25",alpha=0.7,lty=3)+ylim(c(0,41))+
  geom_line(data=popgrowth_full_open2 ,aes(y=median, x=Biomass), col="black",lwd=1.5,lty=1)+xlim(c(0,600))+
  theme(panel.background = element_rect(fill="white",colour="black"))+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))

windows()
juris_refpoints<-ggarrange(combinedMMSY,combinedBMMSY,countryspcfig,nrow=3,ncol=1,labels = c("a","b","c"),heights = c(1.1,1,1))
ggarrange(juris_refpoints,fig2b,nrow=1,ncol=2,widths=c(1,3))


#figure with size scale
ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),  fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(-35, 35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data3, aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data3$fisherystatus, size=jurisdictionscale_data3$overall_count),col="black", pch=21,  alpha=0.9)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  scale_fill_manual (name="Fishery status",values=c( "Recovering"="cyan3",
                                                     "In good condition"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="goldenrod1"),
                     labels = c("Recovering"="Potentially
recovering",
                                "In good condition"="In good condition",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous("Sampled sites",breaks = round(seq(round(min(jurisdictionscale_data3$overall_count),0), round(max(jurisdictionscale_data3$overall_count),0), length.out = 5),0))+
  xlab("")+
  ylab("")+theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),panel.border = element_rect(colour = "darkgrey", fill=NA),panel.background = element_rect(fill="white",color="darkgrey")) 



#correlation among fishing and biomass status
#KOBE PLOT
cor(log(jurisdictionscale_data3$larger_Bstatus_B0),log(jurisdictionscale_data3$medianFstatus),method="pearson")
kobe<-ggplot(jurisdictionscale_data3, aes(x=sqrt(larger_Bstatus_B0), y=sqrt(medianFstatus)))+ geom_pointrange(aes(ymax=sqrt(Fstatus_high),ymin=sqrt(Fstatus_low)),alpha=0.5,col="grey")+
  geom_errorbarh(aes(xmax=sqrt(larger_Bstatus_B0_high),xmin=sqrt(larger_Bstatus_B0_low)),alpha=0.5,col="grey")+
  geom_point(aes(fill=jurisdictionscale_data3$fisherystatus,size=jurisdictionscale_data3$overall_count), pch=21)+
  scale_fill_manual (name="Fishery status",values=c( "Recovering"="turquoise2",
                                                     "In good condition"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="yellow"),
                     labels = c("Recovering"="Recovering",
                                "In good condition"="In good condition",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  theme_classic()+xlab("sqrt(weightedB/BMMSY)")+ylab("sqrt(C/MMSY)")+geom_hline(yintercept=1,lty=2,col="black")+geom_vline(xintercept=1,lty=2,col="black")+geom_text_repel(aes(label=jurisdictionscale_data3$Larger),size=2.8)

#supplementary figure of status
windows()
ggarrange(prob_aboveMMSY_fig,prob_belowBMMSY_fig_all,kobe,nrow=1,ncol=3,labels=c("a","","d"))


#######################################################################################################################################

#high gravity locations have difefrent environmental conditions an dthus B0 (ongoing human induced environmental change?)
f<- ggplot(ordereddata,aes(x=sgravtot,y=site_BMMSY))+geom_point(pch=21,fill="grey",alpha=0.1)+geom_smooth(col="black",method="gam")+xlab("log(Total gravity+1)")+ylab("Site-specific BMMSY (t/km2)")+theme(text = element_text(family="Helvetica"),panel.background = element_rect(fill="white",color="black"))+ labs(y =  expression ("Site-specific "*B["MMSY "]*"("~t/km^2*")"))+xlim(c(0,11))
g<- ggplot(ordereddata,aes(x=sgravtot,y=site_MMSY))+geom_point(pch=21,fill="grey",alpha=0.1)+geom_smooth(col="black",method="gam")+xlab("log(Total gravity+1)")+ylab("Site-specific MMSY (t/km2/y)")+theme(text = element_text(family="Helvetica"),panel.background = element_rect(fill="white",color="black"))+ labs(y =  expression ("Site-specific MMSY ("~t/km^2/y*")"))+xlim(c(0,11))
windows()
ggarrange(f,g,nrow=1,ncol=2,labels=c("a","b"))

#######################################################################################################################################
#relation between biomass and fishing status
fisheddata2<- merge(extracteddata,jurisdictionscale_data2, by="Larger", all.x=T)

site_kobe<- ggplot(fisheddata2, aes(x=sqrt(site_status), y=sqrt(site_catchstatusMMSY)))+guides(fill=F)+ 
  geom_pointrange(aes(ymax=sqrt(site_catchstatusMMSY_prec),ymin=sqrt(site_catchstatusMMSY_op)),alpha=0.1,col="grey")+
  geom_errorbarh(aes(xmax=sqrt(site_status_op),xmin=sqrt(site_status_prec)),alpha=0.1,col="grey")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha=0.4)+scale_fill_gradient(low="plum",  high="darkmagenta")+geom_point( col="black",alpha=0.1)+
  theme_classic()+xlab("Site B/BMMSY")+ylab("Site C/MMSY")+geom_hline(yintercept=1,lty=2,col="black",lwd=1)+geom_vline(xintercept=1,lty=2,col="black",lwd=1)


kobe_notext<- ggplot(jurisdictionscale_data3, aes(x=sqrt(larger_Bstatus_B0), y=sqrt(medianFstatus)))+ 
  geom_pointrange(aes(ymax=sqrt(Fstatus_high),ymin=sqrt(Fstatus_low)),alpha=0.5,col="grey")+
  geom_errorbarh(aes(xmax=sqrt(larger_Bstatus_B0_high),xmin=sqrt(larger_Bstatus_B0_low)),alpha=0.5,col="grey")+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha=0.5)+scale_fill_gradient(low="white",  high="darkmagenta", guide=FALSE)+
  geom_point(aes(size=jurisdictionscale_data3$overall_count),col="black",alpha=0.3)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  theme_classic()+xlab("Jurisdiction weightedB/BMMSY")+ylab("Jurisdiction C/MMSY")+geom_hline(yintercept=1,lty=2,lwd=1.1,col="black")+geom_vline(xintercept=1,lty=2,lwd=1.1,col="black")
windows()
ggarrange(site_kobe,kobe_notext,labels=c("a","b"))

cor(sqrt(jurisdictionscale_data3$larger_Bstatus_B0),sqrt(jurisdictionscale_data3$medianFstatus),method="pearson")

fisheddata2$belowBMMSY_numeric<- ifelse(fisheddata2$site_status<1,1,0 )
extracteddata2<- fisheddata2

#trying a combined measure (F/fmsy patterns): 
fisheddata2$f_fmmsy<-fisheddata2$site_catchstatusMMSY/fisheddata2$site_status
jurisdictionscale_data2$f_fmmsy_j<-jurisdictionscale_data2$medianFstatus/jurisdictionscale_data2$larger_Bstatus_B0
jurisdictionscale_data3$f_fmmsy_j<-jurisdictionscale_data3$medianFstatus/jurisdictionscale_data3$larger_Bstatus_B0
site_kobe2<- ggplot(fisheddata2, aes(x=ifelse(site_status>3,3,site_status), y=ifelse(f_fmmsy>3,3,f_fmmsy)))+guides(fill=F)+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha=0.4)+scale_fill_gradient(low="plum",  high="darkmagenta")+geom_point( col="black",alpha=0.1)+
  theme_classic()+xlab("Site B/BMMSY")+ylab("Site F/FMMSY")+geom_hline(yintercept=1,lty=2,col="black",lwd=1)+geom_vline(xintercept=1,lty=2,col="black",lwd=1)+xlim(c(0,3))+ylim(c(0,3))

kobe_notext2<- ggplot(jurisdictionscale_data3, aes(x=ifelse(larger_Bstatus_B0>3,3,larger_Bstatus_B0), y=ifelse(f_fmmsy_j>3,3,f_fmmsy_j)))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha=0.5)+scale_fill_gradient(low="white",  high="darkmagenta", guide=FALSE)+xlim(c(0,3))+ylim(c(0,3))+
  geom_point(aes(size=jurisdictionscale_data3$overall_count),col="black",alpha=0.3)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+ geom_text_repel(aes(label=Larger2))+
  theme_classic()+xlab("Jurisdiction B/BMMSY")+ylab("Jurisdiction F/FMMSY")+geom_hline(yintercept=1,lty=2,lwd=1.1,col="black")+geom_vline(xintercept=1,lty=2,lwd=1.1,col="black")
windows()
ggarrange(site_kobe2,kobe_notext2, widths = c(1,2))

#######################################################################################################################################
#Trade-offs between production and other ecosystem metrics ............................................................................
#only fished ddata
alldata<- extracteddata
alldata<- droplevels(alldata)

#Tranform (if neccesary) ecosystem metrics and run bayesian hierarchical models using the same covariates as we used for fished biomass 
#Mean fish length .....................................................................................................................
hist(log(alldata$meanSize_cm))
alldata$lSize<- log(alldata$meanSize_cm)

model_Size<- brm(lSize~ 
                   sDepth+sOcean_prod+Atoll+sSST+
                   ReefHabitat+
                   sSampArea+
                   SampMethod+sgravtot+
                   (1|Larger),data=alldata[!is.na(alldata$lSize),],family=gaussian,
                 iter=10000,  warmup=9000,
                 chains=4)

model_Size_Null<- brm(lSize~ 1,data=alldata[!is.na(alldata$lSize),],family=gaussian,
                      iter=10000,  warmup=9000,
                      chains=4)

#model comparison
size_loo<- loo(model_Size,model_Size_Null)
emetrics_loo<- as.data.frame(size_loo$diffs)

#model fit
residualsS<- resid(model_Size)[,1]
a<- ggplot(NULL)+geom_histogram(aes(x=residualsS))+theme_classic()+xlab ("Residuals mean length model")
b<- ggplot(NULL)+geom_point(aes(x=fitted(model_Size)[,1],y=resid(model_Size)[,1]),pch=21,fill="darkgrey",alpha=0.4)+theme_classic()+xlab ("Fitted mean length (log)")+ylab ("Residuals")
c<- pp_check(model_Size,nsamples=4000)+xlab("Mean length cm (log)")

#marginalized response variable
alldata$Size_marg<- alldata$meanSize_cm/exp((fixef(model_Size)["sSampArea",1]*alldata$sSampArea+
                                               fixef(model_Size)["sSST",1]*alldata$sSST+
                                               ifelse(alldata$ReefHabitat=="Crest",fixef(model_Size)[c("ReefHabitatCrest"),1],
                                                      ifelse(alldata$ReefHabitat=="Flat",fixef(model_Size)[c("ReefHabitatFlat"),1],
                                                             ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                    fixef(model_Size)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                               fixef(model_Size)["sDepth",1]*alldata$sDepth+ fixef(model_Size)["sOcean_prod",1]*alldata$sOcean_prod+
                                               fixef(model_Size)["Atoll",1]*alldata$Atoll+
                                               ifelse(alldata$SampMethod=="Distance sampling",fixef(model_Size)[c("SampMethodDistancesampling")],
                                                      ifelse(alldata$SampMethod=="Point intercept",fixef(model_Size)[c("SampMethodPointintercept")],0))))

alldata$Size_marg_meth<- alldata$meanSize_cm/exp((fixef(model_Size)["sSampArea",1]*alldata$sSampArea+
                                                    ifelse(alldata$ReefHabitat=="Crest",fixef(model_Size)[c("ReefHabitatCrest"),1],
                                                           ifelse(alldata$ReefHabitat=="Flat",fixef(model_Size)[c("ReefHabitatFlat"),1],
                                                                  ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                         fixef(model_Size)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                                    fixef(model_Size)["sDepth",1]*alldata$sDepth+ 
                                                    ifelse(alldata$SampMethod=="Distance sampling",fixef(model_Size)[c("SampMethodDistancesampling")],
                                                           ifelse(alldata$SampMethod=="Point intercept",fixef(model_Size)[c("SampMethodPointintercept")],0))))

#Probability of observing top predators ...............................................................................................
hist(log(alldata$TPBiomass_tkm2+1))
length(alldata$TPBiomass_tkm2[alldata$TPBiomass_tkm2==0])/length(alldata$TPBiomass_tkm2)
length(alldata$TPBiomass_tkm2[!is.na(alldata$TPBiomass_tkm2)])

alldata$PA_tp<- ifelse(is.na(alldata$TPBiomass_tkm2),NA, ifelse(alldata$TPBiomass_tkm2>0,1,0))
reefscale_data$PA_tp<- ifelse(is.na(reefscale_data$TPBiomass_tkm2),NA, ifelse(reefscale_data$TPBiomass_tkm2>0,1,0))

model_tp<- brm(PA_tp~ 
                 sDepth+sOcean_prod+Atoll+sSST+
                 ReefHabitat+
                 sSampArea+
                 SampMethod+sgravtot+
                 (1|Larger),data=alldata[!is.na(alldata$PA_tp),],
               family=bernoulli("logit"), iter=10000,  warmup=9000,
               chains=4)

model_tp_Null<- brm(PA_tp~ 1,data=alldata[!is.na(alldata$PA_tp),],
                    family=bernoulli("logit"), iter=10000,  warmup=9000,
                    chains=4)

#model selection
tp_loo<- loo(model_tp,model_tp_Null)
emetrics_loo<- rbind(emetrics_loo,tp_loo$diffs)

#model fit
yrep_PA_tp<-as.data.frame(posterior_predict(model_tp))
x <- createDHARMa(simulatedResponse=t(yrep_PA_tp), observedResponse=alldata$PA_tp[!is.na(alldata$PA_tp)])
d<-ggplot(NULL)+geom_histogram(aes(x=x$scaledResiduals))+xlab("Scaled residuals top predator model")+theme_classic()
e<-ggplot(data=NULL)+
  geom_point(aes(y=x$scaledResiduals, x=x$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")
f<- pp_check(model_tp,nsamples=4000)+xlab("Presence/absence top predators")

#marginalized response variable
alldata$tp_marg<- alldata$PA_tp-(fixef(model_tp)["sSampArea",1]*alldata$sSampArea+
                                   fixef(model_tp)["sSST",1]*alldata$sSST+
                                   ifelse(alldata$ReefHabitat=="Crest",fixef(model_tp)[c("ReefHabitatCrest"),1],
                                          ifelse(alldata$ReefHabitat=="Flat",fixef(model_tp)[c("ReefHabitatFlat"),1],
                                                 ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                        fixef(model_tp)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                   fixef(model_tp)["sDepth",1]*alldata$sDepth+ fixef(model_tp)["sOcean_prod",1]*alldata$sOcean_prod+
                                   fixef(model_tp)["Atoll",1]*alldata$Atoll+
                                   ifelse(alldata$SampMethod=="Distance sampling",fixef(model_tp)[c("SampMethodDistancesampling")],
                                          ifelse(alldata$SampMethod=="Point intercept",fixef(model_tp)[c("SampMethodPointintercept")],0)))
alldata$tp_marg_meth<- alldata$PA_tp-(fixef(model_tp)["sSampArea",1]*alldata$sSampArea+
                                        ifelse(alldata$ReefHabitat=="Crest",fixef(model_tp)[c("ReefHabitatCrest"),1],
                                               ifelse(alldata$ReefHabitat=="Flat",fixef(model_tp)[c("ReefHabitatFlat"),1],
                                                      ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                             fixef(model_tp)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                        fixef(model_tp)["sDepth",1]*alldata$sDepth+ 
                                        ifelse(alldata$SampMethod=="Distance sampling",fixef(model_tp)[c("SampMethodDistancesampling")],
                                               ifelse(alldata$SampMethod=="Point intercept",fixef(model_tp)[c("SampMethodPointintercept")],0)))
#change from log-odds to probability and then back to presence absence
alldata$prob_tpmarg<- exp(alldata$tp_marg)/(1+exp(alldata$tp_marg))
alldata$PA_tpmarg<- ifelse(alldata$prob_tpmarg>0.5,1,0)
alldata$prob_tpmarg_meth<- exp(alldata$tp_marg_meth)/(1+exp(alldata$tp_marg_meth))
alldata$PA_tpmarg_meth<- ifelse(alldata$prob_tpmarg_meth>0.5,1,0)

#Total fish species richness ..........................................................................................................
hist(log(alldata$total_sp))
alldata$ltotalsp<- log(alldata$total_sp)
length(alldata$total_sp[!is.na(alldata$total_sp)])
model_totalsp<- brm(ltotalsp~ 
                      sDepth+sOcean_prod+Atoll+sSST+
                      ReefHabitat+
                      sSampArea+
                      SampMethod+sgravtot+
                      (1|Larger),data=alldata[!is.na(alldata$ltotalsp),],family=gaussian,
                    iter=10000,  warmup=9000,
                    chains=4)

model_totalsp_Null<- brm(ltotalsp~ 1,data=alldata[!is.na(alldata$ltotalsp),],family=gaussian,
                         iter=10000,  warmup=9000,
                         chains=4)

#model selection
sp_loo<- loo(model_totalsp,model_totalsp_Null)
emetrics_loo<- rbind(emetrics_loo,sp_loo$diffs)

#model fit
residualsts<- resid(model_totalsp)[,1]
g<- ggplot(NULL)+geom_histogram(aes(x=residualsts))+theme_classic()+ xlab ("Residuals total species richness model")
h<- ggplot(NULL)+geom_point(aes(x=fitted(model_totalsp)[,1],y=resid(model_totalsp)[,1]),pch=21,fill="darkgrey",alpha=0.4)+theme_classic()+ ylab ("Residuals")+ xlab ("Fitted total species richness model (log)")
i<- pp_check(model_totalsp,nsamples=4000)+xlab("Total fish sp. richness (log)")

#Marginalized response variable
alldata$totalsp_marg<- alldata$total_sp/exp((fixef(model_totalsp)["sSampArea",1]*alldata$sSampArea+
                                               fixef(model_totalsp)["sSST",1]*alldata$sSST+
                                               ifelse(alldata$ReefHabitat=="Crest",fixef(model_totalsp)[c("ReefHabitatCrest"),1],
                                                      ifelse(alldata$ReefHabitat=="Flat",fixef(model_totalsp)[c("ReefHabitatFlat"),1],
                                                             ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                    fixef(model_totalsp)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                               fixef(model_totalsp)["sDepth",1]*alldata$sDepth+ fixef(model_totalsp)["sOcean_prod",1]*alldata$sOcean_prod+
                                               fixef(model_totalsp)["Atoll",1]*alldata$Atoll+
                                               ifelse(alldata$SampMethod=="Distance sampling",fixef(model_totalsp)[c("SampMethodDistancesampling")],
                                                      ifelse(alldata$SampMethod=="Point intercept",fixef(model_totalsp)[c("SampMethodPointintercept")],0))))

alldata$totalsp_marg_meth<- alldata$total_sp/exp((fixef(model_totalsp)["sSampArea",1]*alldata$sSampArea+
                                                    ifelse(alldata$ReefHabitat=="Crest",fixef(model_totalsp)[c("ReefHabitatCrest"),1],
                                                           ifelse(alldata$ReefHabitat=="Flat",fixef(model_totalsp)[c("ReefHabitatFlat"),1],
                                                                  ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                         fixef(model_totalsp)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                                    fixef(model_totalsp)["sDepth",1]*alldata$sDepth+ 
                                                    ifelse(alldata$SampMethod=="Distance sampling",fixef(model_totalsp)[c("SampMethodDistancesampling")],
                                                           ifelse(alldata$SampMethod=="Point intercept",fixef(model_totalsp)[c("SampMethodPointintercept")],0))))


#Parrotfish scapring potential .........................................................................................................
length(alldata$scraping_potential[!is.na(alldata$scraping_potential)])
model_herb<- brm(scraping_potential~ 
                   sDepth+sOcean_prod+Atoll+sSST+
                   ReefHabitat+
                   sSampArea+
                   SampMethod+sgravtot+
                   (1|Larger),data=alldata[!is.na(alldata$scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                                                                                                        link_hu = "logit"),
                 iter=10000,  warmup=9000,
                 chains=4)

model_herb_Null<- brm(scraping_potential~ 
                        1,data=alldata[!is.na(alldata$scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                                                                                                    link_hu = "logit"),
                      iter=10000,  warmup=9000,
                      chains=4)

#model selection
her_loo<- loo(model_herb,model_herb_Null)
emetrics_loo<- rbind(emetrics_loo,her_loo$diffs)
#write.csv(emetrics_loo, "loo_ecosystemmetrics.csv")

#model fit
residualsC<- resid(model_herb)[,1]
yrep_her<-as.data.frame(posterior_predict(model_herb))
x <- createDHARMa(simulatedResponse=t(yrep_her), observedResponse=alldata$scraping_potential[!is.na(alldata$scraping_potential)])
j<-ggplot(NULL)+geom_histogram(aes(x=x$scaledResiduals))+xlab("Scaled residuals parrotfish scraping model")+theme_classic()
k<-ggplot(data=NULL)+
  geom_point(aes(y=x$scaledResiduals, x=x$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")+theme_classic()
l<- pp_check(model_herb,nsamples=4000)+xlab("Parrotfish scraping potential")+xlim(c(0,3000))

#marginalized response variable
alldata$herb_marg<- alldata$scraping_potential/exp(fixef(model_herb)["sSampArea",1]*alldata$sSampArea+
                                                     fixef(model_herb)["sSST",1]*alldata$sSST+
                                                     ifelse(alldata$ReefHabitat=="Crest",fixef(model_herb)[c("ReefHabitatCrest"),1],
                                                            ifelse(alldata$ReefHabitat=="Flat",fixef(model_herb)[c("ReefHabitatFlat"),1],
                                                                   ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                          fixef(model_herb)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                                     fixef(model_herb)["sDepth",1]*alldata$sDepth+ fixef(model_herb)["sOcean_prod",1]*alldata$sOcean_prod+
                                                     fixef(model_herb)["Atoll",1]*alldata$Atoll+
                                                     ifelse(alldata$SampMethod=="Distance sampling",fixef(model_herb)[c("SampMethodDistancesampling")],
                                                            ifelse(alldata$SampMethod=="Point intercept",fixef(model_herb)[c("SampMethodPointintercept")],0)))

alldata$herb_marg_meth<- alldata$scraping_potential/exp(fixef(model_herb)["sSampArea",1]*alldata$sSampArea+
                                                          ifelse(alldata$ReefHabitat=="Crest",fixef(model_herb)[c("ReefHabitatCrest"),1],
                                                                 ifelse(alldata$ReefHabitat=="Flat",fixef(model_herb)[c("ReefHabitatFlat"),1],
                                                                        ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                                               fixef(model_herb)[c("ReefHabitatLagoon_Backreef"),1],0)))+
                                                          fixef(model_herb)["sDepth",1]*alldata$sDepth+ 
                                                          
                                                          ifelse(alldata$SampMethod=="Distance sampling",fixef(model_herb)[c("SampMethodDistancesampling")],
                                                                 ifelse(alldata$SampMethod=="Point intercept",fixef(model_herb)[c("SampMethodPointintercept")],0)))

#supporting figure of model fits of ecosystem metric models
windows()
ggarrange(a,b,c,d,e,f,g,h,i,j,k,l,nrow=4,ncol=3,labels=c("a","b","c","d","e","f","g","h","i","j","k","l"))


#correlation of metrics (and also correlation with mean_length)
windows()
pairs(~log(herb_marg+1)+PA_tpmarg+
        log(totalsp_marg)+Size_marg ,  data= alldata,lower.panel=panel.cor, pch = 21, bg = "darkgrey",labels=c("log(Parrotfish
         scraping potential+1)", "Presence/absence 
         top predators","log(Total species 
         richness)", "Mean length"),cex.labels=1.5,font.labels=2,diag.panel =panel.hist, hist.col="grey" )


#Trade-offs ...........................................................................................................................

#with surplus production at most common and average environmental conditions
#median conditions
MMSY_median<- median(rstan::extract(Fit_full_gf, pars="MMSY")$MMSY)
BMMSY_median<- median(rstan::extract(Fit_full_gf, pars="BMMSY")$BMMSY)
B0_median<- median(rstan::extract(Fit_full_gf, pars="B0")$B0)

#"Pretty good multispecies yield" for median conditions
PGY<- 0.8*MMSY_median
B_pgy<- max(popgrowth_full_open$Biomass[popgrowth_full_open$median>PGY], na.rm=T)
B_pgy_l<- min(popgrowth_full_open$Biomass[popgrowth_full_open$median>PGY], na.rm=T)

#transform the biomass that is marginalized to be comparable with surplus production and ecosystem metrics
alldata$lB_marg<- log(alldata$site_Bmarg_all)

#individual metric relationships 
si_1<- ggplot(alldata)+
  geom_point(aes(lB_marg,Size_marg,shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="darkgreen") +
  scale_alpha_manual(values=c("Openly Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,Size_marg),formula=y~s(x,k=3),method="gam", col="darkgreen",fill="green",alpha=0.3)+
  theme_classic()+xlab("")+ylab("Mean length
 (cm)")+scale_y_continuous(breaks = seq(10, 40, by = 5),limits=c(10,40),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))


si_2<- ggplot(alldata)+
  geom_point(aes(lB_marg,PA_tpmarg,shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="navyblue") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,PA_tpmarg),formula=y~s(x,k=3),method="gam", col="navyblue",fill="blue",alpha=0.3)+
  theme_classic()+xlab("")+ylab("Prob. top 
 predators")+scale_y_continuous(breaks = seq(0,1, by =1),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))



si_3<- ggplot(alldata)+
  geom_point(aes(lB_marg,log(totalsp_marg),shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="goldenrod") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,log(totalsp_marg)),formula=y~s(x,k=3),method="gam", col="goldenrod",fill="yellow",alpha=0.3)+
  theme_classic()+xlab("")+ylab("log(Sp richness)")+scale_y_continuous(breaks = seq(0,8, by = 2),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))


si_4<- ggplot(alldata)+
  geom_point(aes(lB_marg,log(herb_marg+1),shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="violetred4") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,log(herb_marg+1)),formula=y~s(x,k=3),method="gam", col="violetred4",fill="violetred",alpha=0.3)+
  theme_classic()+xlab("log(Biomass (t/km2))")+ylab("log(Scraping potential
 (cm^2/min)+1)")+labs(x= expression ("log(Biomass ("~t/km^2*"))"))+scale_y_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10),labels = scales::number_format(accuracy = 1),position = "left")#+coord_cartesian( expand=c(0,0))

#Ecosystem metric distributions by management
# we overlaid the observed distributions for remote reefs
a<- ggplot(NULL) + 
  geom_density(aes(x=reefscale_data$meanSize_cm[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0]),lty=2)+
  geom_density(aes(x=alldata$Size_marg, fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="darkgreen",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="green","Restricted"="darkgreen"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+ scale_x_continuous(breaks = seq(10, 40, by = 5),limits=c(10,40))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
b<- ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$total_sp[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0])),lty=2)+
  geom_density(aes(x=log(alldata$totalsp_marg), fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="goldenrod",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="yellow","Restricted"="goldenrod"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(1, 8, by = 2))+xlim(c(1,8))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
c<- ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$scraping_potential[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0])+1),lty=2)+
  geom_density(aes(x=log(alldata$herb_marg+1), fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="violetred4",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="violetred2","Restricted"="violetred4"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
summary(alldata$herb_marg[as.factor(alldata$status_management)=="Openly fished"])
summary(alldata$herb_marg[as.factor(alldata$status_management)=="Restricted"])
e<- ggplot(NULL) + geom_density(aes(x=reefscale_data$PA_tp[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0]),lty=2)+
  
  geom_density(aes(x=alldata$PA_tpmarg,fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="navyblue",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="blue","Restricted"="navyblue"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 1, by = 1))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
f<- ggplot(NULL) + geom_density(aes(x=log(reefscale_data$FamBiomass_tkm2[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0] )),lty=2)+
  
  geom_density(aes(log(alldata$site_Bmarg_all),fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="darkred",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="red","Restricted"="darkred"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+#xlim(0,300)+
  labs(x="log(Biomass (t/km2))", y="") +theme_classic()+coord_flip()

#supplementary figure
windows()
ggarrange(si_1,a,si_2,e,si_3,b,si_4,c, nrow=4,ncol=2,widths=c(1.25,1), labels = c("a","e", "b","f","c","g", "d","h"))

#do trends with site-status (i.e., use metrics only marginalized by methodology)
si_1s<- ggplot(alldata)+
  geom_point(aes(log(site_status),Size_marg_meth,shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="darkgreen") +
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(log(site_status),Size_marg_meth),formula=y~s(x,k=3),method="gam", col="darkgreen",fill="green",alpha=0.3)+
  theme_classic()+xlab("")+ylab("Mean length
 (cm)")+scale_y_continuous(breaks = seq(10, 35, by = 5),limits=c(10,35),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))+
  geom_vline(xintercept = 0,lty=2)
si_2s<- ggplot(alldata)+
  geom_point(aes(log(site_status),PA_tpmarg_meth,shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="navyblue") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(log(site_status),PA_tpmarg_meth),formula=y~s(x,k=3),method="gam", col="navyblue",fill="blue",alpha=0.3)+
  theme_classic()+xlab("")+ylab("Prob. top 
 predators")+scale_y_continuous(breaks = seq(0,1, by =1),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))+
  geom_vline(xintercept = 0,lty=2)
si_3s<- ggplot(alldata)+
  geom_point(aes(log(site_status),log(totalsp_marg_meth),shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="goldenrod") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(log(site_status),log(totalsp_marg_meth)),formula=y~s(x,k=3),method="gam", col="goldenrod",fill="yellow",alpha=0.3)+
  theme_classic()+xlab("")+ylab("log(Sp richness)")+scale_y_continuous(breaks = seq(0,8, by = 2),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))+
  geom_vline(xintercept = 0,lty=2)
si_4s<- ggplot(alldata)+
  geom_point(aes(log(site_status),log(herb_marg_meth+1),shape=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="violetred4") +
  scale_shape_manual(values=c("Openly fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Openly fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(log(site_status),log(herb_marg_meth+1)),formula=y~s(x,k=3),method="gam", col="violetred4",fill="violetred",alpha=0.3)+
  theme_classic()+scale_y_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10),labels = scales::number_format(accuracy = 1),position = "left")+#coord_cartesian( expand=c(0,0))+
  geom_vline(xintercept = 0,lty=2)+ylab("log(Scraping potential
 (cm^2/min)+1)")+xlab("log(B/BMMSY)")

a<- ggplot(NULL) + 
  geom_density(aes(x=reefscale_data$meanSize_cm[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0]),lty=2)+
  geom_density(aes(x=alldata$Size_marg_meth, fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="darkgreen",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="green","Restricted"="darkgreen"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+ scale_x_continuous(breaks = seq(10, 35, by = 5),limits=c(10,35))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
b<- ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$total_sp[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0])),lty=2)+
  geom_density(aes(x=log(alldata$totalsp_marg_meth), fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="goldenrod",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="yellow","Restricted"="goldenrod"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(1, 8, by = 2))+xlim(c(1,8))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
c<- ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$scraping_potential[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0])+1),lty=2)+
  geom_density(aes(x=log(alldata$herb_marg_meth+1), fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="violetred4",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="violetred2","Restricted"="violetred4"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
e<- ggplot(NULL) + geom_density(aes(x=reefscale_data$PA_tp[reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0]),lty=2)+
  
  geom_density(aes(x=alldata$PA_tpmarg_meth,fill=as.factor(alldata$status_management),alpha=as.factor(alldata$status_management)),col="navyblue",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="blue","Restricted"="navyblue"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 1, by = 1))+
  labs(x="", y="") +
  coord_flip()+theme_classic()

windows()
ggarrange(si_1s,a,si_2s,e,si_3s,b,si_4s,c, nrow=4,ncol=2,widths=c(1.25,1), labels = c("a","e", "b","f","c","g", "d","h"))


#now we get the distributions covering the range in our suplus production curves for the main (up to 300 t/km2)
alldata2<- alldata[alldata$site_Bmarg_all<300,]

#distributions
a<- ggplot(NULL) + 
  geom_density(aes(x=alldata2$Size_marg_meth, fill=as.factor(alldata2$status_management),alpha=as.factor(alldata2$status_management)),col="darkgreen",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="green","Restricted"="darkgreen"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+xlim(c(10,40))+
  guides(alpha=F,fill=F)+ # scale_x_continuous(breaks = seq(10, 40, by = 5))+
  labs(x="mean Length (cm)", y="") +
  coord_flip()+theme_classic()
b<- ggplot(NULL) + 
  geom_density(aes(x=alldata2$totalsp_marg_meth, fill=as.factor(alldata2$status_management),alpha=as.factor(alldata2$status_management)),col="goldenrod",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="yellow","Restricted"="goldenrod"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  xlim(c(0,200))+
  labs(x="Total species richness", y="") +
  coord_flip()+theme_classic()
d<- ggplot(NULL) + 
  geom_density(aes(x=alldata2$herb_marg_meth, fill=as.factor(alldata2$status_management),alpha=as.factor(alldata2$status_management)),col="violetred4",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="violetred2","Restricted"="violetred4"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+xlim(c(0,750))+
  guides(alpha=F,fill=F)+  #scale_x_continuous(breaks = seq(0, 10, by = 2))+
  labs(x="Parrotfish scrapping potential (cm^2/min)", y="") +
  coord_flip()+theme_classic()

median(alldata$herb_marg_meth[alldata$status_management=="Restricted"],na.rm=T)/ median(alldata$herb_marg_meth[alldata$status_management=="Openly fished"],na.rm=T)
median(alldata2$herb_marg_meth[alldata2$status_management=="Restricted"],na.rm=T)/ median(alldata2$herb_marg_meth[alldata2$status_management=="Openly fished"],na.rm=T)

e<- ggplot(NULL) + 
  geom_density(aes(x=alldata2$PA_tpmarg_meth,fill=as.factor(alldata2$status_management),alpha=as.factor(alldata2$status_management)),col="navyblue",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="blue","Restricted"="navyblue"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 1, by = 1))+
  labs(x="Prob. encountering top predators", y="") +
  coord_flip()+theme_classic()
f<- ggplot(NULL) + 
  geom_density(aes(alldata2$site_Bmarg,fill=as.factor(alldata2$status_management),alpha=as.factor(alldata2$status_management)),col="darkred",  
               na.rm=T) + scale_fill_manual(values=c("Openly fished"="red","Restricted"="darkred"))+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+#xlim(0,300)+
  labs(x="Biomass (t/km2)", y="") +theme_classic()+coord_flip()

#part of figure 3
windows()
ggarrange(a,b,d,e,f,nrow=2, ncol=3)

#Generalized additive models to show the relationships between ecosystem metrics and biomass
#Mean fish length ........................................ ............................................................................
gam_size<- gam(Size_marg~s(lB_marg,k=3),data=alldata)
#model fit
ggplot(NULL)+geom_histogram(aes(x=resid(gam_size)))+xlab ("Residuals mean length gam")  
#predicted length
size_pred<- as.data.frame(predict(gam_size,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
size_pred$lB_marg<- sort(alldata$lB_marg)
size_pred$upr<- size_pred$fit+(2*size_pred$se.fit)
size_pred$low<- size_pred$fit-(2*size_pred$se.fit)
#biomass arithmetic scale (for the trade-offs graph)
size_pred$Biomass<- exp(size_pred$lB_marg)
#Length at the different reference points
currentB<- median(alldata$site_Bmarg)
ScurrentB<- size_pred$fit[size_pred$Biomass>currentB][1]
SBMSY<- size_pred$fit[size_pred$Biomass>BMMSY_median][1]
SPGY<- size_pred$fit[size_pred$Biomass>B_pgy][1]
SPGY_l<- size_pred$fit[size_pred$Biomass>B_pgy_l][1]
SBO<- size_pred$fit[size_pred$Biomass>B0_median][1]
#change with respect to unOpenly fished biomass conditions (unOpenly fished biomass 100%)
100-((SBMSY*100)/SBO)
(100-(SPGY_l*100)/SBO)-(100-((SBMSY*100)/SBO))
(100-(SBMSY*100)/SBO)-(100-(SPGY*100)/SBO)

#Probability of observing top predators .................. ............................................................................
gam_tp<- gam(PA_tpmarg~s(lB_marg,k=3),data=alldata)
plot(fitted(gam_tp), resid(gam_tp))  
tp_pred<- as.data.frame(predict(gam_tp,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
tp_pred$lB_marg<- sort(alldata$lB_marg)
tp_pred$upr<- tp_pred$fit+(2*tp_pred$se.fit)
tp_pred$low<- tp_pred$fit-(2*tp_pred$se.fit)
tp_pred$Biomass<- exp(tp_pred$lB_marg)
pBMSY<- tp_pred$fit[tp_pred$Biomass>BMMSY_median][1]
pPGY<- tp_pred$fit[tp_pred$Biomass>B_pgy][1]
pBO<- tp_pred$fit[tp_pred$Biomass>B0_median][1]
pcurrentB<- tp_pred$fit[tp_pred$Biomass>currentB][1]
#percentages
pPGY_l<- tp_pred$fit[tp_pred$Biomass>B_pgy_l][1]

#change with respect to unOpenly fished biomass conditions (unOpenly fished biomass 100%)
100-((pBMSY*100)/pBO)
(100-(pPGY_l*100)/pBO)-(100-((pBMSY*100)/pBO))
(100-(pBMSY*100)/pBO)-(100-(pPGY*100)/pBO)


#Total fish species richness .............................. ............................................................................
gam_tr<- gam(log(totalsp_marg)~s(lB_marg,k=3),data=alldata)
hist(resid(gam_tr))  
tr_pred<- as.data.frame(predict(gam_tr,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
tr_pred$lB_marg<- sort(alldata$lB_marg)
tr_pred$fit<- exp(tr_pred$fit)
tr_pred$upr<- tr_pred$fit+(2*exp(tr_pred$se.fit))
tr_pred$low<- tr_pred$fit-(2*exp(tr_pred$se.fit))
tr_pred$Biomass<- exp(tr_pred$lB_marg)
rBMSY<- tr_pred$fit[tr_pred$Biomass>BMMSY_median][1]
rPGY<- tr_pred$fit[tr_pred$Biomass>B_pgy][1]
rPGY_l<- tr_pred$fit[tr_pred$Biomass>B_pgy_l][1]
rBO<- tr_pred$fit[tr_pred$Biomass>B0_median][1]
rcurrentB<- tr_pred$fit[tr_pred$Biomass>currentB][1]
#percentages
100-((rBMSY*100)/rBO)
(100-(rPGY_l*100)/rBO)-(100-((rBMSY*100)/rBO))
(100-(rBMSY*100)/rBO)-(100-(rPGY*100)/rBO)

#Parrotfish scraping potential ............................ ............................................................................
gam_her<- gam(log(herb_marg+1)~s(lB_marg,k=3),data=alldata)
hist(resid(gam_her))  
her_pred<- as.data.frame(predict(gam_her,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
her_pred$lB_marg<- sort(alldata$lB_marg)
her_pred$fit<- exp(her_pred$fit)
her_pred$upr<- her_pred$fit+(2*exp(her_pred$se.fit))
her_pred$low<- her_pred$fit-(2*exp(her_pred$se.fit))
her_pred$Biomass<- exp(her_pred$lB_marg)
hBMSY<- her_pred$fit[her_pred$Biomass>BMMSY_median][1]
hPGY<- her_pred$fit[her_pred$Biomass>B_pgy][1]
hPGY_l<- her_pred$fit[her_pred$Biomass>B_pgy_l][1]
hBO<- her_pred$fit[her_pred$Biomass>B0_median][1]
hcurrentB<- her_pred$fit[her_pred$Biomass>currentB][1]
#percenatges
100-((hBMSY*100)/hBO)
(100-(hPGY_l*100)/hBO)-(100-((hBMSY*100)/hBO))
(100-(hBMSY*100)/hBO)-(100-(hPGY*100)/hBO)

#trade-offs figure matching with the distributions
a<- ggplot(NULL)+geom_line(aes(y=popgrowth_full_open2$median[!is.na(popgrowth_full_open2$median)], x=popgrowth_full_open2$Biomass[!is.na(popgrowth_full_open2$median)]), col="black", lwd=2)+
  geom_ribbon(aes(x=popgrowth_full_open$Biomass,ymin=popgrowth_full_open$`05%`[!is.na(popgrowth_full_open$`05%`)],ymax=popgrowth_full_open$`95%`[!is.na(popgrowth_full_open$`95%`)]), fill="grey", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ylim(c(0,20))+xlim(c(0,200))
b<- ggplot(size_pred[size_pred$Biomass<200,])+
  geom_line(aes(y=fit, x=Biomass), col='darkgreen',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='green', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("mean Size")+  coord_cartesian(xlim=c(0,200),ylim=c(10,30))+scale_y_continuous(position = "right")
c<- ggplot(tp_pred[tp_pred$Biomass<200,])+
  geom_line(aes(y=fit, x=Biomass), col='navyblue',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='blue', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Prb. top predators")+  coord_cartesian(xlim=c(0,200),ylim=c(0,1))+scale_y_continuous(position = "right")
d<- ggplot(tr_pred[tr_pred$Biomass<200,])+
  geom_line(aes(y=fit, x=Biomass), col='goldenrod',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='yellow', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Total species richness")+  coord_cartesian(xlim=c(0,200),ylim=c(0,200))+scale_y_continuous(position = "right")
e<- ggplot(her_pred[her_pred$Biomass<200,])+
  geom_line(aes(y=fit, x=Biomass), col='violetred4',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='violetred3', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Parrotfish scraping potential (cm^2/min)")+  coord_cartesian(xlim=c(0,200),ylim=c(0,750))+scale_y_continuous(position = "right")

#second-part of figure 3
windows()
ggarrange(a,b,c,d,e,nrow=2,ncol=3)

#status and ecosystem metrics
alldata$site_status_cat<-ifelse(alldata$site_status>1, "B>BMMSY","B<BMMSY")
windows()
a<-ggplot(alldata)+geom_jitter(aes(x=Size_marg_meth,y=site_status_cat,fill=as.factor(status_management)),alpha=0.05,pch=21)+
  geom_violin(aes(x=Size_marg_meth,y=site_status_cat,fill=as.factor(status_management),alpha=as.factor(status_management)),draw_quantiles = 0.5)+
  scale_fill_manual(values=c("Openly fished"="green","Restricted"="darkgreen"))+guides(fill=F,alpha=F)+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.5))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("")+xlab("Mean length (cm)")
b<-ggplot(alldata)+geom_jitter(aes(x=totalsp_marg_meth,y=site_status_cat,fill=as.factor(status_management)),alpha=0.05,pch=21)+
  geom_violin(aes(x=totalsp_marg_meth,y=site_status_cat,fill=as.factor(status_management),alpha=as.factor(status_management)),draw_quantiles = 0.5)+
  scale_fill_manual(values=c("Openly fished"="yellow","Restricted"="goldenrod"))+guides(fill=F,alpha=F)+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.5))+theme(axis.text.y = element_blank(), panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("")+xlab("Fish species richness")+coord_cartesian(xlim=c(0,250))
c<-ggplot(alldata[!is.na(alldata$PA_tpmarg_meth),])+geom_jitter(aes(x=PA_tpmarg_meth,y=site_status_cat,fill=as.factor(status_management)),alpha=0.05,pch=21)+
  geom_violin(aes(x=PA_tpmarg_meth,y=site_status_cat,fill=as.factor(status_management),alpha=as.factor(status_management)),draw_quantiles = 0.5)+
  scale_fill_manual(values=c("Openly fished"="blue","Restricted"="navyblue"))+guides(fill=F,alpha=F)+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.5))+theme(axis.text.y = element_blank(), panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("")+xlab("Presence/absence top predators")+scale_x_continuous(breaks = seq(0, 1, by = 1))
d<-ggplot(alldata)+geom_jitter(aes(x=herb_marg_meth,y=site_status_cat,fill=as.factor(status_management)),alpha=0.05,pch=21)+
  geom_violin(aes(x=herb_marg_meth,y=site_status_cat,fill=as.factor(status_management),alpha=as.factor(status_management)),draw_quantiles = 0.5)+
  scale_fill_manual(values=c("Openly fished"="violetred2","Restricted"="violetred4"))+guides(fill=F,alpha=F)+
  scale_alpha_manual(values=c("Openly fished"=0.3,"Restricted"=0.5))+theme(axis.text.y = element_blank(), panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("")+xlab("Parrotfish scrapping potential (cm^2/min)")+coord_cartesian(xlim=c(0,750))

ggarrange(a,b,c,d,nrow=1,ncol=4,widths=c(1.3,1,1,1),labels=c("f","g","h","i"))


#######################################################################################################################################

#if we did not include remote: parameters dependent on priors (i.e., not enough infor in only reserve recovery to on the later part of curve)
Fit_full_gf_nre<- stan(file = "Full_gompertzfox_noremote_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
ggarrange(plot(Fit_full_gf,pars="log_B0")+xlim(c(4,6)),plot(Fit_full_gf_nre,pars="log_B0")+xlim(c(4,6)),ncol=1,nrow=2)
exp(median(rstan::extract(Fit_full_gf_nre,pars=c("log_B0"))$log_B0))
exp(median(rstan::extract(Fit_full_gf,pars=c("log_B0"))$log_B0))
sd(rstan::extract(Fit_full_gf_nre,pars=c("log_B0"))$log_B0)
sd(rstan::extract(Fit_full_gf,pars=c("log_B0"))$log_B0)
contraction_full_nre<-as.data.frame(cbind(round(c((var_prior_logB0-(sd(rstan::extract(Fit_full_gf_nre,pars=c("log_B0"))$log_B0)^2))/var_prior_logB0,(var_prior_logr-(sd(rstan::extract(Fit_full_gf_nre,pars=c("log_r"))$log_r)^2))/var_prior_logr,(var_prior_logbmin-(sd(rstan::extract(Fit_full_gf_nre,pars=c("log_bmin"))$log_bmin)^2))/var_prior_logbmin,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,1])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,2])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,3])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,4])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,5])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,6])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,7])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,8])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,9])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,10])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,11])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,12])^2))/var_prior_beta,(var_prior_beta-(sd(rstan::extract(Fit_full_gf_nre,pars=c("beta"))$beta[,13])^2))/var_prior_beta,(var_priorI_fished-(sd(rstan::extract(Fit_full_gf_nre,pars=c("I_fished"))$I_fished)^2))/var_priorI_fished,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf_nre,pars=c("sigma_e"))$sigma_e)^2))/var_prior_sigma,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf_nre,pars=c("sigma_f"))$sigma_f)^2))/var_prior_sigma,(var_prior_sigma-(sd(rstan::extract(Fit_full_gf_nre,pars=c("sigma_u"))$sigma_u)^2))/var_prior_sigma),2),c("log_B0","log_r","log_bmin","b_op","b_gft","b_hc","b_at","b_rhc","b_rhf","b_rhbr","b_cmpc","b_sa","b_rs","b_cmds","b_d","b_gr","I_fished", "sd_res","sd_fis","sd_u")))
colnames( contraction_full_nre)<-c("posterior contraction","parameter")
print(contraction_full_nre)


#######################################################################################################################################

#model without gravity and with exports parametrized in different ways

Fit_full_exp<- stan(file = "Full_gompertzfox_exportrate.stan", data = stanDat_full_country,iter=10000,warmup = 9000, chains = 4,control = list(adapt_delta = 0.9999))
Fit_full_exppropr<- stan(file = "Full_gompertzfox_dig.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))

#check divergences
check_divergences(Fit_full_exp)


#compare against the version trying to estimate export rate (but not identifiable)
full_exp_loglik<- extract_log_lik(Fit_full_exp, merge_chains = F)
r_eff_full_exp <- relative_eff(exp(full_exp_loglik)) 
loo_full_exp <- loo(full_exp_loglik,r_eff =r_eff_full_exp)
full_exppropr_loglik<- extract_log_lik(Fit_full_exppropr, merge_chains = F)
r_eff_full_exppropr <- relative_eff(exp(full_exppropr_loglik)) 
loo_full_exppropr <- loo(full_exppropr_loglik,r_eff =r_eff_full_exppropr)

comp2 <- loo_compare(loo_full, loo_full_exp,loo_full_exppropr)
print(comp2, simplify=F)

#as another step, we test different fixed percentages of density-independent growth exported: (0%, 5%,10%,15%,20%,25%,30%)
stanDat_full_country$p<-0
Fit_full_exppropddr_0<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.05
Fit_full_exppropddr_5<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.10
Fit_full_exppropddr_10<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.15
Fit_full_exppropddr_15<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.20
Fit_full_exppropddr_20<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.25
Fit_full_exppropddr_25<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
stanDat_full_country$p<-0.30
Fit_full_exppropddr_30<- stan(file = "Full_gompertzfox_dig_fixed.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))

full_exppropddr0_loglik<- extract_log_lik(Fit_full_exppropddr_0, merge_chains = F)
r_eff_full_exppropddr0 <- relative_eff(exp(full_exppropddr0_loglik)) 
loo_full_exppropddr0 <- loo(full_exppropddr0_loglik,r_eff =r_eff_full_exppropddr0)

full_exppropddr5_loglik<- extract_log_lik(Fit_full_exppropddr_5, merge_chains = F)
r_eff_full_exppropddr5 <- relative_eff(exp(full_exppropddr5_loglik)) 
loo_full_exppropddr5 <- loo(full_exppropddr5_loglik,r_eff =r_eff_full_exppropddr5)

full_exppropddr10_loglik<- extract_log_lik(Fit_full_exppropddr_10, merge_chains = F)
r_eff_full_exppropddr10 <- relative_eff(exp(full_exppropddr10_loglik)) 
loo_full_exppropddr10 <- loo(full_exppropddr10_loglik,r_eff =r_eff_full_exppropddr10)

full_exppropddr15_loglik<- extract_log_lik(Fit_full_exppropddr_15, merge_chains = F)
r_eff_full_exppropddr15 <- relative_eff(exp(full_exppropddr15_loglik)) 
loo_full_exppropddr15 <- loo(full_exppropddr15_loglik,r_eff =r_eff_full_exppropddr15)

full_exppropddr20_loglik<- extract_log_lik(Fit_full_exppropddr_20, merge_chains = F)
r_eff_full_exppropddr20 <- relative_eff(exp(full_exppropddr20_loglik)) 
loo_full_exppropddr20 <- loo(full_exppropddr20_loglik,r_eff =r_eff_full_exppropddr20)

full_exppropddr25_loglik<- extract_log_lik(Fit_full_exppropddr_25, merge_chains = F)
r_eff_full_exppropddr25 <- relative_eff(exp(full_exppropddr25_loglik)) 
loo_full_exppropddr25 <- loo(full_exppropddr25_loglik,r_eff =r_eff_full_exppropddr25)

full_exppropddr30_loglik<- extract_log_lik(Fit_full_exppropddr_30, merge_chains = F)
r_eff_full_exppropddr30 <- relative_eff(exp(full_exppropddr30_loglik)) 
loo_full_exppropddr30 <- loo(full_exppropddr30_loglik,r_eff =r_eff_full_exppropddr30)

comp3 <- loo_compare(loo_full, loo_full_exppropddr0, loo_full_exppropddr5, loo_full_exppropddr10, loo_full_exppropddr15, loo_full_exppropddr20, loo_full_exppropddr25, loo_full_exppropddr30)
print(comp3, simplify=F)

################################################################################
#Choice of surplus production model

#originally we seeked to fit the generalized pella-tomlinson model
Fit_full_pt<- stan(file = "Full_pellatomlinson_free_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.9))
#pella tomllinson  with an additional parameter that defines where the surplus curve peaks (not enough info to estimate it)
#I always get divergences
#so I fit special cases of the pella-tomlinson

#gompertz fox (main preferred)
median(rstan::extract(Fit_full_gf,pars=c("MMSY"))$MMSY)
median(rstan::extract(Fit_full_gf,pars=c("BMMSY"))$BMMSY)

#graham scahefer
Fit_full_gs<- stan(file = "Full_schaefer_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.999))
median(rstan::extract(Fit_full_gs,pars=c("MMSY"))$MMSY)
median(rstan::extract(Fit_full_gs,pars=c("BMMSY"))$BMMSY)

stanDat_full_country$pe<-3
Fit_full_pt3<- stan(file = "Full_pellatomlinson_fixed_grav.stan", data = stanDat_full_country, chains = 4,control = list(adapt_delta = 0.99))
median(rstan::extract(Fit_full_pt3,pars=c("MMSY"))$MMSY)
median(rstan::extract(Fit_full_pt3,pars=c("BMMSY"))$BMMSY)

stanDat_full_country$pe<-4
Fit_full_pt4<- stan(file = "Full_pellatomlinson_fixed_grav.stan", data = stanDat_full_country,chains = 4,control = list(adapt_delta = 0.99))
median(rstan::extract(Fit_full_pt4,pars=c("MMSY"))$MMSY)
median(rstan::extract(Fit_full_pt4,pars=c("BMMSY"))$BMMSY)

#model selection
full_gs_loglik<- extract_log_lik(Fit_full_gs, merge_chains = F)
r_eff_full_gs <- relative_eff(exp(full_gs_loglik)) 
loo_full_gs <- loo(full_gs_loglik,r_eff =r_eff_full_gs)

full_pt3_loglik<- extract_log_lik(Fit_full_pt3, merge_chains = F)
r_eff_full_pt3 <- relative_eff(exp(full_pt3_loglik)) 
loo_full_pt3 <- loo(full_pt3_loglik,r_eff =r_eff_full_pt3)

full_pt4_loglik<- extract_log_lik(Fit_full_pt4, merge_chains = F)
r_eff_full_pt4 <- relative_eff(exp(full_pt4_loglik)) 
loo_full_pt4 <- loo(full_pt4_loglik,r_eff =r_eff_full_pt4)

comp_surplus <- loo_compare( loo_full,loo_full_gs,loo_full_pt3,loo_full_pt4)
rownames<-row.names(comp_surplus)
rownames2<-ifelse(rownames=="model1","Gompertz-Fox",ifelse(rownames=="model2","Graham-Schaefer",ifelse(rownames=="model3","Pella-Tomlinson 3","Pella-Tomlinson 4")))
row.names(comp_surplus)<-rownames2
print(comp_surplus, simplify=F)

#check those few with high pareto
highparetot<-unique(c(pareto_k_ids(loo_full, threshold = 0.7) , pareto_k_ids(loo_full_gs, threshold = 0.7),pareto_k_ids(loo_full_pt3, threshold = 0.7), pareto_k_ids(loo_full_pt4, threshold = 0.7)))

full_loglik<- extract_log_lik(Fit_full_gf, merge_chains = F)[,,-highparetot]
r_eff_full <- relative_eff(exp(full_loglik)) 
loo_full <- loo(full_loglik,r_eff =r_eff_full)

full_gs_loglik<- extract_log_lik(Fit_full_gs, merge_chains = F)[,,-highparetot]
r_eff_full_gs <- relative_eff(exp(full_gs_loglik)) 
loo_full_gs <- loo(full_gs_loglik,r_eff =r_eff_full_gs)

full_pt3_loglik<- extract_log_lik(Fit_full_pt3, merge_chains = F)[,,-highparetot]
r_eff_full_pt3 <- relative_eff(exp(full_pt3_loglik)) 
loo_full_pt3 <- loo(full_pt3_loglik,r_eff =r_eff_full_pt3)

full_pt4_loglik<- extract_log_lik(Fit_full_pt4, merge_chains = F)[,,-highparetot]
r_eff_full_pt4 <- relative_eff(exp(full_pt4_loglik)) 
loo_full_pt4 <- loo(full_pt4_loglik,r_eff =r_eff_full_pt4)

comp_surplus <- loo_compare( loo_full,loo_full_gs,loo_full_pt3,loo_full_pt4)
rownames<-row.names(comp_surplus)
rownames2<-ifelse(rownames=="model1","Gompertz-Fox",ifelse(rownames=="model2","Graham-Schaefer",ifelse(rownames=="model3","Pella-Tomlinson 3","Pella-Tomlinson 4")))
row.names(comp_surplus)<-rownames2
print(comp_surplus, simplify=F)
#write.csv(as.data.frame(print(comp_surplus, simplify=F)),"surplusmodel_comp.csv")


#extract posteriors
list_of_draws_gs <- as.data.frame(Fit_full_gs) 
list_of_draws_pt3 <- as.data.frame(Fit_full_pt3) 
list_of_draws_pt4 <- as.data.frame(Fit_full_pt4) 
#effect sizes consistent
windows()
ggarrange(plot(Fit_full_ss,pars="beta"),plot(Fit_full_gs,pars="beta"),plot(Fit_full_pt3,pars="beta"),plot(Fit_full_pt4,pars="beta"))

#extract model posterior
everything_gs<- rstan::extract(Fit_full_gs)
everything_pt3<- rstan::extract(Fit_full_pt3)
everything_pt4<- rstan::extract(Fit_full_pt4)

#surplus for each model
#sites_surplus_gs<-list_of_draws_gs$r* everything_gs$site_Bmarg * (1-(everything_gs$site_Bmarg/everything_gs$site_B0))
#sites_surplus_gs[sites_surplus_gs<0] <- 0


#MMSY and BMMSY at average environmental conditions
median(everything$MMSY)
median(everything_gs$MMSY)
median(everything_pt3$MMSY)
median(everything_pt4$MMSY)

median(everything$BMMSY)
median(everything_gs$BMMSY)
median(everything_pt3$BMMSY)
median(everything_pt4$BMMSY)

#status for each reef site (note this includes unfished sites as well, but we are doing status of fished reefs)
ordereddata$site_status_gs<- matrixStats::colMedians(everything_gs$site_status)
ordereddata$site_status_pt3<- matrixStats::colMedians(everything_pt3$site_status)
ordereddata$site_status_pt4<- matrixStats::colMedians(everything_pt4$site_status)
ordereddata$site_MMSY_gs<- matrixStats::colMedians(everything_gs$site_MMSY)
ordereddata$site_MMSY_pt3<- matrixStats::colMedians(everything_pt3$site_MMSY)
ordereddata$site_MMSY_pt4<- matrixStats::colMedians(everything_pt4$site_MMSY)
ordereddata$site_BMMSY_gs<- matrixStats::colMedians(everything_gs$site_BMMSY)
ordereddata$site_BMMSY_pt3<- matrixStats::colMedians(everything_pt3$site_BMMSY)
ordereddata$site_BMMSY_pt4<- matrixStats::colMedians(everything_pt4$site_BMMSY)

#ranges of MMSY
summary(ordereddata$site_MMSY)
summary(ordereddata$site_MMSY_gs)
summary(ordereddata$site_MMSY_pt3)
summary(ordereddata$site_MMSY_pt4)

#ranges of BMMSY
summary(ordereddata$site_BMMSY)
summary(ordereddata$site_BMMSY_gs)
summary(ordereddata$site_BMMSY_pt3)
summary(ordereddata$site_BMMSY_pt4)

#Status of exploited sites
ordereddata$site_catchstatusMMSY_gs<-ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY_gs
ordereddata$site_catchstatusMMSY_pt3<-ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY_pt3
ordereddata$site_catchstatusMMSY_pt4<-ordereddata$mean_tonnes_km2cell/ordereddata$site_MMSY_pt4
extracteddata2<-ordereddata[ordereddata$model_component==3,]

(length(extracteddata2$site_status_gs[extracteddata2$site_status_gs<1])/length(extracteddata2$site_status_gs))*100
(length(extracteddata2$site_status[extracteddata2$site_status<1])/length(extracteddata2$site_status))*100
(length(extracteddata2$site_status_pt3[extracteddata2$site_status_pt3<1])/length(extracteddata2$site_status_pt3))*100
(length(extracteddata2$site_status_pt4[extracteddata2$site_status_pt4<1])/length(extracteddata2$site_status_pt4))*100

(length(extracteddata2$site_catchstatusMMSY_gs[extracteddata2$site_catchstatusMMSY_gs<1])/length(extracteddata2$site_catchstatusMMSY_gs))*100
(length(extracteddata2$site_catchstatusMMSY[extracteddata2$site_catchstatusMMSY<1])/length(extracteddata2$site_catchstatusMMSY))*100
(length(extracteddata2$site_catchstatusMMSY_pt3[extracteddata2$site_catchstatusMMSY_pt3<1])/length(extracteddata2$site_catchstatusMMSY_pt3))*100
(length(extracteddata2$site_catchstatusMMSY_pt4[extracteddata2$site_catchstatusMMSY_pt4<1])/length(extracteddata2$site_catchstatusMMSY_pt4))*100

#conservation concern
(length(extracteddata2$site_catchstatusMMSY[extracteddata2$site_catchstatusMMSY>1| extracteddata2$site_status<1 ])/length(extracteddata2$site_catchstatusMMSY))*100
(length(extracteddata2$site_catchstatusMMSY_gs[extracteddata2$site_catchstatusMMSY_gs>1| extracteddata2$site_status_gs<1 ])/length(extracteddata2$site_catchstatusMMSY_gs))*100
(length(extracteddata2$site_catchstatusMMSY_pt3[extracteddata2$site_catchstatusMMSY_pt3>1| extracteddata2$site_status_pt3<1 ])/length(extracteddata2$site_catchstatusMMSY_pt3))*100
(length(extracteddata2$site_catchstatusMMSY_pt4[extracteddata2$site_catchstatusMMSY_pt4>1| extracteddata2$site_status_pt4<1 ])/length(extracteddata2$site_catchstatusMMSY_pt4))*100


#difference in parameter distributions
MMSYdata<-as.data.frame(cbind(c(list_of_draws_full$B0,list_of_draws_gs$B0,list_of_draws_pt3$B0,list_of_draws_pt4$B0),c(list_of_draws_full$MMSY,list_of_draws_gs$MMSY,list_of_draws_pt3$MMSY,list_of_draws_pt4$MMSY), c(list_of_draws_full$BMMSY,list_of_draws_gs$BMMSY,list_of_draws_pt3$BMMSY,list_of_draws_pt4$BMMSY),c(rep("Gompertz-Fox", length.out=length(list_of_draws_full$MMSY)),rep("Graham-Schaefer", length.out=length(list_of_draws_gs$MMSY)),rep("Pella-Tomlinson 3", length.out=length(list_of_draws_pt3$MMSY)),rep("Pella-Tomlinson 4", length.out=length(list_of_draws_pt4$MMSY)))))
colnames(MMSYdata)<-c("B0","MMSY","BMMSY","model")
MMSYdata$MMSY<-as.numeric(as.character(MMSYdata$MMSY))
MMSYdata$BMMSY<-as.numeric(as.character(MMSYdata$BMMSY))
MMSYdata$B0<-as.numeric(as.character(MMSYdata$B0))

a<-ggplot(MMSYdata, aes(x=MMSY,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+xlab("MMSY (t/km2/y)")+ylab("Posterior density")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))+xlim(c(0,30))

b<-ggplot(MMSYdata, aes(x=BMMSY,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+
  theme_classic()+xlab("BMMSY (t/km2)")+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"),fill="Model")+guides(fill=FALSE)
c<-ggplot(MMSYdata, aes(x=B0,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+
  theme_classic()+xlab("B0 (t/km2)")+ylab("")+ labs(x =  expression ("B0 ("~t/km^2*")"),fill="Model")

refpoints_pt<-ggarrange(a,b,c, widths=c(1,1,1.5), nrow=1, ncol=3, labels=c("a","b","c"))


#posteriors
Fit_gs_ss_summary <- summary(Fit_full_gs,probs=c(0.1,0.25,0.5,0.75, 0.9))
output_Fit_gs_ss<- as.data.frame(Fit_gs_ss_summary$summary)
Fit_pt3_ss_summary <- summary(Fit_full_pt3,probs=c(0.1,0.25,0.5,0.75, 0.9))
output_Fit_pt3_ss<- as.data.frame(Fit_pt3_ss_summary$summary)
Fit_pt4_ss_summary <- summary(Fit_full_pt4,probs=c(0.1,0.25,0.5,0.75, 0.9))
output_Fit_pt4_ss<- as.data.frame(Fit_pt4_ss_summary$summary)

#effect sizes from covariates
betas$model<-rep("Gompertz-Fox",nrow(betas))

betas_gs<- output_Fit_gs_ss[1:13,c("50%","10%","90%")]
betas_gs$variable<- c("Ocean productivity", "SST", "Hard coral","Atoll","Crest", "Flat", "Backreef/lagoon","Point count","Sampling area","Reserve size", "Distance sampling","Depth","Gravity")
betas_gs$order<- c(1,4,2,3,6,7,8,9,10,11,12,5,13)
betas_gs$variable <- factor(betas_gs$variable, levels = betas_gs$variable[order(betas_gs$order)])
betas_gs$model<-rep("Graham-Schaefer",nrow(betas_gs))

betas_pt3<- output_Fit_pt3_ss[1:13,c("50%","10%","90%")]
betas_pt3$variable<- c("Ocean productivity", "SST", "Hard coral","Atoll","Crest", "Flat", "Backreef/lagoon","Point count","Sampling area","Reserve size", "Distance sampling","Depth","Gravity")
betas_pt3$order<- c(1,4,2,3,6,7,8,9,10,11,12,5,13)
betas_pt3$variable <- factor(betas_pt3$variable, levels = betas_pt3$variable[order(betas_pt3$order)])
betas_pt3$model<-rep("P-T 3",nrow(betas_pt3))

betas_pt4<- output_Fit_pt4_ss[1:13,c("50%","10%","90%")]
betas_pt4$variable<- c("Ocean productivity", "SST", "Hard coral","Atoll","Crest", "Flat", "Backreef/lagoon","Point count","Sampling area","Reserve size", "Distance sampling","Depth","Gravity")
betas_pt4$order<- c(1,4,2,3,6,7,8,9,10,11,12,5,13)
betas_pt4$variable <- factor(betas_pt4$variable, levels = betas_pt4$variable[order(betas_pt4$order)])
betas_pt4$model<-rep("P-T 4",nrow(betas_pt4))

betas_pt<- ggplot()+
  geom_pointrange(data=betas_pt3,aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`),size=0.7,col="turquoise",alpha=0.5)+
   geom_pointrange(data=betas_pt4,aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`),size=0.7,col="purple",alpha=0.5)+
  geom_pointrange(data=betas_gs,aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`),size=0.7,col="darkred",alpha=0.5)+
  geom_pointrange(data=betas,aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`),size=0.7,col="darkgreen",alpha=0.5)+
 theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("Effect size")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")

betas_pt<- ggplot(rbind(betas,betas_gs,betas_pt3,betas_pt4),aes(x=variable,y=`50%`,ymin=`10%`,ymax=`90%`))+
  geom_pointrange(aes(col=model),size=0.7,alpha=0.5,position = position_dodge(width = 1))+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+
  ylab("Effect size")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")

windows()
ggarrange(refpoints_pt,betas_pt,nrow=2,ncol=1,labels = c("","d"),heights = c(1.5,1))


#fit under different surplus production models
pred_full_gs<- c(matrixStats::colMedians(rstan::extract(Fit_full_gs,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_gs,pars=c("mu2"))$mu2),matrixStats::colMedians(rstan::extract(Fit_full_gs,pars=c("mu3"))$mu3))
resid_full_gs<- c(log(reserves_complete$FamBiomass_tkm2),log(remote_complete$FamBiomass_tkm2), log(fished_data$FamBiomass_tkm2))-pred_full_gs
a_full_open_gs<- ggplot(data=NULL,aes(x=pred_full_gs,y=resid_full_gs))+geom_point()+theme_classic()+ggtitle("Graham-Schaefer")+xlab("Fitted ")+ylab("Residuals ")
b_full_open_gs<- ggplot(NULL, aes(x = resid_full_gs)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("Graham-Schaefer")+ylab("Count")+xlab("Residuals ")

pred_full_pt3<- c(matrixStats::colMedians(rstan::extract(Fit_full_pt3,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_pt3,pars=c("mu2"))$mu2),matrixStats::colMedians(rstan::extract(Fit_full_pt3,pars=c("mu3"))$mu3))
resid_full_pt3<- c(log(reserves_complete$FamBiomass_tkm2),log(remote_complete$FamBiomass_tkm2), log(fished_data$FamBiomass_tkm2))-pred_full_pt3
a_full_open_pt3<- ggplot(data=NULL,aes(x=pred_full_pt3,y=resid_full_pt3))+geom_point()+theme_classic()+ggtitle("P-T 3")+xlab("Fitted ")+ylab("Residuals ")
b_full_open_pt3<- ggplot(NULL, aes(x = resid_full_pt3)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("P-T 3")+ylab("Count")+xlab("Residuals ")

pred_full_pt4<- c(matrixStats::colMedians(rstan::extract(Fit_full_pt4,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_pt4,pars=c("mu2"))$mu2),matrixStats::colMedians(rstan::extract(Fit_full_pt4,pars=c("mu3"))$mu3))
resid_full_pt4<- c(log(reserves_complete$FamBiomass_tkm2),log(remote_complete$FamBiomass_tkm2), log(fished_data$FamBiomass_tkm2))-pred_full_pt4
a_full_open_pt4<- ggplot(data=NULL,aes(x=pred_full_pt4,y=resid_full_pt4))+geom_point()+theme_classic()+ggtitle("P-T 4")+xlab("Fitted ")+ylab("Residuals ")
b_full_open_pt4<- ggplot(NULL, aes(x = resid_full_pt4)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("P-T 4")+ylab("Count")+xlab("Residuals ")
windows()
ggarrange(a_full_open,b_full_open,a_full_open_gs,b_full_open_gs,a_full_open_pt3,b_full_open_pt3,a_full_open_pt4,b_full_open_pt4,nrow=4,ncol=2,labels=c("a","","b","","c","","d",""))

#status per jurisdiction.........................................................
#GS
sites_B0_gs<- everything_gs$site_B0
colnames(sites_B0_gs)<-ordereddata$Larger
sites_B0_gs<-melt(sites_B0_gs)
sites_B0_gs$iter_country<-paste(sites_B0_gs$iterations,sites_B0_gs$Var.2,sep="::")
sites_B0_gs$model_component<-rep(ordereddata$model_component,each=4000)
gf_B0distributions<-ddply(sites_B0_gs,.(iter_country),summarise,country=Var.2[1],iter=iterations[1],meanlogB0=mean(log(value)),n_sites=length(Var.2))
fac2 <- with( gf_B0distributions, reorder(country, meanlogB0, median, order = TRUE, na.rm=T))
gf_B0distributions<- within(gf_B0distributions, 
                               country<- factor(country, 
                                                levels=levels(fac2)))
gf_B0distributions$B0<-exp(gf_B0distributions$meanlogB0)
gf_B0distributions$meanBMMSY<-gf_B0distributions$B0/2
gf_B0distributions_h<-reshape2::dcast(gf_B0distributions, iter ~ country, value.var="B0")
gf_MMSYdistributions<-melt((list_of_draws_gs$r*gf_B0distributions_h[,-1])/4)
gf_MMSYdistributions$iter<-rep(seq(1:4000),length(unique(gf_MMSYdistributions$variable)))
gf_MMSYdistributions$iter_country<-paste(gf_MMSYdistributions$iter,gf_MMSYdistributions$variable,sep="::")
gf_MMSYdistributions$meanMMSY<-gf_MMSYdistributions$value
sites_B_gs<- everything_gs$site_Bmarg
colnames(sites_B_gs)<-ordereddata$Larger
sites_B_gs<-melt(sites_B_gs)
colnames(sites_B_gs)<-c("iterations","jurisdiction","Bmarg")
sites_B_gs$iter_country<-paste(sites_B_gs$iterations,sites_B_gs$jurisdiction,sep="::")
sites_B_gs$model_component<-rep(ordereddata$model_component,each=4000)
sites_B_gs$B0<-sites_B0_gs$value
sites_B_gs$definedprotection<-rep(ordereddata$definedprotection ,each=4000)
sites_B_gs_extracted<-sites_B_gs[sites_B_gs$model_component==3,]
gf_Bdistributions<-ddply(sites_B_gs_extracted,.(iter_country),summarise,jurisdiction=jurisdiction[1],iterations=iterations[1],meanBmarg=exp(mean(log(Bmarg))))
gf_Bdistributions<- merge(gf_Bdistributions,jurisdictionscale_data2[,c("Larger","mpa_perc")],by.x="jurisdiction",by.y="Larger",all.x=T)
gf_Bdistributions<-merge(gf_Bdistributions,gf_B0distributions[c("iter_country","meanBMMSY","B0")],by="iter_country",all.x=T)
gf_Bdistributions<-merge(gf_Bdistributions,gf_MMSYdistributions[,c("iter_country","meanMMSY")],by="iter_country",all.x=T)
gf_Bdistributions$weightedB_B0<-gf_Bdistributions$meanBmarg*(1-(gf_Bdistributions$mpa_perc/100))+gf_Bdistributions$B0*(gf_Bdistributions$mpa_perc/100)
gf_Bdistributions$Bstatus_onlyfishedreefs<-gf_Bdistributions$meanBmarg/gf_Bdistributions$meanBMMSY
gf_Bdistributions$Bstatus_B0<-gf_Bdistributions$weightedB_B0/gf_Bdistributions$meanBMMSY
country_refpoints_gs<- ddply(gf_Bdistributions,.(jurisdiction),summarize, larger_MMSY_gs=median(meanMMSY,na.rm=T),larger_B0_gs=median(B0,na.rm=T),larger_BMMSY_gs=median(meanBMMSY,na.rm=T),larger_Bstatus_gs=median(Bstatus_onlyfishedreefs),larger_Bstatus_B0_gs=median(Bstatus_B0,na.rm=T), larger_weightedB_gs=median(weightedB_B0,na.rm=T), larger_Bmarg_gs=median(meanBmarg,na.rm=T))
jurisdictionscale_data_gs<- merge(jurisdictionscale_data2,country_refpoints_gs,by.x="Larger", by.y="jurisdiction",all.x=T )
jurisdictionscale_data_gs$MMSY_median_gs<- ifelse(is.na(jurisdictionscale_data_gs$larger_MMSY_gs),median(gf_Bdistributions$meanMMSY),jurisdictionscale_data_gs$larger_MMSY_gs)
jurisdictionscale_data_gs$BMMSY_median_gs<- ifelse(is.na(jurisdictionscale_data_gs$larger_BMMSY_gs),median(gf_Bdistributions$meanBMMSY),jurisdictionscale_data_gs$larger_BMMSY_gs)
jurisdictionscale_data_gs$B0_median_gs<- ifelse(is.na(jurisdictionscale_data_gs$larger_B0_gs),median(gf_Bdistributions$B0),jurisdictionscale_data_gs$larger_B0_gs)
jurisdictions_nob_gs<- jurisdictionscale_data_gs[is.na(jurisdictionscale_data_gs$larger_MMSY_gs),]
jurisdictions_nob_gs<- droplevels(jurisdictions_nob_gs)
country_proboverfishing_gs<- matrix(NA,nrow=length(gf_Bdistributions$meanMMSY),ncol=length(jurisdictions_nob_gs$Larger))
for (i in 1:length(jurisdictions_nob_gs$catch_tkm2)){
  country_proboverfishing_gs[,i]<- jurisdictions_nob_gs$catch_tkm2[i]/gf_Bdistributions$meanMMSY
}
colnames(country_proboverfishing_gs)<-jurisdictions_nob_gs$area_name
country_proboverfishing_gs<-melt(country_proboverfishing_gs)
colnames(country_proboverfishing_gs)<-c("iterations","jurisdiction","fishingstatus")
gf_Bdistributions<-merge(gf_Bdistributions,jurisdictionscale_data2[,c("Larger","catch_tkm2")],by.x="jurisdiction",by.y="Larger",all.x=T)
gf_Bdistributions$fishingstatus<-gf_Bdistributions$catch_tkm2/gf_Bdistributions$meanMMSY
country_proboverfishing_gs<-rbind(country_proboverfishing_gs,gf_Bdistributions[,c("iterations","jurisdiction","fishingstatus")])
fac2 <- with( country_proboverfishing_gs, reorder(jurisdiction, fishingstatus, median, order = TRUE, na.rm=T))
country_proboverfishing_gsord<- within(country_proboverfishing_gs, 
                                      jurisdiction<- factor(jurisdiction, 
                                                            levels=levels(fac2)))
overfishing_country_gs<- ddply(country_proboverfishing_gsord,.(jurisdiction),summarize,medianFstatus_gs=median(fishingstatus,na.rm=T))
jurisdictionscale_data_gs<-merge(jurisdictionscale_data_gs,overfishing_country_gs,by.x="Larger2",by.y="jurisdiction",all.x=T)
jurisdictionscale_data_gs$overfishing_gs<- ifelse(is.na(jurisdictionscale_data_gs$medianFstatus_gs),NA,ifelse(jurisdictionscale_data_gs$medianFstatus_gs>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data_gs$overfishing_gs[!is.na(jurisdictionscale_data_gs$overfishing_gs)&jurisdictionscale_data_gs$overfishing_gs=="Overfishing"])/length(jurisdictionscale_data_gs$overfishing_gs[!is.na(jurisdictionscale_data_gs$overfishing_gs)])
jurisdictionscale_data_gs$belowBMMSY_gs<- ifelse(is.na(jurisdictionscale_data_gs$larger_Bstatus_gs),NA,ifelse(jurisdictionscale_data_gs$larger_Bstatus_gs>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_gs$belowBMMSY_gs[!is.na(jurisdictionscale_data_gs$belowBMMSY_gs)&jurisdictionscale_data_gs$belowBMMSY_gs=="belowBMMSY"])/length(jurisdictionscale_data_gs$belowBMMSY_gs[!is.na(jurisdictionscale_data_gs$belowBMMSY_gs)])
jurisdictionscale_data_gs$belowBMMSY_B0_gs<- ifelse(is.na(jurisdictionscale_data_gs$larger_Bstatus_B0_gs),NA,ifelse(jurisdictionscale_data_gs$larger_Bstatus_B0_gs>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_gs$belowBMMSY_B0_gs[!is.na(jurisdictionscale_data_gs$belowBMMSY_B0_gs)&jurisdictionscale_data_gs$belowBMMSY_B0_gs=="belowBMMSY"])/length(jurisdictionscale_data_gs$belowBMMSY_B0_gs[!is.na(jurisdictionscale_data_gs$belowBMMSY_B0_gs)])
jurisdictionscale_data_gs$concervationconcern<-ifelse(is.na(jurisdictionscale_data_gs$belowBMMSY_B0_gs)|is.na(jurisdictionscale_data_gs$overfishing_gs),NA, ifelse(jurisdictionscale_data_gs$belowBMMSY_B0_gs=="Not belowBMMSY" &jurisdictionscale_data_gs$overfishing_gs=="Not overfishing","In good condition","Conservation concern"))
length(jurisdictionscale_data_gs$concervationconcern[!is.na(jurisdictionscale_data_gs$concervationconcern)&jurisdictionscale_data_gs$concervationconcern=="Conservation concern"])/length(jurisdictionscale_data_gs$concervationconcern[!is.na(jurisdictionscale_data_gs$concervationconcern)])

#pt3
sites_B0_pt3<- everything_pt3$site_B0
colnames(sites_B0_pt3)<-ordereddata$Larger
sites_B0_pt3<-melt(sites_B0_pt3)
sites_B0_pt3$iter_country<-paste(sites_B0_pt3$iterations,sites_B0_pt3$Var.2,sep="::")
sites_B0_pt3$model_component<-rep(ordereddata$model_component,each=4000)
pt3_B0distributions<-ddply(sites_B0_pt3,.(iter_country),summarise,country=Var.2[1],iter=iterations[1],meanlogB0=mean(log(value)),n_sites=length(Var.2))
fac2 <- with( pt3_B0distributions, reorder(country, meanlogB0, median, order = TRUE, na.rm=T))
pt3_B0distributions<- within(pt3_B0distributions, 
                            country<- factor(country, 
                                             levels=levels(fac2)))
pt3_B0distributions$B0<-exp(pt3_B0distributions$meanlogB0)
pt3_B0distributions$meanBMMSY<-3^(1/(1-3))*pt3_B0distributions$B0
pt3_B0distributions_h<-reshape2::dcast(pt3_B0distributions, iter ~ country, value.var="B0")
pt3_MMSYdistributions<-melt((list_of_draws_pt3$r*pt3_B0distributions_h[,-1])*(1/(1+(3/2)))^((1/(3/2))+1))
pt3_MMSYdistributions$iter<-rep(seq(1:4000),length(unique(pt3_MMSYdistributions$variable)))
pt3_MMSYdistributions$iter_country<-paste(pt3_MMSYdistributions$iter,pt3_MMSYdistributions$variable,sep="::")
pt3_MMSYdistributions$meanMMSY<-pt3_MMSYdistributions$value
sites_B_pt3<- everything_pt3$site_Bmarg
colnames(sites_B_pt3)<-ordereddata$Larger
sites_B_pt3<-melt(sites_B_pt3)
colnames(sites_B_pt3)<-c("iterations","jurisdiction","Bmarg")
sites_B_pt3$iter_country<-paste(sites_B_pt3$iterations,sites_B_pt3$jurisdiction,sep="::")
sites_B_pt3$model_component<-rep(ordereddata$model_component,each=4000)
sites_B_pt3$B0<-sites_B0_pt3$value
sites_B_pt3$definedprotection<-rep(ordereddata$definedprotection ,each=4000)
sites_B_pt3_extracted<-sites_B_pt3[sites_B_pt3$model_component==3,]
pt3_Bdistributions<-ddply(sites_B_pt3_extracted,.(iter_country),summarise,jurisdiction=jurisdiction[1],iterations=iterations[1],meanBmarg=exp(mean(log(Bmarg))))
pt3_Bdistributions<- merge(pt3_Bdistributions,jurisdictionscale_data2[,c("Larger","mpa_perc")],by.x="jurisdiction",by.y="Larger",all.x=T)
pt3_Bdistributions<-merge(pt3_Bdistributions,pt3_B0distributions[c("iter_country","meanBMMSY","B0")],by="iter_country",all.x=T)
pt3_Bdistributions<-merge(pt3_Bdistributions,pt3_MMSYdistributions[,c("iter_country","meanMMSY")],by="iter_country",all.x=T)
pt3_Bdistributions$weightedB_B0<-pt3_Bdistributions$meanBmarg*(1-(pt3_Bdistributions$mpa_perc/100))+pt3_Bdistributions$B0*(pt3_Bdistributions$mpa_perc/100)
pt3_Bdistributions$Bstatus_onlyfishedreefs<-pt3_Bdistributions$meanBmarg/pt3_Bdistributions$meanBMMSY
pt3_Bdistributions$Bstatus_B0<-pt3_Bdistributions$weightedB_B0/pt3_Bdistributions$meanBMMSY
country_refpoints_pt3<- ddply(pt3_Bdistributions,.(jurisdiction),summarize, larger_MMSY_pt3=median(meanMMSY,na.rm=T),larger_B0_pt3=median(B0,na.rm=T),larger_BMMSY_pt3=median(meanBMMSY,na.rm=T),larger_Bstatus_pt3=median(Bstatus_onlyfishedreefs),larger_Bstatus_B0_pt3=median(Bstatus_B0,na.rm=T), larger_weightedB_pt3=median(weightedB_B0,na.rm=T), larger_Bmarg_pt3=median(meanBmarg,na.rm=T))
jurisdictionscale_data_pt3<- merge(jurisdictionscale_data2,country_refpoints_pt3,by.x="Larger", by.y="jurisdiction",all.x=T )
jurisdictionscale_data_pt3$MMSY_median_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$larger_MMSY_pt3),median(pt3_Bdistributions$meanMMSY),jurisdictionscale_data_pt3$larger_MMSY_pt3)
jurisdictionscale_data_pt3$BMMSY_median_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$larger_BMMSY_pt3),median(pt3_Bdistributions$meanBMMSY),jurisdictionscale_data_pt3$larger_BMMSY_pt3)
jurisdictionscale_data_pt3$B0_median_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$larger_B0_pt3),median(pt3_Bdistributions$B0),jurisdictionscale_data_pt3$larger_B0_pt3)
jurisdictions_nob_pt3<- jurisdictionscale_data_pt3[is.na(jurisdictionscale_data_pt3$larger_MMSY_pt3),]
jurisdictions_nob_pt3<- droplevels(jurisdictions_nob_pt3)
country_proboverfishing_pt3<- matrix(NA,nrow=length(pt3_Bdistributions$meanMMSY),ncol=length(jurisdictions_nob_pt3$Larger))
for (i in 1:length(jurisdictions_nob_pt3$catch_tkm2)){
  country_proboverfishing_pt3[,i]<- jurisdictions_nob_pt3$catch_tkm2[i]/pt3_Bdistributions$meanMMSY
}
colnames(country_proboverfishing_pt3)<-jurisdictions_nob_pt3$area_name
country_proboverfishing_pt3<-melt(country_proboverfishing_pt3)
colnames(country_proboverfishing_pt3)<-c("iterations","jurisdiction","fishingstatus")
pt3_Bdistributions<-merge(pt3_Bdistributions,jurisdictionscale_data2[,c("Larger","catch_tkm2")],by.x="jurisdiction",by.y="Larger",all.x=T)
pt3_Bdistributions$fishingstatus<-pt3_Bdistributions$catch_tkm2/pt3_Bdistributions$meanMMSY
country_proboverfishing_pt3<-rbind(country_proboverfishing_pt3,pt3_Bdistributions[,c("iterations","jurisdiction","fishingstatus")])
fac2 <- with( country_proboverfishing_pt3, reorder(jurisdiction, fishingstatus, median, order = TRUE, na.rm=T))
country_proboverfishing_pt3ord<- within(country_proboverfishing_pt3, 
                                       jurisdiction<- factor(jurisdiction, 
                                                             levels=levels(fac2)))
overfishing_country_pt3<- ddply(country_proboverfishing_pt3ord,.(jurisdiction),summarize,medianFstatus_pt3=median(fishingstatus,na.rm=T))
jurisdictionscale_data_pt3<-merge(jurisdictionscale_data_pt3,overfishing_country_pt3,by.x="Larger2",by.y="jurisdiction",all.x=T)
jurisdictionscale_data_pt3$overfishing_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$medianFstatus_pt3),NA,ifelse(jurisdictionscale_data_pt3$medianFstatus_pt3>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data_pt3$overfishing_pt3[!is.na(jurisdictionscale_data_pt3$overfishing_pt3)&jurisdictionscale_data_pt3$overfishing_pt3=="Overfishing"])/length(jurisdictionscale_data_pt3$overfishing_pt3[!is.na(jurisdictionscale_data_pt3$overfishing_pt3)])
jurisdictionscale_data_pt3$belowBMMSY_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$larger_Bstatus_pt3),NA,ifelse(jurisdictionscale_data_pt3$larger_Bstatus_pt3>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_pt3$belowBMMSY_pt3[!is.na(jurisdictionscale_data_pt3$belowBMMSY_pt3)&jurisdictionscale_data_pt3$belowBMMSY_pt3=="belowBMMSY"])/length(jurisdictionscale_data_pt3$belowBMMSY_pt3[!is.na(jurisdictionscale_data_pt3$belowBMMSY_pt3)])
jurisdictionscale_data_pt3$belowBMMSY_B0_pt3<- ifelse(is.na(jurisdictionscale_data_pt3$larger_Bstatus_B0_pt3),NA,ifelse(jurisdictionscale_data_pt3$larger_Bstatus_B0_pt3>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3[!is.na(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3)&jurisdictionscale_data_pt3$belowBMMSY_B0_pt3=="belowBMMSY"])/length(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3[!is.na(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3)])
jurisdictionscale_data_pt3$concervationconcern<-ifelse(is.na(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3)|is.na(jurisdictionscale_data_pt3$overfishing_pt3),NA, ifelse(jurisdictionscale_data_pt3$belowBMMSY_B0_pt3=="Not belowBMMSY" &jurisdictionscale_data_pt3$overfishing_pt3=="Not overfishing","In good condition","Conservation concern"))
length(jurisdictionscale_data_pt3$concervationconcern[!is.na(jurisdictionscale_data_pt3$concervationconcern)&jurisdictionscale_data_pt3$concervationconcern=="Conservation concern"])/length(jurisdictionscale_data_pt3$concervationconcern[!is.na(jurisdictionscale_data_pt3$concervationconcern)])

#pt4
sites_B0_pt4<- everything_pt4$site_B0
colnames(sites_B0_pt4)<-ordereddata$Larger
sites_B0_pt4<-melt(sites_B0_pt4)
sites_B0_pt4$iter_country<-paste(sites_B0_pt4$iterations,sites_B0_pt4$Var.2,sep="::")
sites_B0_pt4$model_component<-rep(ordereddata$model_component,each=4000)
pt4_B0distributions<-ddply(sites_B0_pt4,.(iter_country),summarise,country=Var.2[1],iter=iterations[1],meanlogB0=mean(log(value)),n_sites=length(Var.2))
fac2 <- with( pt4_B0distributions, reorder(country, meanlogB0, median, order = TRUE, na.rm=T))
pt4_B0distributions<- within(pt4_B0distributions, 
                             country<- factor(country, 
                                              levels=levels(fac2)))
pt4_B0distributions$B0<-exp(pt4_B0distributions$meanlogB0)
pt4_B0distributions$meanBMMSY<-4^(1/(1-4))*pt4_B0distributions$B0
pt4_B0distributions_h<-reshape2::dcast(pt4_B0distributions, iter ~ country, value.var="B0")
pt4_MMSYdistributions<-melt((list_of_draws_pt4$r*pt4_B0distributions_h[,-1])*(1/(1+(4/2)))^((1/(4/2))+1))
pt4_MMSYdistributions$iter<-rep(seq(1:4000),length(unique(pt4_MMSYdistributions$variable)))
pt4_MMSYdistributions$iter_country<-paste(pt4_MMSYdistributions$iter,pt4_MMSYdistributions$variable,sep="::")
pt4_MMSYdistributions$meanMMSY<-pt4_MMSYdistributions$value
sites_B_pt4<- everything_pt4$site_Bmarg
colnames(sites_B_pt4)<-ordereddata$Larger
sites_B_pt4<-melt(sites_B_pt4)
colnames(sites_B_pt4)<-c("iterations","jurisdiction","Bmarg")
sites_B_pt4$iter_country<-paste(sites_B_pt4$iterations,sites_B_pt4$jurisdiction,sep="::")
sites_B_pt4$model_component<-rep(ordereddata$model_component,each=4000)
sites_B_pt4$B0<-sites_B0_pt4$value
sites_B_pt4$definedprotection<-rep(ordereddata$definedprotection ,each=4000)
sites_B_pt4_extracted<-sites_B_pt4[sites_B_pt4$model_component==3,]
pt4_Bdistributions<-ddply(sites_B_pt4_extracted,.(iter_country),summarise,jurisdiction=jurisdiction[1],iterations=iterations[1],meanBmarg=exp(mean(log(Bmarg))))
pt4_Bdistributions<- merge(pt4_Bdistributions,jurisdictionscale_data2[,c("Larger","mpa_perc")],by.x="jurisdiction",by.y="Larger",all.x=T)
pt4_Bdistributions<-merge(pt4_Bdistributions,pt4_B0distributions[c("iter_country","meanBMMSY","B0")],by="iter_country",all.x=T)
pt4_Bdistributions<-merge(pt4_Bdistributions,pt4_MMSYdistributions[,c("iter_country","meanMMSY")],by="iter_country",all.x=T)
pt4_Bdistributions$weightedB_B0<-pt4_Bdistributions$meanBmarg*(1-(pt4_Bdistributions$mpa_perc/100))+pt4_Bdistributions$B0*(pt4_Bdistributions$mpa_perc/100)
pt4_Bdistributions$Bstatus_onlyfishedreefs<-pt4_Bdistributions$meanBmarg/pt4_Bdistributions$meanBMMSY
pt4_Bdistributions$Bstatus_B0<-pt4_Bdistributions$weightedB_B0/pt4_Bdistributions$meanBMMSY
country_refpoints_pt4<- ddply(pt4_Bdistributions,.(jurisdiction),summarize, larger_MMSY_pt4=median(meanMMSY,na.rm=T),larger_B0_pt4=median(B0,na.rm=T),larger_BMMSY_pt4=median(meanBMMSY,na.rm=T),larger_Bstatus_pt4=median(Bstatus_onlyfishedreefs),larger_Bstatus_B0_pt4=median(Bstatus_B0,na.rm=T), larger_weightedB_pt4=median(weightedB_B0,na.rm=T), larger_Bmarg_pt4=median(meanBmarg,na.rm=T))
jurisdictionscale_data_pt4<- merge(jurisdictionscale_data2,country_refpoints_pt4,by.x="Larger", by.y="jurisdiction",all.x=T )
jurisdictionscale_data_pt4$MMSY_median_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$larger_MMSY_pt4),median(pt4_Bdistributions$meanMMSY),jurisdictionscale_data_pt4$larger_MMSY_pt4)
jurisdictionscale_data_pt4$BMMSY_median_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$larger_BMMSY_pt4),median(pt4_Bdistributions$meanBMMSY),jurisdictionscale_data_pt4$larger_BMMSY_pt4)
jurisdictionscale_data_pt4$B0_median_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$larger_B0_pt4),median(pt4_Bdistributions$B0),jurisdictionscale_data_pt4$larger_B0_pt4)
jurisdictions_nob_pt4<- jurisdictionscale_data_pt4[is.na(jurisdictionscale_data_pt4$larger_MMSY_pt4),]
jurisdictions_nob_pt4<- droplevels(jurisdictions_nob_pt4)
country_proboverfishing_pt4<- matrix(NA,nrow=length(pt4_Bdistributions$meanMMSY),ncol=length(jurisdictions_nob_pt4$Larger))
for (i in 1:length(jurisdictions_nob_pt4$catch_tkm2)){
  country_proboverfishing_pt4[,i]<- jurisdictions_nob_pt4$catch_tkm2[i]/pt4_Bdistributions$meanMMSY
}
colnames(country_proboverfishing_pt4)<-jurisdictions_nob_pt4$area_name
country_proboverfishing_pt4<-melt(country_proboverfishing_pt4)
colnames(country_proboverfishing_pt4)<-c("iterations","jurisdiction","fishingstatus")
pt4_Bdistributions<-merge(pt4_Bdistributions,jurisdictionscale_data2[,c("Larger","catch_tkm2")],by.x="jurisdiction",by.y="Larger",all.x=T)
pt4_Bdistributions$fishingstatus<-pt4_Bdistributions$catch_tkm2/pt4_Bdistributions$meanMMSY
country_proboverfishing_pt4<-rbind(country_proboverfishing_pt4,pt4_Bdistributions[,c("iterations","jurisdiction","fishingstatus")])
fac2 <- with( country_proboverfishing_pt4, reorder(jurisdiction, fishingstatus, median, order = TRUE, na.rm=T))
country_proboverfishing_pt4ord<- within(country_proboverfishing_pt4, 
                                        jurisdiction<- factor(jurisdiction, 
                                                              levels=levels(fac2)))
overfishing_country_pt4<- ddply(country_proboverfishing_pt4ord,.(jurisdiction),summarize,medianFstatus_pt4=median(fishingstatus,na.rm=T))
jurisdictionscale_data_pt4<-merge(jurisdictionscale_data_pt4,overfishing_country_pt4,by.x="Larger2",by.y="jurisdiction",all.x=T)
jurisdictionscale_data_pt4$overfishing_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$medianFstatus_pt4),NA,ifelse(jurisdictionscale_data_pt4$medianFstatus_pt4>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data_pt4$overfishing_pt4[!is.na(jurisdictionscale_data_pt4$overfishing_pt4)&jurisdictionscale_data_pt4$overfishing_pt4=="Overfishing"])/length(jurisdictionscale_data_pt4$overfishing_pt4[!is.na(jurisdictionscale_data_pt4$overfishing_pt4)])
jurisdictionscale_data_pt4$belowBMMSY_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$larger_Bstatus_pt4),NA,ifelse(jurisdictionscale_data_pt4$larger_Bstatus_pt4>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_pt4$belowBMMSY_pt4[!is.na(jurisdictionscale_data_pt4$belowBMMSY_pt4)&jurisdictionscale_data_pt4$belowBMMSY_pt4=="belowBMMSY"])/length(jurisdictionscale_data_pt4$belowBMMSY_pt4[!is.na(jurisdictionscale_data_pt4$belowBMMSY_pt4)])
jurisdictionscale_data_pt4$belowBMMSY_B0_pt4<- ifelse(is.na(jurisdictionscale_data_pt4$larger_Bstatus_B0_pt4),NA,ifelse(jurisdictionscale_data_pt4$larger_Bstatus_B0_pt4>1,"Not belowBMMSY","belowBMMSY"))
length(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4[!is.na(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4)&jurisdictionscale_data_pt4$belowBMMSY_B0_pt4=="belowBMMSY"])/length(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4[!is.na(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4)])
jurisdictionscale_data_pt4$concervationconcern<-ifelse(is.na(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4)|is.na(jurisdictionscale_data_pt4$overfishing_pt4),NA, ifelse(jurisdictionscale_data_pt4$belowBMMSY_B0_pt4=="Not belowBMMSY" &jurisdictionscale_data_pt4$overfishing_pt4=="Not overfishing","In good condition","Conservation concern"))
length(jurisdictionscale_data_pt4$concervationconcern[!is.na(jurisdictionscale_data_pt4$concervationconcern)&jurisdictionscale_data_pt4$concervationconcern=="Conservation concern"])/length(jurisdictionscale_data_pt4$concervationconcern[!is.na(jurisdictionscale_data_pt4$concervationconcern)])

#example of surplus production comparison

popgrowth_gf<-rstan::extract(Fit_full_gf,pars=c("sustyield"))[[1]]
popgrowth_gf[is.nan(popgrowth_gf)] <- 0
popgrowth_gf<-broom.mixed::tidyMCMC(coda::as.mcmc(popgrowth_gf),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
popgrowth_gf$median=ifelse(popgrowth_gf$estimate<0,0.00000001,popgrowth_gf$estimate)
popgrowth_gf$Biomass=Bio
#places that should have no surplus production assign very small (to not get inf)
popgrowth_gf$`5%`<- ifelse(popgrowth_gf$conf.low<0,0.000001,popgrowth_gf$conf.low)
popgrowth_gf$`95%`<- ifelse(popgrowth_gf$conf.high<0,0.000001,popgrowth_gf$conf.high)

#for plotting purposes
popgrowth_full_open2<- popgrowth_full_open
popgrowth_full_open2$median<- ifelse(popgrowth_full_open2$median==0.00000001,NA,popgrowth_full_open2$median)
popgrowth_full_open2$`5%`<- ifelse(popgrowth_full_open2$`5%`==0.000001,NA,popgrowth_full_open2$`5%`)
popgrowth_full_open2$`95%`<- ifelse(popgrowth_full_open2$`95%`==0.000001,NA,popgrowth_full_open2$`95%`)
windows()
ggplot(NULL)+geom_line(aes(y=popgrowth_gf$median[!is.na(popgrowth_gf$median)], x=popgrowth_gf$Biomass[!is.na(popgrowth_gf$median)]), col="navyblue", lwd=2)+
  geom_line(aes(y=popgrowth_full_open$median[!is.na(popgrowth_full_open$median)], x=popgrowth_full_open$Biomass[!is.na(popgrowth_full_open$median)]), col="black", lwd=2)+
  geom_ribbon(aes(x=popgrowth_full_open$Biomass,ymin=popgrowth_full_open$`5%`,ymax=popgrowth_full_open2$`95%`), fill="black", alpha=0.3)+
  geom_ribbon(aes(x=popgrowth_gf$Biomass,ymin=popgrowth_gf$`5%`[!is.na(popgrowth_gf$`5%`)],ymax=popgrowth_gf$`95%`[!is.na(popgrowth_gf$`95%`)]), fill="blue", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))+ylim(c(0,20))+xlim(c(0,200))+
  geom_text(aes(x=150,y=10),label="Graham-Schaefer",col="black")+
  geom_text(aes(x=130,y=15),label="Gompertz-Fox",col="navyblue")

#######################################################################################################################################
#results with alternate catch statistics
catchdataoptions=jurisdictionscale_data2[,c("mean_tonnes",
                                            "mean_tonnes_reported",
                                            "mean_totalcatch_t","mean_tcatch_notindus_t","mean_tcatch_onlyreported_t", "mean_tcatch_notindus_onlyrep_t")]
#tnnes per km2 reef
catch_km2_options=catchdataoptions/jurisdictionscale_data2$CoralReefArea_km2

#fishing status
fishingstatus_catch=catch_km2_options/jurisdictionscale_data2$MMSY_median
pairs(~., data=log(fishingstatus_catch))
overfishing_catch=ifelse(is.na(fishingstatus_catch),NA,ifelse(fishingstatus_catch>1,"C>MMSY","C<=MMSY"))
#get mean fishing status per jurisdiction
fishingstatus_catch_mean<-fishingstatus_catch %>% select(mean_tonnes,mean_tonnes_reported:mean_tcatch_notindus_onlyrep_t)
fishingstatus_catch_mean$mean_fishingstatus<-rowMeans(fishingstatus_catch_mean,na.rm=T)
hist(fishingstatus_catch_mean$mean_fishingstatus)
fishingstatus_catch_mean$overfishing_catch=ifelse(is.na(fishingstatus_catch_mean$mean_fishingstatus),NA,ifelse(fishingstatus_catch_mean$mean_fishingstatus>1,"C>MMSY","C<=MMSY"))
summary(as.factor(fishingstatus_catch_mean$overfishing_catch))
length(fishingstatus_catch_mean$overfishing_catch[fishingstatus_catch_mean$overfishing_catch=="C>MMSY"])/length(fishingstatus_catch_mean$overfishing_catch)

#fishing status based on surplus curve
fishingstatus_curve_catch=catch_km2_options/jurisdictionscale_data2$surplus
overfishing_curve_catch=ifelse(is.na(fishingstatus_curve_catch),NA,ifelse(fishingstatus_curve_catch>1,1,0))
fishingstatus_curve_catch2=fishingstatus_curve_catch
fishingstatus_curve_catch2$biomassstatus=jurisdictionscale_data2$larger_Bstatus_B0
fisherystatus_catch=ifelse(is.na(fishingstatus_curve_catch),NA, ifelse(fishingstatus_curve_catch<1 & fishingstatus_curve_catch2$biomassstatus>1, "In good condition", ifelse(fishingstatus_curve_catch>1 & fishingstatus_curve_catch2$biomassstatus<1, "Unsustainable", ifelse(fishingstatus_curve_catch>1 & fishingstatus_curve_catch2$biomassstatus>1, "Warning", "Rebuilding"))))
#add quasi
fishingstatus_catch2=catch_km2_options/jurisdictionscale_data2$MMSY_median
fisherystatus_catch=ifelse(is.na(fisherystatus_catch),NA,ifelse(fisherystatus_catch=="Warning" &fishingstatus_catch2<1, "In good condition",as.character(fisherystatus_catch)))

#sumamry results
catch_status_sensitivity=matrix(NA,nrow=9,ncol=ncol(overfishing_catch))

for (i in 1:ncol(overfishing_catch)){
  catch_status_sensitivity[,i]=c(length(fishingstatus_catch[,i][!is.na(fishingstatus_catch[,i])]),
                                 (length(overfishing_catch[,i][!is.na(overfishing_catch[,i])&overfishing_catch[,i]=="C>MMSY"])/length(overfishing_catch[,i][!is.na(overfishing_catch[,i])]))*100,
                                 length(fishingstatus_curve_catch[,i][!is.na(fishingstatus_curve_catch[,i])]),
                                 (length(overfishing_curve_catch[,i][!is.na(overfishing_curve_catch[,i])&overfishing_curve_catch[,i]>0.9])/length(overfishing_curve_catch[,i][!is.na(overfishing_curve_catch[,i])]))*100,
                                 (length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="In good condition"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))*100,
                                 (length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Unsustainable"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))*100,
                                 (length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Warning"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))*100,
                                 (length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Rebuilding"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))*100,
                                 (length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&!fisherystatus_catch[,i]=="In good condition"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))*100)
}

catch_status_sensitivity=as.data.frame( catch_status_sensitivity)
colnames(catch_km2_options)
colnames( catch_status_sensitivity)=c("mean_spatial_reconstructed","mean_spatial_reported","non-spatial reconstructed","non-spatial non-industrial reconstructed","non-spatial reported","non-spatial reported non-industrial")
row.names(catch_status_sensitivity)=c("n_fishingstatus","perc_catching>MMSY","n_fisherystatus","perc_catching>surplus", "perc_sustainable","perc_unsustainable","perc_warning","perc_allowedrebuild","perc_conservation_concern")
#write.csv(catch_status_sensitivity,"sensitivity_catchdata.csv")


#######################################################################################################################################
#load(file='ZamborainMasonetal_2022_ReefSustainability_revised.RData')
#save.image(file='ZamborainMasonetal_2022_ReefSustainability_revised.RData')
