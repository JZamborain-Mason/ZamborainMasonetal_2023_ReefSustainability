#Code for Zamborain-Mason et al. 2020: Sustainability of the world's coral reef fisheries
#"R version 3.5.3 (2019-03-11)"

#This script implements the analyses in the manuscript Zamborain-Mason et al 2020

#First, it calculates MMSY and BMMSY for coral reef associated fish.  
#please note that results may vary slightly from the manuscript given that we are using a bayesian probabilistic framework


#Secondly, it uses catch (t/y) from the SAUP project, reef area (UNDP) and standing stock biomass (t/km2)  to estimate the status of the world's coral reefs from a production perpective at a jurisdiction scale.
# As observed reef fish biomass is obtained at a reef scale, for this section observed reef fish biomass is marginalized to account for methodological effects and we calculate the median marginalized biomass for a given jurisdiction. 

#Thirdly, we examine the jurisdiction level socio-economic factors that could be asscoiated with biomass and fishing status of a jurisdiction by fitting linear models to these relationships (and by also using generalized linear modelS).
# As not all factors are available for each country, we deal with missingness by predictve mean matching. 


#Next, we compare five ecosystem metrics (Total species richness, presence/absence of top predators, hard coral cover, parrotfish scraping potential and mean observed fish size) along a gradient of biomass to examine the potential trafe-offs between production and ecosystem objectives. 
#As with observed biomass, we also marginalize the ecosystem metrics by accounting for methodological effects

#Next, it perfoms the analyses from the shifting baseline box: How do lower unfished biomass impact our estimated surplus production. 

#Additionally, this script evaluates the representativennes of the our biomass samples with respect to its jurisdiction by looking at the relationship between our samples total gravity and their jurisditions mean total gravity. 

#Finally, I perform different sensitivity tests to explore how the robustness of our results change under differnet assumptions


#clear R
rm(list=ls())

#set working directory (where the data and code of this repository is stored)
setwd("c:/Users/jzamb/Documents/AUSTRALIA/PhD/Chapter 1/ALL ANALYSES/FINAL CODE AND DATA/FINAL AFTER REVISIONS (GITHUB)")


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


#to use or disconnect all cores
#rstan_options(auto_write = T)
#options(mc.cores = parallel::detectCores())
#backup_options <- options()
#options(backup_options)

#upload data (data is available in the paper's supporting information)
reefscale_data=read.csv("reefscale_data_submitted.csv", head=T)
jurisdictionscale_data=read.csv("jurisdictionscale_data_submitted.csv", head=T)

#functions used in this script
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

#functions to do model selection through k-fold cross-validation
#functions slightly modified from: https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd

stan_kfold <- function(file, list_of_datas, chains, cores,...){
  library(pbmcapply)
  badRhat <- 1.1 # don't know why we need this?
  n_fold <- length(list_of_datas)
  model <- stan_model(file=file)
  # First parallelize all chains:
  sflist <- 
    pbmclapply(1:(n_fold*chains), mc.cores = cores, 
               function(i){
                 # Fold number:
                 k <- ceiling(i / chains)
                 s <- sampling(model, data = list_of_datas[[k]], 
                               chains = 1, chain_id = i,...)
                 return(s)
               })
  
  # Then merge the K * chains to create K stanfits:
  stanfit <- list()
  for(k in 1:n_fold){
    inchains <- (chains*k - (chains - 1)):(chains*k)
    #  Merge `chains` of each fold
    stanfit[[k]] <- sflist2stanfit(sflist[inchains])
  }  
  return(stanfit) 
}

#extract log-likelihoods of held-out data
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  require(loo)
  K <- length(list_of_stanfits)
  list_of_log_liks <- plyr::llply(1:K, function(k){
    extract_log_lik(list_of_stanfits[[k]],...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  log_lik_heldout <- list_of_log_liks[[1]] * NA
  for(k in 1:K){
    log_lik <- list_of_log_liks[[k]]
    samples <- dim(log_lik)[1] 
    N_obs <- dim(log_lik)[2]
    # This is a matrix with the same size as log_lik_heldout
    # with 1 if the data was held out in the fold k
    heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
    # Sanity check that the previous log_lik is not being overwritten:
    if(any(!is.na(log_lik_heldout[heldout==1]))){
      warning("Heldout log_lik has been overwritten!!!!")
    }
    # We save here the log_lik of the fold k in the matrix:
    log_lik_heldout[heldout==1] <- log_lik[heldout==1]
  }
  return(log_lik_heldout)
}

#compute ELPD
kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise) #log predictive density
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  out <- list(
    pointwise = pointwise,
    elpd_kfold = elpd_kfold,
    se_elpd_kfold = se_elpd_kfold)
  #structure(out, class = "loo")  
  (out)
}

#modified from https://rdrr.io/cran/loo/src/R/loo_compare.R
#following Vehtari's comment: https://discourse.mc-stan.org/t/compare-models-with-k-fold-cv/9042/2
#kfold_a is the outout from kfold()
#pointwise expected log predictive density difference
elpd_diffs <- function(kfold_a, kfold_b) {
  pt_a <- kfold_a$pointwise
  pt_b <- kfold_b$pointwise
  elpd <- grep("^elpd", colnames(pt_a))
  pt_b[, elpd] - pt_a[, elpd]
}

#Compute standard error of the elpd difference
# diffs: Vector of pointwise elpd differences
se_elpd_diff <- function(diffs) {
  N <- length(diffs)
  sqrt(N) * sd(diffs)
}


#number of sites sampled per jurisidiction
sites=ddply(reefscale_data,.(Larger),summarize, reefsites=length(UniqueSite)) #bonaire is separated from neteherlands antilles but for reef area and catch it is the same
reefscale_data$Larger=as.factor(ifelse(reefscale_data$Larger=="Bonaire","Netherlands Antilles", as.character(reefscale_data$Larger)))
sites=ddply(reefscale_data,.(Larger),summarize, reefsites=length(UniqueSite)) #bonaire is separated from neteherlands antilles but for reef area and catch it is the same

#####
#SUSTAINABLE REFERENCE POINTS FOR CORAL REEF FISH...............................................................................................
#This section separates the unfished reef scale data (from high complainace marine reserves and remote uninhabited islands)
#then it checks for colinearity among explanatory variables
#it estimates sustainable reference points for coral reef fish


#standardize and relevel categorical variables
reefscale_data$indexr=ifelse(reefscale_data$Region=="C",1,ifelse(reefscale_data$Region=="I",2,3))
reefscale_data$Region=as.factor(reefscale_data$Region)
reefscale_data$Atoll=ifelse(reefscale_data$Atoll=="1",1,0)
reefscale_data=reefscale_data[!reefscale_data$ReefHabitat=="Bank",]
reefscale_data=droplevels(reefscale_data)
reefscale_data=reefscale_data %>%
  mutate(SampMethod=relevel(SampMethod, ref="Standard belt transect"),
         ReefHabitat=relevel(ReefHabitat,ref="Slope"),
         sClosure.size=standardise(log(Closure.size)),
         sSampArea=standardise(log(reefscale_data$SampArea)),
         sDepth=standardise(sqrt(Depth)),
         sOcean_prod=standardise(log(Prod_mgCm2d_finer)),
         sHardCoral=standardise(sqrt(HardCoral)))


#data potentially for benchmar (including 10 h baseline)
data_benchmark=reefscale_data[reefscale_data$section=="benchmark",]
data_benchmark=droplevels(data_benchmark)

#map remote and reserve locations used
newmap <- getMap(resolution = "high")

#jitter points 
reefscale_data$Site_Lat2=reefscale_data$Site_Lat+runif(length(reefscale_data$Site_Lat), min=0, max=3)
reefscale_data$Site_Long2=reefscale_data$Site_Long+runif(length(reefscale_data$Site_Long), min=0, max=4)
reefscale_data$Site_Lat2=ifelse(reefscale_data$Site_Lat2>23.5, reefscale_data$Site_Lat,reefscale_data$Site_Lat2)

#map of defined protection
reefscale_data$definedprotection=ifelse(reefscale_data$section=="benchmark"& reefscale_data$reserves==1,"HC Marine reserves", ifelse(reefscale_data$section=="benchmark"& reefscale_data$remote_20h==1 & reefscale_data$remote_20h_inhabited==0,"Remote",ifelse(reefscale_data$section=="status"& reefscale_data$Management=="Restricted","Restricted","Openly fished")))
windows()
ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=reefscale_data,aes(x=Site_Long2, y=Site_Lat2, fill = as.factor(reefscale_data$definedprotection)),colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="HC Marine reserves",],aes(x=Site_Long2, y=Site_Lat2), fill ="red",colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=reefscale_data[reefscale_data$definedprotection=="Remote",],aes(x=Site_Long2, y=Site_Lat2), fill ="turquoise1",colour="black", pch=21,size=3, alpha=0.7)+
  scale_fill_manual (name="Protection",values=c( "Openly fished"="darkorchid1","Restricted"="seagreen1","HC Marine reserves"="red", "Remote"="turquoise1"))+geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("Longitude")+
  ylab("Latitude")+theme_classic()  

#create dummy varibles for the different categorical variables
data_benchmark$pointintercept=ifelse(data_benchmark$SampMethod=="Point intercept", 1,0)
data_benchmark$backreef=ifelse(data_benchmark$ReefHabitat=="Lagoon_Back reef", 1,0)
data_benchmark$crest=ifelse(data_benchmark$ReefHabitat=="Crest", 1,0)
data_benchmark$flat=ifelse(data_benchmark$ReefHabitat=="Flat", 1,0)

#separate reserves and remote
reserves_complete=data_benchmark[data_benchmark$reserves==1,]
remote_complete=data_benchmark[data_benchmark$remote_20h==1 & data_benchmark$remote_20h_inhabited==0,]
reserves_complete=droplevels(reserves_complete)
remote_complete=droplevels(remote_complete)

#plot reserve and remote data
a=ggplot(reserves_complete, aes(x=Closure.age, y=FamBiomass_tkm2))+geom_point(aes(shape=reserves_complete$Region),col="darkgrey",alpha=0.5)+guides(shape=F)+theme_classic()+ labs(y = expression ("Biomass ("~t/km^2*")"))
b=ggplot(remote_complete, aes(x=Locality, y=log(FamBiomass_tkm2)))+geom_boxplot(fill="grey", col="black",alpha=0.5)+geom_jitter(width = 0.2, height = 0,col="darkgrey",alpha=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill="white",colour="black"))+ labs(y = expression ("Biomass (log("~t/km^2*"))"))
windows()
ggarrange(a,b, nrow=1, ncol=2)

#check correlation of continuous covariates
windows()
pairs(~SampMethod+sOcean_prod+
        sClosure.size+sSampArea+sDepth+sHardCoral+Closure.age,  data= reefscale_data,lower.panel=panel.cor )

windows()
pairs(~SampMethod+sOcean_prod+
        sClosure.size+sSampArea+sDepth+sHardCoral+Closure.age,  data= data_benchmark,lower.panel=panel.cor )
#reserve size and sampling are colinear (so we use one of them)

#non colinear continuous covariates
pairs(~sOcean_prod+sClosure.size+sDepth+sHardCoral+Closure.age,  data= data_benchmark,lower.panel=panel.cor )

#check correlation of continuous covariates with categorical ones
#sampling method, reef habitat and atoll
a=ggplot(data_benchmark,aes(x=as.factor(Atoll), y=sOcean_prod))+geom_boxplot()
b=ggplot(data_benchmark,aes(x=as.factor(ReefHabitat), y=sOcean_prod))+geom_boxplot()
c=ggplot(data_benchmark,aes(x=as.factor(SampMethod), y=sOcean_prod))+geom_boxplot()
d=ggplot(data_benchmark,aes(x=as.factor(Atoll), y=sHardCoral))+geom_boxplot()
e=ggplot(data_benchmark,aes(x=as.factor(ReefHabitat), y=sHardCoral))+geom_boxplot()
f=ggplot(data_benchmark,aes(x=as.factor(SampMethod), y=sHardCoral))+geom_boxplot()
g=ggplot(data_benchmark,aes(x=as.factor(Atoll), y=sDepth))+geom_boxplot()
h=ggplot(data_benchmark,aes(x=as.factor(ReefHabitat), y=sDepth))+geom_boxplot()
i=ggplot(data_benchmark,aes(x=as.factor(SampMethod), y=sDepth))+geom_boxplot()
j=ggplot(reserves_complete,aes(x=as.factor(Atoll), y=sClosure.size))+geom_boxplot()
k=ggplot(reserves_complete,aes(x=as.factor(ReefHabitat), y=sClosure.size))+geom_boxplot()
l=ggplot(reserves_complete,aes(x=as.factor(SampMethod), y=sClosure.size))+geom_boxplot()
m=ggplot(reserves_complete,aes(x=as.factor(Atoll), y=Closure.age))+geom_boxplot()
n=ggplot(reserves_complete,aes(x=as.factor(ReefHabitat), y=Closure.age))+geom_boxplot()
o=ggplot(reserves_complete,aes(x=as.factor(SampMethod), y=Closure.age))+geom_boxplot()
windows()
ggarrange(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,nrow=5,ncol=3)

#multicolineariy:
#reef habitat colinear with depth and size (e.g., slopes have deeper depths), 
#atoll colinear with ocean productivity (e.g.,sampled atolls have lower productivity)
#and sampling method colinear with size and coral cover (likely an artifact though)
#sampling area vs closure size (places with larger reserves, sampling area is usually larger)

#Full model we use coral cover, ocean productivity, depth and coral cover, age and size but care should be taken attributing the estimated effect to the specific covariates (it might be those they are colinear to)
#We take it out colinear covariates and keep the best options (crossvalidation we test different choices (e.g., Ocean productivity vs atoll))

#variance inflation factors
VIF.table=as.data.frame(vif(lm(log(FamBiomass_tkm2)~sOcean_prod+
                                 sClosure.size+sDepth+sHardCoral+Closure.age, data=data_benchmark)))
colnames(VIF.table)="VIF"
print(VIF.table)
#write.csv(VIF.table,"VIF.table.csv")

#create new data to plot model predictions
#As we have standardized all continuous variables, we create newdata  based only on reserve age (i.e., at average conditions)
newdata=with(reserves_complete, data.frame(Closure.age=seq(min(Closure.age), max(Closure.age), len=length(reserves_complete$Closure.age)))) 
Bio=seq(0,350,1) #biomass data for surplus production curve
B=length(Bio)

#null model: assuming closed populations that follow a schaeffer surplus-production model
stanDat_null <- list(p=0, b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                     N = nrow(reserves_complete),
                     br=log(remote_complete$FamBiomass_tkm2),
                     RE=nrow(remote_complete),
                     newage=newdata$Closure.age,
                     Bio=Bio, B=B)
Fit_null <- stan(file = "Null_exports_schaeffer.stan", data = stanDat_null, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))

#check model diagnostics
modeldiagplot1=stan_trace(Fit_null, pars=c("B0","bmin","r","sigma_e", "sigma_r"))
modeldiagplot2=stan_rhat(Fit_null)
pairs(Fit_null, pars=c("B0","bmin","r","sigma_e", "sigma_r","lp__"))
ggarrange(modeldiagplot1,modeldiagplot2, nrow=1,ncol=2, labels=c("a","b"), widths = c(1.8,1))

#check model fit
Fit_null_summary <- summary(Fit_null,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_null=as.data.frame(Fit_null_summary$summary)
#write.csv(output_Fit_null, "output_Fit_null.csv", row.names=T)

#examining model fit
row.names(output_Fit_null)
pred_null1=output_Fit_null[6:75,1]
resid_null1=log(reserves_complete$FamBiomass_tkm2)-pred_null1
pred_null2=output_Fit_null[76:155,1]
resid_null2=log(remote_complete$FamBiomass_tkm2)-pred_null2

a_null=ggplot(data=NULL,aes(x=pred_null1,y=resid_null1))+geom_point()+theme_classic()+ggtitle("reserves")+xlab("fitted ")+ylab("residuals ")
b_null=ggplot(NULL, aes(x = resid_null1)) +
  geom_histogram(colour = "white", fill = "black", bins=5) +theme_classic()+ggtitle("")+xlab("residuals ")
c_null=ggplot(data=NULL,aes(x=pred_null2,y=resid_null2))+geom_point()+theme_classic()+ggtitle("remote")+xlab("fitted ")+ylab("residuals ")
d_null=ggplot(NULL, aes(x = resid_null2)) +
  geom_histogram(colour = "white", fill = "black", bins=5) +theme_classic()+ggtitle("")+xlab("residuals ")

#posterior predictive checks
joined_sim <- rstan::extract(Fit_null)
n_sims <- length(joined_sim $lp__)
y_rep_reserves <- array(NA, c(n_sims, nrow(reserves_complete)))
y_rep_remote <- array(NA, c(n_sims, nrow(remote_complete)))

for (s in 1:n_sims){
  y_rep_reserves[s,] <- rnorm(nrow(reserves_complete), joined_sim$mu[s,], joined_sim$sigma_e[s])
  y_rep_remote[s,] <- rnorm(nrow(remote_complete), joined_sim$mu2[s,], joined_sim$sigma_r[s])
}
bayesplot::color_scheme_set(scheme = "gray")
a=bayesplot::ppc_dens_overlay(log(reserves_complete$FamBiomass_tkm2),y_rep_reserves[1:400,])+ggtitle("  ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
b=bayesplot::ppc_dens_overlay(log(remote_complete$FamBiomass_tkm2),y_rep_remote[1:400,])+ggtitle("")+ labs(y="Density",x = expression ("Biomass (log("~t/km^2*"))"))
windows()
ggarrange(a_null,c_null,b_null,d_null,a,b,nrow=3, ncol=2, labels=c("a","b","c","d","e","f"), heights=c(1,1,1.6))

#posterior draws
list_of_draws_null <- as.data.frame(Fit_null) #rows are iterations and columns are parameters

#show priors were not highly informative
sigma_e_prior=rcauchy(10000, 0,1)
sigma_r_prior=rcauchy(10000, 0,1)
BO_prior=rnorm(10000, 120,200)
bmin_prior=rnorm(10000, 0,200)
r_prior=rnorm(10000, 0.2,1)

length(list_of_draws_null$sigma_e)
a=ggplot(NULL)+geom_histogram(aes(x=list_of_draws_null$sigma_e),fill="blue",col="blue",alpha=0.5)+geom_histogram(aes(x=sigma_e_prior[sigma_e_prior>0]), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("sd reserve biomass (log)")+xlim(c(0,10))
b=ggplot(NULL)+geom_histogram(aes(x=list_of_draws_null$sigma_r),fill="blue",col="blue",alpha=0.5)+geom_histogram(aes(x=sigma_r_prior[sigma_r_prior>0]), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("sd remote biomass (log)")+xlim(c(0,10))
c=ggplot(NULL)+geom_histogram(aes(x=list_of_draws_null$B0),fill="blue",col="blue",alpha=0.5)+geom_histogram(aes(x=BO_prior[BO_prior>0]), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Unfished Biomass")
d=ggplot(NULL)+geom_histogram(aes(x=list_of_draws_null$bmin),fill="blue",col="blue",alpha=0.5)+geom_histogram(aes(x=bmin_prior[bmin_prior>0]), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Biomass reserve age 0")
e=ggplot(NULL)+geom_histogram(aes(x=list_of_draws_null$r),fill="blue",col="blue",alpha=0.5)+geom_histogram(aes(x=r_prior[r_prior>0]), lwd=2, fill="grey",alpha=0.5)+theme_classic()+ xlab("Intrinsic growth rate")
ggarrange(a,b,c,d,e,nrow=1,ncol=5, labels=c("a","b","c","d","e"))

#full model including covariates and spatial structure (i.e., locality as a random effect)
#add an index for locality (patially clustered sites)
remote_complete$indexl=as.numeric(remote_complete$Locality)
reserves_complete$indexl=as.numeric(reserves_complete$Locality)

##do the full version: with random effects and covariates
stanDat_full <- list(p=0,b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                     N = nrow(reserves_complete), d=reserves_complete$sDepth, hc=reserves_complete$sHardCoral, si=reserves_complete$sClosure.size,
                     op=reserves_complete$sOcean_prod,
                     R=nlevels(reserves_complete$Locality),
                     pr=reserves_complete$indexl,
                     br=log(remote_complete$FamBiomass_tkm2),d2=remote_complete$sDepth, hc2=remote_complete$sHardCoral, op2=remote_complete$sOcean_prod,
                     RE=nrow(remote_complete),
                     R2=nlevels(remote_complete$Locality),
                     pr2=remote_complete$indexl,
                     newage=newdata$Closure.age,
                     Bio=Bio, B=B)
#Fit_full_centered <- stan(file = "Full_exports_schaeffer.stan", data = stanDat_full,iter=100000,warmup=90000,chains = 4,control = list(adapt_delta = 0.999999999,stepsize = 0.0001,max_treedepth = 20))#decreased step size and increased adapt delta  to check for geometric ergodicity (divergences))
Fit_full <- stan(file = "Full_exports_schaeffer_noncentered.stan", data = stanDat_full,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9999,stepsize = 0.001,max_treedepth = 20))#non-centered version 

#check model diagnostics
stan_trace(Fit_full, pars=c("B0","bmin","r","sigma_e", "sigma_r"))
stan_rhat(Fit_full)
pairs(Fit_full, pars=c("B0","bmin","r","sigma_e", "sigma_r", "lp__"))
plot(Fit_full, pars=c("beta")) #covariates for ocean productivity, depth, coral cover and reserve size
plot(Fit_full, pars=c("u")) #random effects for reserve localities
plot(Fit_full, pars=c("u2")) #random effects for remote localities

#posterior draws
list_of_draws_full <- as.data.frame(Fit_full) #rows are iterations and columns are parameters

#difference in sustainable reference points between the null model, the full model, and the full model making bmin also a function of covariates
a=ggplot(NULL)+geom_density(aes(x=list_of_draws_null$MMSY),  alpha=0.5, lwd=2,lty=2)+geom_density(aes(x=list_of_draws_full$MMSY), fill="red", alpha=0.5)+theme_classic()+labs(x = expression ("MMSY ("~t/km^2/y*")"))+ylab("Posterior density")
b=ggplot(NULL)+geom_density(aes(x=list_of_draws_null$BMMSY),  alpha=0.5, lwd=2,lty=2)+geom_density(aes(x=list_of_draws_full$BMMSY), fill="red", alpha=0.5)+theme_classic()+labs(x = expression ("BMMSY ("~t/km^2*")"))+ylab("")
ggarrange(a,b, nrow=1, ncol=2, labels=c("a","b"), widths=c(1.1,1))

#model selection taking spatially clustered sites out at each time (defined by locality)
data_benchmark2=rbind(reserves_complete,remote_complete)
ldata=length(data_benchmark2$FamBiomass_tkm2)
locality=data_benchmark2$Locality
nlocality=length(unique(data_benchmark2$Locality))
hh <- kfold_split_grouped(nlocality, locality) 
holdout_data<- matrix(0, nrow = ldata, ncol = nlocality)
for(i in 1:ldata) holdout_data[i, hh[i]] <- 1

#match fold with locality
countryfold=as.data.frame(hh)
countryfold$Locality=data_benchmark2$Locality
countryfold$variable=paste("V",countryfold$hh,sep = "")

#now we separate it for reserves and remote locations
holdout_data2=as.data.frame(holdout_data)
holdout_data2$ManagementCategory=data_benchmark2$reserves
holdout_reserves=holdout_data2[holdout_data2$ManagementCategory==1,]
holdout_remote=holdout_data2[holdout_data2$ManagementCategory==0,]
holdout_remote$ManagementCategory=NULL
holdout_reserves$ManagementCategory=NULL
holdout_remote=as.matrix(holdout_remote)
holdout_reserves=as.matrix(holdout_reserves)

#turn into a list
holdout_reserves <- split(holdout_reserves,rep(1:ncol(holdout_reserves),each=nrow(holdout_reserves)))
holdout_remote <- split(holdout_remote ,rep(1:ncol(holdout_remote ),each=nrow(holdout_remote )))
holdout_data=split(holdout_data,rep(1:ncol(holdout_data),each=nrow(holdout_data)))

#null model holding out some data
stanDat_null_l <- rep(list(stanDat_null),nlocality)
#add the holdout index to it
for(i in 1:nlocality) stanDat_null_l[[i]]$holdout <- holdout_reserves[[i]]
for(i in 1:nlocality) stanDat_null_l[[i]]$holdout2 <- holdout_remote[[i]]

#full model holding out some data
stanDat_full_l= rep(list(stanDat_full),nlocality)
#add the holdout index to it
for(i in 1:nlocality) stanDat_full_l[[i]]$holdout <- holdout_reserves[[i]]
for(i in 1:nlocality) stanDat_full_l[[i]]$holdout2 <- holdout_remote[[i]]


#fit null model with cross validation
null_kcross <- stan_kfold(file="Null_exports_schaeffer_kfold.stan",stanDat_null_l ,chains=4,cores=2)
#extract log likelihood
loglik_null_kcross <- extract_log_lik_K(null_kcross,holdout_data)
#calculate the expected log pointwise predictive density elpd (i.e.,height (density) of the probability distribution, given the model parameters, at the data point (pointwise) that were held-out (predictive).
kcross_null <- kfold(loglik_null_kcross) 
elpd_null=kcross_null$elpd_kfold
elpd_se_null=kcross_null$se_elpd_kfold
logpreddensity_null=kcross_null$pointwise

#fit full model with cross validation
full_kcross <- stan_kfold(file="Full_exports_schaeffer_noncentered_kfold.stan",stanDat_full_l ,chains=4,cores=2)
loglik_full_kcross <- extract_log_lik_K(full_kcross,holdout_data)
kcross_full <- kfold(loglik_full_kcross) 
elpd_full=kcross_full$elpd_kfold
elpd_se_full=kcross_full$se_elpd_kfold
logpreddensity_full=kcross_full$pointwise


#compare the models (negative elpd_diff favors  the first model)
#poinwise difference in expected log predictive density: difference in their expected predictive accuracy
elpd_diff_null_full=sum(elpd_diffs(kcross_null,kcross_full))
se_elpd_diff_null_full=se_elpd_diff(elpd_diffs(kcross_null,kcross_full))
windows()
ggplot(NULL, aes(x=elpd_diffs(kcross_null,kcross_full)))+geom_density(fill="grey", alpha=0.5)+geom_vline(xintercept=0)+xlab("pointwise expected log predictive density difference (null minus full)")+
  geom_text(aes(x=-50,y=0.05,label=paste("sum elpd diff=",round(elpd_diff_null_full,2),"(SE:",round(se_elpd_diff_null_full,2),")")))+theme_classic()

#the null model perfomrs best in terms of preditive accurancy
#thus we proceed with that model

#obtain median 50% a¡and 90% credible intervals for parameters of interest
r_null=c(median(list_of_draws_null$r), mean(list_of_draws_null$r), quantile(list_of_draws_null$r,  probs=c(0.05,0.25,0.75, 0.95)))
B0_null=c(median(list_of_draws_null$B0), mean(list_of_draws_null$B0), quantile(list_of_draws_null$B0,  probs=c(0.05,0.25,0.75, 0.95)))
bmin_null=c(median(list_of_draws_null$bmin), mean(list_of_draws_null$bmin), quantile(list_of_draws_null$bmin,  probs=c(0.05,0.25,0.75, 0.95)))
bmmsy_null=c(median(list_of_draws_null$BMMSY), mean(list_of_draws_null$BMMSY), quantile(list_of_draws_null$BMMSY,  probs=c(0.05,0.25,0.75, 0.95)))
mmsy_null=c(median(list_of_draws_null$MMSY), mean(list_of_draws_null$MMSY), quantile(list_of_draws_null$MMSY,  probs=c(0.05,0.25,0.75, 0.95)))
ummsy_null=c(median(list_of_draws_null$uMMSY), mean(list_of_draws_null$uMMSY), quantile(list_of_draws_null$uMMSY,  probs=c(0.05,0.25,0.75, 0.95)))

Fit_null_pars=cbind(r_null,B0_null,bmin_null,bmmsy_null,mmsy_null,ummsy_null)
rownames(Fit_null_pars)=c("median","mean","5%","25%","75%","95%")
print(Fit_null_pars)
#write.csv(Fit_null_pars, "null_pars.csv", row.names=F)

#plot median reserve trajectory with 50 and 90% credible intervals
row.names(output_Fit_null)
newpred_null=output_Fit_null[159:228,]
colnames(newpred_null)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
newpred_null$Closure.age=newdata$Closure.age
null_tra=ggplot(NULL)+
  geom_point(data=reserves_complete, aes(y=FamBiomass_tkm2, x=Closure.age),col="darkgrey",fill="darkgrey", alpha=0.35)+
  geom_line(aes(y=newpred_null$median, x=newpred_null$Closure.age), col="navyblue",lwd=2)+geom_ribbon(aes(x=newpred_null$Closure.age,ymin=newpred_null$`5%`, ymax=newpred_null$`95%`), alpha=0.3, fill="blue")+geom_ribbon(aes(x=newpred_null$Closure.age,ymin=newpred_null$`25%`, ymax=newpred_null$`75%`), alpha=0.3, fill="navyblue")+theme_classic()+xlab("MPA age (years)")+ylab("Biomass (t/km2)")+ labs(y = expression ("Biomass ("~t/km^2*")"))

#plot posterior unfished biomass with 90% credible intervals on remote boxplot
null_remote=ggplot(remote_complete, aes(x=Locality, y=log(FamBiomass_tkm2)))+geom_boxplot(fill="darkgrey", alpha=0.35)+geom_jitter(width = 0.2, height = 0,col="darkgrey",alpha=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+geom_abline(intercept=log(Fit_null_pars[1,2]), slope=0, col="navyblue",lwd=2)+geom_abline(intercept=log(Fit_null_pars[3,2]),slope=0, lty=2, col="navyblue")+geom_abline(intercept=log(Fit_null_pars[6,2]),slope=0, lty=2, col="navyblue")+ labs(y = expression ("Biomass (log("~t/km^2*"))"))

#plot surplus production curve along a gradient of biomass
popgrowth_null=output_Fit_null[309:659,]
colnames(popgrowth_null)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
popgrowth_null$mean_estimate=ifelse(popgrowth_null$mean_estimate<0,0,popgrowth_null$mean_estimate)
popgrowth_null$median=ifelse(popgrowth_null$median<0,0.000001,popgrowth_null$median)
popgrowth_null$Biomass=Bio

#places that should have no surplus production
#assign very small (to not get inf)
popgrowth_null$`5%`=ifelse(popgrowth_null$`5%`<0,0.000001,popgrowth_null$`5%`)
popgrowth_null$`95%`=ifelse(popgrowth_null$`95%`<0,0.000001,popgrowth_null$`95%`)
popgrowth_null$`25%`=ifelse(popgrowth_null$`25%`<0,0.000001,popgrowth_null$`25%`)
popgrowth_null$`75%`=ifelse(popgrowth_null$`75%`<0,0.000001,popgrowth_null$`75%`)

#for plotting purposes
popgrowth_null2=popgrowth_null
popgrowth_null2$median=ifelse(popgrowth_null2$median==0.000001,NA,popgrowth_null2$median)
popgrowth_null2$`5%`=ifelse(popgrowth_null2$`5%`==0.000001,NA,popgrowth_null2$`5%`)
popgrowth_null2$`95%`=ifelse(popgrowth_null2$`95%`==0.000001,NA,popgrowth_null2$`95%`)
popgrowth_null2$`25%`=ifelse(popgrowth_null2$`25%`==0.000001,NA,popgrowth_null2$`25%`)
popgrowth_null2$`75%`=ifelse(popgrowth_null2$`75%`==0.000001,NA,popgrowth_null2$`75%`)

#shading regions
shade1<- rbind(c(0,0), subset(popgrowth_null2, Biomass >= 
                                0 & Biomass <= median(list_of_draws_null$BMMSY)), c(median(list_of_draws_null$BMMSY), 0))
shade3<- rbind(c(0,Inf), subset(popgrowth_null2, Biomass >= 
                                  0 & Biomass <= median(list_of_draws_null$BMMSY)), c(median(list_of_draws_null$BMMSY), Inf))
null_grow1=ggplot(NULL)+geom_polygon(data=shade1,aes(Biomass, median), fill="cyan3",alpha=0.7)+
  geom_rect(aes(xmin=median(list_of_draws_null$BMMSY), xmax=Inf,ymin=0,ymax=median(list_of_draws_null$MMSY)),fill="navyblue",alpha=0.6)+geom_polygon(data=shade3,aes(Biomass, median), fill="red3",alpha=0.7)+
  geom_rect(aes(xmin=median(list_of_draws_null$BMMSY), xmax=Inf,ymin=median(list_of_draws_null$MMSY)),ymax=Inf, fill="goldenrod1",alpha=0.7)+
  geom_line(aes(y=popgrowth_null2$median[!is.na(popgrowth_null2$median)], x=popgrowth_null2$Biomass[!is.na(popgrowth_null2$median)]), col="black", lwd=2)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression (atop("Potential sustainable", paste("yield ("~t/km^2/y*")"))))+geom_vline(xintercept=median(list_of_draws_null$BMMSY), lwd=1.2)+scale_x_continuous(expand = c(0, 0),limits=c(0,230)) + scale_y_continuous(expand = c(0, 0),limits=c(0,3))

null_grow=ggplot(NULL)+geom_line(aes(y=popgrowth_null2$median[!is.na(popgrowth_null2$median)], x=popgrowth_null2$Biomass[!is.na(popgrowth_null2$median)]), col="navyblue", lwd=2)+
  geom_ribbon(aes(x=popgrowth_null$Biomass,ymin=popgrowth_null$`5%`[!is.na(popgrowth_null$`5%`)],ymax=popgrowth_null$`95%`[!is.na(popgrowth_null$`95%`)]), fill="blue", alpha=0.3)+
  geom_ribbon(aes(x=popgrowth_null$Biomass,ymin=popgrowth_null$`25%`[!is.na(popgrowth_null$`25%`)],ymax=popgrowth_null$`75%`[!is.na(popgrowth_null$`75%`)]), fill="navyblue", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))+ylim(c(0,3.5))
modelfit1=ggarrange(null_tra,null_remote, null_grow, nrow=1, ncol=3, labels=c("a","b","c"))

#plot MMSY and BMMSY distribution
dmmsy_null=ggplot(list_of_draws_null)+geom_density(aes(x=MMSY), col="black",fill="darkgrey",alpha=0.7)+theme_classic()+xlab("MMSY (t/km2/y)")+ylab("Posterior density")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))#+geom_vline(xintercept=median(list_of_draws_null$MMSY), lwd=2)
dbmmsy_null=ggplot(list_of_draws_null)+geom_density(aes(x=BMMSY), col="black",fill="darkgrey",alpha=0.7)+theme_classic()+xlab("B_MMSY (t/km2)")+ylab("Posterior density")+ labs(x =  expression (B["MMSY "]*"("~t/km^2*")"))#+geom_vline(xintercept=median(list_of_draws_null$BMMSY),lwd=2)
fig1a=ggarrange(dmmsy_null,dbmmsy_null,null_grow1,nrow=3, ncol=1, labels=c("a","b","c"))

#unfished biomass, intrinsic growth rate and biomass at reserve age 0 posterior distributions
B0_null=ggplot(list_of_draws_null)+geom_histogram(aes(x=B0, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("Unfished Biomass(t/km2)")+ylab("")+ labs(x = expression ("Unfished Biomass ("~t/km^2*")"))
r_null=ggplot(list_of_draws_null)+geom_histogram(aes(x=r, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("Intrinsic growth rate")+ylab("Posterior density")
bmin_null=ggplot(list_of_draws_null)+geom_histogram(aes(x=bmin, y=..count../sum(..count..)), col="black",fill="navyblue",alpha=0.7)+theme_classic()+xlab("Biomass age 0 (t/km2)")+ylab("")+ labs(x = expression ("Biomass age 0 ("~t/km^2*")"))
logistic_par_post=ggarrange(r_null,B0_null,bmin_null,nrow=1, ncol=3, labels=c("d","e","f"))
modelfit=ggarrange(modelfit1, logistic_par_post,nrow=2, ncol=1, heights=c(1.5,1))
windows()
modelfit


#For the biomass and fishing status, we use the median and 90% credible intervals of the BMMSY and MMSY respectively to define places as oversihed or overfishing.
#we will also use the distributions of MMSY and BMMSY to show the probabilities of a jurisdiction being overfished/overfishing
#For the sustainability status we will use the surplus production curve median and 90% credible intervals.

#Thus, we get assign the estimated B0, MMSY and BMMSY values
MMSY_median=Fit_null_pars[1,5]
MMSY_low=Fit_null_pars[3,5]
MMSY_high=Fit_null_pars[6,5]
BMMSY_median=Fit_null_pars[1,4]
BMMSY_low=Fit_null_pars[3,4]
BMMSY_high=Fit_null_pars[6,4]
B0_median=Fit_null_pars[1,2]
B0_low=Fit_null_pars[3,2]
B0_high=Fit_null_pars[6,2]

#####
#STATUS OF CORAL REEF FISH STOCKS...............................................................................................

#Status of the world's coral reef fish fisheries from both a production and ecosystem perspective
#select only sites that are open to fishing in some way (i.e., not reserves)

alldata=reefscale_data[!(reefscale_data$Management=="UnfishedHigh"),]
length(alldata$FamBiomass_tkm2)

#Correct for methodological factors---------------------------------
#First we get the reef scale response variables and account for the potential effects from the way they were sampled
#we calculate marginalized response variables using "Slope", "Standard Belt transect", "4-10m" and mean sampling area as a reference. 

#correlation among methodological covariates
pairs(~DepthCategory+
        ReefHabitat+
        SampArea+
        SampMethod, data= alldata,lower.panel=panel.cor )

#relevel and standardize methodological factors
summary(alldata$DepthCategory)
alldata$DepthCategory<-relevel(alldata$DepthCategory,ref="4-10m")
summary(alldata$ReefHabitat)
alldata$ReefHabitat<-relevel(alldata$ReefHabitat,ref="Slope")
summary(alldata$SampMethod)
alldata$SampMethod<-relevel(alldata$SampMethod,ref="Standard belt transect")


#mixed effects models with methodological covaraites to calculate reef scale marginalized response variables
#BIOMASS
hist(log(alldata$FamBiomass_tkm2))
alldata$lBiomass<-log(alldata$FamBiomass_tkm2)
model_Biomass<-lmer(lBiomass~ 
                      DepthCategory+
                      ReefHabitat+
                      sSampArea+
                      SampMethod+
                      (1|Locality),data=alldata)
#model fit
residualsB<-resid(model_Biomass)
a=ggplot(NULL)+geom_histogram(aes(x=residualsB))+theme_classic()+xlab ("Residuals biomass model")
b=ggplot(NULL)+geom_point(aes(x=predict(model_Biomass),y=resid(model_Biomass)))+theme_classic()+xlab ("Fitted biomass (log)")+ylab ("Residuals")

#marginalized
alldata$lB_marg=alldata$lBiomass-(fixef(model_Biomass)["sSampArea"]*alldata$sSampArea+
                                    ifelse(alldata$ReefHabitat=="Crest",fixef(model_Biomass)[c("ReefHabitatCrest")],
                                           ifelse(alldata$ReefHabitat=="Flat",fixef(model_Biomass)[c("ReefHabitatFlat")],
                                                  ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                         fixef(model_Biomass)[c("ReefHabitatLagoon_Back reef")],0)))+
                                    ifelse(alldata$DepthCategory=="0-4m",fixef(model_Biomass)[c("DepthCategory0-4m")],
                                           ifelse(alldata$DepthCategory==">10m",fixef(model_Biomass)[c("DepthCategory>10m")],0))+
                                    ifelse(alldata$SampMethod=="Distance sampling",fixef(model_Biomass)[c("SampMethodDistance sampling")],
                                           ifelse(alldata$SampMethod=="Point intercept",fixef(model_Biomass)[c("SampMethodPoint intercept")],0)))
alldata$Biomass_marg=exp(alldata$lB_marg)

#OBSERVED FISH LENGTH
hist(log(alldata$meanSize_cm))
alldata$lSize<-log(alldata$meanSize_cm)
model_Size<-lmer(lSize~ 
                   DepthCategory+
                   ReefHabitat+
                   sSampArea+
                   SampMethod+
                   (1|Locality),data=alldata)
residualsS<-resid(model_Size)
c=ggplot(NULL)+geom_histogram(aes(x=residualsS))+theme_classic()+xlab ("Residuals mean length model")
d=ggplot(NULL)+geom_point(aes(x=predict(model_Size),y=resid(model_Size)))+theme_classic()+xlab ("Fitted mean length (log)")+ylab ("Residuals")
alldata$lS_marg=alldata$lSize-(fixef(model_Size)["sSampArea"]*alldata$sSampArea+
                                 ifelse(alldata$ReefHabitat=="Crest",fixef(model_Size)[c("ReefHabitatCrest")],
                                        ifelse(alldata$ReefHabitat=="Flat",fixef(model_Size)[c("ReefHabitatFlat")],
                                               ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                      fixef(model_Size)[c("ReefHabitatLagoon_Back reef")],0)))+
                                 ifelse(alldata$DepthCategory=="0-4m",fixef(model_Size)[c("DepthCategory0-4m")],
                                        ifelse(alldata$DepthCategory==">10m",fixef(model_Size)[c("DepthCategory>10m")],0))+
                                 ifelse(alldata$SampMethod=="Distance sampling",fixef(model_Size)[c("SampMethodDistance sampling")],
                                        ifelse(alldata$SampMethod=="Point intercept",fixef(model_Size)[c("SampMethodPoint intercept")],0)))
alldata$Size_marg=exp(alldata$lS_marg)

#PROBBAILITY OF ENCOUNTERING TOP PREDATORS
hist(log(alldata$TPBiomass_tkm2+1))
length(alldata$TPBiomass_tkm2[alldata$TPBiomass_tkm2==0])/length(alldata$TPBiomass_tkm2)
#for 67% of our fished sites observed top predator biomass was 0, thus we look at presence/absence (binomial)
reefscale_data$PA_tp=ifelse(reefscale_data$TPBiomass_tkm2>0,1,0)
alldata$PA_tp=ifelse(alldata$TPBiomass_tkm2>0,1,0)

model_tp<-glmer(PA_tp~ 
                  DepthCategory+
                  ReefHabitat+
                  sSampArea+
                  SampMethod+
                  (1|Locality),data=alldata,
                family=binomial(link = "logit"),
                control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000000)))
e=ggplot(NULL)+geom_histogram(aes(x=resid(model_tp)))+theme_classic()+ xlab ("Residuals top predator model")
f=ggplot(NULL)+geom_point(aes(x=predict(model_tp),y=resid(model_tp)))+theme_classic()+xlab ("Fitted")+ylab ("Residuals")
alldata$tp_marg=alldata$PA_tp-(fixef(model_tp)["sSampArea"]*alldata$sSampArea+
                                 ifelse(alldata$ReefHabitat=="Crest",fixef(model_tp)[c("ReefHabitatCrest")],
                                        ifelse(alldata$ReefHabitat=="Flat",fixef(model_tp)[c("ReefHabitatFlat")],
                                               ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                      fixef(model_tp)[c("ReefHabitatLagoon_Back reef")],0)))+
                                 ifelse(alldata$DepthCategory=="0-4m",fixef(model_tp)[c("DepthCategory0-4m")],
                                        ifelse(alldata$DepthCategory==">10m",fixef(model_tp)[c("DepthCategory>10m")],0))+
                                 ifelse(alldata$SampMethod=="Distance sampling",fixef(model_tp)[c("SampMethodDistance sampling")],
                                        ifelse(alldata$SampMethod=="Point intercept",fixef(model_tp)[c("SampMethodPoint intercept")],0)))
#change from log-odds to probability and then back to presence absence
alldata$prob_tpmarg=exp(alldata$tp_marg)/(1+exp(alldata$tp_marg))
alldata$PA_tpmarg=ifelse(alldata$prob_tpmarg>0.5,1,0)
length(alldata$PA_tpmarg[alldata$PA_tpmarg==0])/length(alldata$PA_tpmarg)

#CORAL COVER
hist(sqrt(alldata$HardCoral))
alldata$lCoral<-sqrt(alldata$HardCoral)
model_Coral<-lmer(lCoral~ 
                    DepthCategory+
                    ReefHabitat+
                    sSampArea+
                    (1|Locality),data=alldata)
residualsC<-resid(model_Coral)
g=ggplot(NULL)+geom_histogram(aes(x=resid(model_Coral)))+theme_classic()+xlab ("Residuals coral cover model")
h=ggplot(NULL)+geom_point(aes(x=predict(model_Coral),y=resid(model_Coral)))+theme_classic()+ylab ("Residuals")+xlab("Fitted coral cover (sqrt)")
alldata$lC_marg=alldata$lCoral-(fixef(model_Coral)["sSampArea"]*alldata$sSampArea+
                                  ifelse(alldata$ReefHabitat=="Crest",fixef(model_Coral)[c("ReefHabitatCrest")],
                                         ifelse(alldata$ReefHabitat=="Flat",fixef(model_Coral)[c("ReefHabitatFlat")],
                                                ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                       fixef(model_Coral)[c("ReefHabitatLagoon_Back reef")],0)))+
                                  ifelse(alldata$DepthCategory=="0-4m",fixef(model_Coral)[c("DepthCategory0-4m")],
                                         ifelse(alldata$DepthCategory==">10m",fixef(model_Coral)[c("DepthCategory>10m")],0)))
alldata$Coral_marg=(alldata$lC_marg)^2

#TOTAL SPECIES RICHNESS
#Note that to calculate total species richness from observed species richness, we fitted a possion-lognormal distribution(poilog package) to the 
#reef scale species abundance distributions, and estimated, by maximum likelihood, the mean and sd of the lognormal distribution. That gave us the approximate fraction of species revealed by the sample.
#Thus total species richness is the observed species richness divided by the fraction gives the estimated total species richnessSee Bulmer 1974: On Fitting the Poisson Lognormal Distribution to Species-Abundance Data for more details.
hist(log(alldata$total_sp))
alldata$ltotalsp<-log(alldata$total_sp)

model_totalsp<-lmer(ltotalsp~ 
                      DepthCategory+
                      ReefHabitat+
                      sSampArea+
                      SampMethod+
                      (1|Locality),data=alldata)
residualsts<-resid(model_totalsp)
i=ggplot(NULL)+geom_histogram(aes(x=residualsts))+theme_classic()+ xlab ("Residuals total species richness model")
j=ggplot(NULL)+geom_point(aes(x=predict(model_totalsp),y=resid(model_totalsp)))+theme_classic()+ ylab ("Residuals")+ xlab ("Fitted total species richness model (log)")
alldata$lts_marg=alldata$ltotalsp-(fixef(model_totalsp)["sSampArea"]*alldata$sSampArea+
                                     ifelse(alldata$ReefHabitat=="Crest",fixef(model_totalsp)[c("ReefHabitatCrest")],
                                            ifelse(alldata$ReefHabitat=="Flat",fixef(model_totalsp)[c("ReefHabitatFlat")],
                                                   ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                          fixef(model_totalsp)[c("ReefHabitatLagoon_Back reef")],0)))+
                                     ifelse(alldata$DepthCategory=="0-4m",fixef(model_totalsp)[c("DepthCategory0-4m")],
                                            ifelse(alldata$DepthCategory==">10m",fixef(model_totalsp)[c("DepthCategory>10m")],0))+
                                     ifelse(alldata$SampMethod=="Distance sampling",fixef(model_totalsp)[c("SampMethodDistance sampling")],
                                            ifelse(alldata$SampMethod=="Point intercept",fixef(model_totalsp)[c("SampMethodPoint intercept")],0)))

alldata$totalsp_marg=exp(alldata$lts_marg)


#parrotfish scraping potential when parrotfish are present......................................
hist(log(alldata$scraping_potential))
alldata$lherb<-ifelse(alldata$scraping_potential==0,NA,log(alldata$scraping_potential))
hist(alldata$lherb)
model_herb<-lmer(lherb~ 
                   DepthCategory+
                   ReefHabitat+
                   sSampArea+
                   SampMethod+
                   (1|Locality),data=alldata)
residualsC<-resid(model_herb)
k=ggplot(NULL)+geom_histogram(aes(x=residualsC))+theme_classic()+ xlab ("Residuals parrotfish scraping model")
l=ggplot(NULL)+geom_point(aes(x=predict(model_herb),y=resid(model_herb)))+theme_classic()+ ylab ("Residuals")+ xlab ("Fitted parrotfish scraping (log)")
alldata$lh_marg=alldata$lherb-(fixef(model_herb)["sSampArea"]*alldata$sSampArea+
                                 ifelse(alldata$ReefHabitat=="Crest",fixef(model_herb)[c("ReefHabitatCrest")],
                                        ifelse(alldata$ReefHabitat=="Flat",fixef(model_herb)[c("ReefHabitatFlat")],
                                               ifelse(alldata$ReefHabitat=="Lagoon_Back reef",
                                                      fixef(model_herb)[c("ReefHabitatLagoon_Back reef")],0)))+
                                 ifelse(alldata$DepthCategory=="0-4m",fixef(model_herb)[c("DepthCategory0-4m")],
                                        ifelse(alldata$DepthCategory==">10m",fixef(model_herb)[c("DepthCategory>10m")],0))+
                                 ifelse(alldata$SampMethod=="Distance sampling",fixef(model_herb)[c("SampMethodDistance sampling")],
                                        ifelse(alldata$SampMethod=="Point intercept",fixef(model_herb)[c("SampMethodPoint intercept")],0)))
alldata$herb_marg=exp(alldata$lh_marg)


#model fit all models that were used to correct for methodological effects
windows()
ggarrange(a,b,c,d,e,f,g,h,i,j,k,l, nrow=6,ncol=2, labels=c("a","","b","","c","","d","","e","","f",""))

#...............................................................................................

#STATUS
#Jurisdiction Biomass status ---------------------------------------

#collapsed referencepoint:0.1 of B0
collapsedB=0.1*B0_median
collapsedB_op=0.1*B0_low
collapsedB_prec=0.1*B0_high

#reef-scale biomass status
#collapsed reefs
alldata$collapsed=ifelse(alldata$Biomass_marg==collapsedB|alldata$Biomass_marg<collapsedB,1,0)
length(alldata$collapsed[alldata$collapsed==1])/length(alldata$collapsed)
alldata$collapsed_prec=ifelse(alldata$Biomass_marg==collapsedB_prec|alldata$Biomass_marg<collapsedB_prec,1,0)
length(alldata$collapsed_prec[alldata$collapsed_prec==1])/length(alldata$collapsed_prec)
alldata$collapsed_op=ifelse(alldata$Biomass_marg==collapsedB_op|alldata$Biomass_marg<collapsedB_op,1,0)
length(alldata$collapsed_op[alldata$collapsed_op==1])/length(alldata$collapsed_op)

#biomass status (B/BMMSY) at a reef scale
alldata$biomassstatus=alldata$Biomass_marg/BMMSY_median
alldata$biomassstatus_prec=alldata$Biomass_marg/BMMSY_high
alldata$biomassstatus_op=alldata$Biomass_marg/BMMSY_low
alldata$overfished=ifelse(alldata$biomassstatus>1,0,1)
length(alldata$overfished[alldata$overfished==1])/length(alldata$overfished)
alldata$overfished_prec=ifelse(alldata$biomassstatus_prec>1,0,1)
length(alldata$overfished_prec[alldata$overfished_prec==1])/length(alldata$overfished_prec)
alldata$overfished_op=ifelse(alldata$biomassstatus_op>1,0,1)
length(alldata$overfished_op[alldata$overfished_op==1])/length(alldata$overfished_op)

#calculate median biomass,number of reef sites, median biomass status and proprtion of reefs overfished by country 
Bbycountry=ddply(alldata,.(Larger),summarize,medianB_tkm2=median(Biomass_marg),  reefsites=length(lB_marg),propoverfished=mean(overfished), medianBstatus=median(biomassstatus),propoverfished_prec=mean(overfished_prec), medianBstatus_prec=median(biomassstatus_prec),propoverfished_op=mean(overfished_op), medianBstatus_op=median(biomassstatus_op))


#probability of being overfished
country_proboverfished=matrix(NA,nrow=length(list_of_draws_null$BMMSY),ncol=length(Bbycountry$medianB_tkm2))
for (i in 1:length(Bbycountry$medianB_tkm2)){
  country_proboverfished[,i]=Bbycountry$medianB_tkm2[i]/list_of_draws_null$BMMSY
  Bbycountry$prob_overfished[i]=length(country_proboverfished[,i][country_proboverfished[,i]<1])/length(country_proboverfished[,i])
}
country_proboverfished=as.data.frame(country_proboverfished)
colnames(country_proboverfished)=Bbycountry$Larger
country_proboverfished2=melt(country_proboverfished)

#ridge plot figure
prob_overfished_fig=ggplot(data=country_proboverfished2,aes(x=log(value), y=variable, fill=stat(x)))+geom_density_ridges_gradient()+geom_vline(xintercept = 0, lty=2)+
  scale_fill_gradient2("log(B/BMMSY)",low="darkred",mid="white",high="navyblue",midpoint = 0)+theme_classic()+xlab("log(B/BMMSY)")+ylab("")

#jurisdictions with median biomass values collapsed:
Bbycountry$collapsed=ifelse(Bbycountry$medianB_tkm2==collapsedB|Bbycountry$medianB_tkm2<collapsedB,1,0)
length(Bbycountry$collapsed[Bbycountry$collapsed==1])/length(Bbycountry$collapsed)
Bbycountry$Larger[Bbycountry$collapsed==1]
Bbycountry$collapsed_prec=ifelse(Bbycountry$medianB_tkm2==collapsedB_prec|Bbycountry$medianB_tkm2<collapsedB_prec,1,0)
length(Bbycountry$collapsed_prec[Bbycountry$collapsed_prec==1])/length(Bbycountry$collapsed_prec)
Bbycountry$collapsed_op=ifelse(Bbycountry$medianB_tkm2==collapsedB_op|Bbycountry$medianB_tkm2<collapsedB_op,1,0)
length(Bbycountry$collapsed_op[Bbycountry$collapsed_op==1])/length(Bbycountry$collapsed_op)

#jurisdictions overfished (with median biomass values below BMMSY): 
Bbycountry$overfished=ifelse(Bbycountry$medianBstatus>1,"Not overfished","Overfished")
length(Bbycountry$overfished[Bbycountry$overfished=="Overfished"])/length(Bbycountry$overfished)
Bbycountry$overfished_prec=ifelse(Bbycountry$medianBstatus_prec>1,"Not overfished","Overfished")
length(Bbycountry$overfished_prec[Bbycountry$overfished_prec=="Overfished"])/length(Bbycountry$overfished_prec)
Bbycountry$overfished_op=ifelse(Bbycountry$medianBstatus_op>1,"Not overfished","Overfished")
length(Bbycountry$overfished_op[Bbycountry$overfished_op=="Overfished"])/length(Bbycountry$overfished_op)


#merge with jurisdiction's scale data
jurisdictionscale_data=merge(jurisdictionscale_data,Bbycountry,by="Larger", all.x=T )

#Jurisdiction Fishing status ---------------------------------------

#calculate catch per unit area of reef
jurisdictionscale_data$catch_tkm2=jurisdictionscale_data$mean_spatial_totalcatch_reefs_t/jurisdictionscale_data$CoralReefArea_km2 

#probability of overfishing
country_proboverfishing=matrix(NA,nrow=length(list_of_draws_null$MMSY),ncol=length(jurisdictionscale_data$catch_tkm2))
for (i in 1:length(jurisdictionscale_data$catch_tkm2)){
  country_proboverfishing[,i]=jurisdictionscale_data$catch_tkm2[i]/list_of_draws_null$MMSY
  jurisdictionscale_data$prob_overfishing[i]=length(country_proboverfishing[,i][country_proboverfishing[,i]>1])/length(country_proboverfishing[,i])
}

country_proboverfishing=as.data.frame(country_proboverfishing)
colnames(country_proboverfishing)=jurisdictionscale_data$area_name
jurisdictionscale_data$prob_overfishing=ifelse(is.na(jurisdictionscale_data$catch_tkm2),NA,jurisdictionscale_data$prob_overfishing)
country_proboverfishing2=melt(country_proboverfishing)
country_proboverfishing3=country_proboverfishing2[!is.na(country_proboverfishing2$value),]

#ridge plot
prob_overfishing_fig=ggplot(data=country_proboverfishing3,aes(x=log(value), y=variable, fill=stat(x)))+geom_density_ridges_gradient()+geom_vline(xintercept = 0, lty=2)+
  scale_fill_gradient2("log(C/MMSY)",high="darkred",mid="white",low="navyblue",midpoint = 0)+theme_classic()+xlab("log(C/MMSY)")+ylab("")+theme_classic()+
  theme(axis.text.y = element_text(size=6))

#jurisdictions fishing status: 
jurisdictionscale_data$fishingstatus=jurisdictionscale_data$catch_tkm2/MMSY_median
jurisdictionscale_data$overfishing=ifelse(is.na(jurisdictionscale_data$fishingstatus),NA,ifelse(jurisdictionscale_data$fishingstatus>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$overfishing)&jurisdictionscale_data$overfishing=="Overfishing"])/length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$overfishing)])
jurisdictionscale_data$fishingstatus_prec=jurisdictionscale_data$catch_tkm2/MMSY_low
jurisdictionscale_data$overfishing_prec=ifelse(is.na(jurisdictionscale_data$fishingstatus),NA,ifelse(jurisdictionscale_data$fishingstatus_prec>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data$overfishing_prec[!is.na(jurisdictionscale_data$overfishing_prec)&jurisdictionscale_data$overfishing_prec=="Overfishing"])/length(jurisdictionscale_data$overfishing_prec[!is.na(jurisdictionscale_data$overfishing_prec)])
jurisdictionscale_data$fishingstatus_op=jurisdictionscale_data$catch_tkm2/MMSY_high
jurisdictionscale_data$overfishing_op=ifelse(is.na(jurisdictionscale_data$fishingstatus),NA,ifelse(jurisdictionscale_data$fishingstatus_op>1,"Overfishing","Not overfishing"))
length(jurisdictionscale_data$overfishing_op[!is.na(jurisdictionscale_data$overfishing_op)&jurisdictionscale_data$overfishing_op=="Overfishing"])/length(jurisdictionscale_data$overfishing_op[!is.na(jurisdictionscale_data$overfishing_op)])

#Jurisdiction Sustainability status --------------------------------

#make the biomass overlap the biomass of the surplus production curve by rounding it to the nearest integer
jurisdictionscale_data$roundedmedianB=round(jurisdictionscale_data$medianB_tkm2)
jurisdictionscale_data2<-jurisdictionscale_data
jurisdictionscale_data2$Biomass=jurisdictionscale_data2$roundedmedianB

#merge
jurisdictionscale_data2=merge(jurisdictionscale_data2,popgrowth_null,by="Biomass",all.x=T)

#give almost 0 yield to those biomass values that had no surplus production based on our benchmark model
jurisdictionscale_data2$avyield=ifelse(is.na(jurisdictionscale_data2$median),NA,jurisdictionscale_data2$median)
jurisdictionscale_data2$highyield=ifelse(is.na(jurisdictionscale_data2$`75%`),NA, ifelse(jurisdictionscale_data2$`75%`==0,0.0001,jurisdictionscale_data2$`75%`))
jurisdictionscale_data2$lowyield=ifelse(is.na(jurisdictionscale_data2$`25%`),NA, ifelse(jurisdictionscale_data2$`25%`==0,0.0001,jurisdictionscale_data2$`25%`))

#calculate fishing status based on the surplus production curve
jurisdictionscale_data2$fishingstatus_curve=jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$avyield
length(jurisdictionscale_data2$fishingstatus_curve[!is.na(jurisdictionscale_data2$fishingstatus_curve)])

jurisdictionscale_data2$overfishing_curve=ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve>1,1,0))
jurisdictionscale_data2$fishingstatus_curve_prec=jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$lowyield
jurisdictionscale_data2$overfishing_curve_prec=ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_prec),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1,1,0))
jurisdictionscale_data2$fishingstatus_curve_op=jurisdictionscale_data2$catch_tkm2/jurisdictionscale_data2$highyield
jurisdictionscale_data2$overfishing_curve_op=ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_op),NA,ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1,1,0))

#biomass status the same
jurisdictionscale_data2$biomassstatus_curve=jurisdictionscale_data2$medianB_tkm2/BMMSY_median 
jurisdictionscale_data2$biomassstatus_curve_prec=jurisdictionscale_data2$medianB_tkm2/BMMSY_high 
jurisdictionscale_data2$biomassstatus_curve_op=jurisdictionscale_data2$medianB_tkm2/BMMSY_low 

#sustainability status 
jurisdictionscale_data2$fisherystatus=as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve<1 & jurisdictionscale_data2$biomassstatus_curve>1, "Sustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve>1 & jurisdictionscale_data2$biomassstatus_curve<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve>1 & jurisdictionscale_data2$biomassstatus_curve>1, "Warning", "Rebuilding")))))
jurisdictionscale_data2$fisherystatus_prec=as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_prec),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve_prec<1 & jurisdictionscale_data2$biomassstatus_curve_prec>1, "Sustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1 & jurisdictionscale_data2$biomassstatus_curve_prec<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_prec>1 & jurisdictionscale_data2$biomassstatus_curve_prec>1, "Warning", "Rebuilding")))))
jurisdictionscale_data2$fisherystatus_op=as.factor(ifelse(is.na(jurisdictionscale_data2$fishingstatus_curve_op),NA, ifelse(jurisdictionscale_data2$fishingstatus_curve_op<1 & jurisdictionscale_data2$biomassstatus_curve_op>1, "Sustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1 & jurisdictionscale_data2$biomassstatus_curve_op<1, "Unsustainable", ifelse(jurisdictionscale_data2$fishingstatus_curve_op>1 & jurisdictionscale_data2$biomassstatus_curve_op>1, "Warning", "Rebuilding")))))
#(incorporating quasi-sustainable) as sustainable
length(jurisdictionscale_data3$area_name[jurisdictionscale_data2$fisherystatus=="Warning" &jurisdictionscale_data3$fishingstatus<1])/length(jurisdictionscale_data3$area_name[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data2$fisherystatus=ifelse(is.na(jurisdictionscale_data2$fisherystatus),NA,ifelse(jurisdictionscale_data2$fisherystatus=="Warning" &jurisdictionscale_data2$fishingstatus<1, "Sustainable",as.character(jurisdictionscale_data2$fisherystatus)))
jurisdictionscale_data2$fisherystatus_prec=ifelse(is.na(jurisdictionscale_data2$fisherystatus_prec),NA,ifelse(jurisdictionscale_data2$fisherystatus_prec=="Warning" &jurisdictionscale_data2$fishingstatus<1, "Sustainable",as.character(jurisdictionscale_data2$fisherystatus_prec)))
jurisdictionscale_data2$fisherystatus_op=ifelse(is.na(jurisdictionscale_data2$fisherystatus_op),NA,ifelse(jurisdictionscale_data2$fisherystatus_op=="Warning" &jurisdictionscale_data2$fishingstatus<1, "Sustainable",as.character(jurisdictionscale_data2$fisherystatus_op)))

#separating only those that have both biomass and catch
jurisdictionscale_data3=jurisdictionscale_data2[!is.na(jurisdictionscale_data2$fisherystatus),]
jurisdictionscale_data3=droplevels(jurisdictionscale_data3)

#percentages in each category 
#sustainable
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="Sustainable"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="Sustainable"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="Sustainable"]

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
colnames(jurisdictionscale_data3)

#potentially rebuilding
length(jurisdictionscale_data3$fisherystatus[jurisdictionscale_data3$fisherystatus=="Rebuilding"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus=="Rebuilding"]
length(jurisdictionscale_data3$fisherystatus_prec[jurisdictionscale_data3$fisherystatus_prec=="Rebuilding"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_prec=="Rebuilding"]
length(jurisdictionscale_data3$fisherystatus_op[jurisdictionscale_data3$fisherystatus_op=="Rebuilding"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])
jurisdictionscale_data3$area_name[jurisdictionscale_data3$fisherystatus_op=="Rebuilding"]

#have past one or both MMSY benchmarks
length(jurisdictionscale_data3$fisherystatus[!jurisdictionscale_data3$fisherystatus=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus[!is.na(jurisdictionscale_data3$fisherystatus)])
length(jurisdictionscale_data3$fisherystatus_prec[!jurisdictionscale_data3$fisherystatus_prec=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus_prec[!is.na(jurisdictionscale_data3$fisherystatus_prec)])
length(jurisdictionscale_data3$fisherystatus_op[!jurisdictionscale_data3$fisherystatus_op=="Sustainable"])/length(jurisdictionscale_data3$fisherystatus_op[!is.na(jurisdictionscale_data3$fisherystatus_op)])

#create a data frame to store all the values
status_results=as.data.frame(c(length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$catch_tkm2)]),length(Bbycountry$collapsed),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                               length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$catch_tkm2)&jurisdictionscale_data$overfishing=="Overfishing"])/length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$catch_tkm2)]),
                               length(Bbycountry$collapsed[Bbycountry$collapsed==1])/length(Bbycountry$collapsed),
                               length(Bbycountry$overfished[Bbycountry$overfished=="Overfished"])/length(Bbycountry$overfished),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Sustainable"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Unsustainable"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Warning"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&jurisdictionscale_data2$fisherystatus=="Rebuilding"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)]),
                               length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)&!jurisdictionscale_data2$fisherystatus=="Sustainable"])/length(jurisdictionscale_data2$fisherystatus[!is.na(jurisdictionscale_data2$fisherystatus)])))
row.names(status_results)=c("n_fishingstatus","n_biomassstatus","n_fisherystatus","overfishing","collapsed","overfished","sustainable","unsustainable","warning","rebuilding","conservation_concern")
#write.csv(status_results, "main_status_results_quasi.csv")

#map of status
newmap <- getMap(resolution = "high")

Fig.1a<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data[!is.na(jurisdictionscale_data$overfishing),],aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$overfishing)]),colour="black", pch=21,size=3, alpha=0.9)+
  scale_fill_manual(name = "Fishing status
(C/MMSY)",
                    values = c("Not overfishing" = "navyblue",
                               "Overfishing" = "red"),
                    labels = c("Not overfishing" ="Not overfishing", 
                               "Overfishing" ="Overfishing"))+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic() 

Fig.1b<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group), fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data, aes(x=EEZ_long, y=EEZ_lat, fill =jurisdictionscale_data$propoverfished, size=jurisdictionscale_data$reefsites),col="black", pch=21)+
  scale_fill_gradient2(low="navyblue",mid="white",high= "red",midpoint=0.5,
                       name="Proportion of
reefs where
B/B_MMSY<1
                       ")+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic()

Fig.1c<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),  fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data3, aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data3$fisherystatus, size=jurisdictionscale_data3$reefsites),col="black", pch=21,  alpha=0.9)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  scale_fill_manual (name="Fishery status",values=c( "Rebuilding"="cyan3",
                                                     "Sustainable"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="goldenrod1"),
                     labels = c("Rebuilding"="Potentially
rebuilding",
                                "Sustainable"="Sustainable",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"))+
  guides(fill = guide_legend(override.aes = list(size=3)))+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  xlab("")+
  ylab("")+theme_classic()

windows()
fig1b=ggarrange(Fig.1a,Fig.1b,Fig.1c, ncol=1, nrow=3,heights=c(1,1,1), labels=c("d","e","f"))
ggarrange(fig1a,fig1b,nrow=1,ncol=2,widths=c(1,3))


#optimistic map
Fig.1aO<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data[!is.na(jurisdictionscale_data$overfishing),],aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data$overfishing_op[!is.na(jurisdictionscale_data$overfishing)]),colour="black", pch=21,size=3, alpha=0.9)+
  scale_fill_manual(name = "Fishing status
(C/MMSY)",
                    values = c("Not overfishing" = "navyblue",
                               "Overfishing" = "red"),
                    labels = c("Not overfishing" ="Not overfishing", 
                               "Overfishing" ="Overfishing"), guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic() 

Fig.1bO<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group), fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data, aes(x=EEZ_long, y=EEZ_lat, fill =jurisdictionscale_data$propoverfished_op, size=jurisdictionscale_data$reefsites),col="black", pch=21)+
  scale_fill_gradient2(low="navyblue",mid="white",high= "red",midpoint=0.5,
                       name="Proportion of
reefs where
B < B_MMSY", guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic()

Fig.1cO<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),  fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data3, aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data3$fisherystatus_op, size=jurisdictionscale_data3$reefsites),col="black", pch=21,  alpha=0.9)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  scale_fill_manual (name="Fishery status",values=c( "Rebuilding"="turquoise2",
                                                     "Sustainable"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="yellow"),
                     labels = c("Rebuilding"="Potentially
rebuilding",
                                "Sustainable"="Sustainable",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  xlab("")+
  ylab("")+theme_classic()

#PRECAUTIONARY
Fig.1aP<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data[!is.na(jurisdictionscale_data$overfishing),],aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data$overfishing_prec[!is.na(jurisdictionscale_data$overfishing)]),colour="black", pch=21,size=3, alpha=0.9)+
  scale_fill_manual(name = "Fishing status
(C/MMSY)",
                    values = c("Not overfishing" = "navyblue",
                               "Overfishing" = "red"),
                    labels = c("Not overfishing" ="Not overfishing", 
                               "Overfishing" ="Overfishing"), guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic() 

Fig.1bP<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group), fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data, aes(x=EEZ_long, y=EEZ_lat, fill =jurisdictionscale_data$propoverfished_prec, size=jurisdictionscale_data$reefsites),col="black", pch=21)+
  scale_fill_gradient2(low="navyblue",mid="white",high= "red",midpoint=0.5,
                       name="Proportion of
reefs where
B < B_MMSY
                       ", guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme_classic()

Fig.1cP<-ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),  fill = "grey", color = "grey", alpha=0.7)+
  coord_fixed(xlim = c(180, -180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=jurisdictionscale_data3, aes(x=EEZ_long, y=EEZ_lat, fill = jurisdictionscale_data3$fisherystatus_prec, size=jurisdictionscale_data3$reefsites),col="black", pch=21,  alpha=0.9)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+
  scale_fill_manual (name="Fishery status",values=c( "Rebuilding"="turquoise2",
                                                     "Sustainable"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="yellow"),
                     labels = c("Rebuilding"="Potentially
rebuilding",
                                "Sustainable"="Sustainable",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  xlab("")+
  ylab("")+theme_classic()


windows()
ggarrange(Fig.1aP,Fig.1aO,Fig.1bP,Fig.1bO,Fig.1cP,Fig.1cO, ncol=2, nrow=3,heights=c(1,1,1), labels=c("a","b"))


#correlation among fishing and biomass status
#KOBE PLOT
cor(log(jurisdictionscale_data3$medianBstatus),log(jurisdictionscale_data3$fishingstatus),method="pearson")
kobe=ggplot(jurisdictionscale_data3, aes(x=log(medianBstatus), y=log(fishingstatus)))+geom_point(aes(fill=jurisdictionscale_data3$fisherystatus,size=jurisdictionscale_data3$reefsites), pch=21)+
  scale_fill_manual (name="Fishery status",values=c( "Rebuilding"="turquoise2",
                                                     "Sustainable"="navyblue",
                                                     "Unsustainable"="red", 
                                                     "Warning"="yellow"),
                     labels = c("Rebuilding"="Potentially
                                rebuilding",
                                "Sustainable"="Sustainable",
                                "Unsustainable"="Unsustainable", 
                                "Warning"="Warning"), guide=FALSE)+
  scale_size_continuous(name="Number of biomass sites per region", guide=FALSE)+
  theme_classic()+xlab("log(B/BMMSY)")+ylab("log(C/MMSY)")+geom_hline(yintercept=0,lty=2,lwd=1.5,col="navyblue")+geom_vline(xintercept=0,lty=2,lwd=1.5,col="navyblue")+geom_text_repel(aes(label=jurisdictionscale_data3$Larger))

windows()
ggarrange(prob_overfishing_fig,prob_overfished_fig,kobe,nrow=1,ncol=3,labels=c("a","b","c"))

#Ecosystem status and trade-offs
alldata$Management=as.factor(ifelse(alldata$Management=="Remote","Restricted", as.character(alldata$Management)))
alldata=droplevels(alldata)
reefscale_data$nonzeroherb=ifelse(reefscale_data$scraping_potential==0,NA,reefscale_data$scraping_potential)

#correlation of metrics (and also correlation with mean_length)
windows()
pairs(~log(herb_marg)+sqrt(Coral_marg)+PA_tpmarg+
        log(totalsp_marg)+Size_marg ,  data= alldata,lower.panel=panel.cor, pch = 21, bg = "darkgrey",labels=c("log(Parrotfish
scrapping potential)","sqrt(Hard coral cover)", "Presence/absence 
top predators","log(Total species 
richness)", "Mean length"),cex.labels=1.5,font.labels=2,diag.panel =panel.hist, hist.col="grey" )


#"Pretty good multispecies yield"
PGY=0.8*MMSY_median
B_pgy=max(popgrowth_null$Biomass[popgrowth_null$median>PGY], na.rm=T)

#individual metric relationships 
si_1=ggplot(alldata)+
  geom_point(aes(lB_marg,Size_marg,shape=alldata$Management,alpha=alldata$Management),col="darkgreen") +
  scale_alpha_manual(values=c("Fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,Size_marg),formula=y~s(x,k=3),method="gam", col="darkgreen",fill="green",alpha=0.3)+
  theme_classic()+xlab(")")+ylab("Mean length (cm)")+scale_y_continuous(position = "left")+coord_cartesian( expand=c(0,0))

si_2=ggplot(alldata)+
  geom_point(aes(lB_marg,PA_tpmarg,shape=alldata$Management,alpha=alldata$Management),col="navyblue") +
  scale_shape_manual(values=c("Fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,PA_tpmarg),formula=y~s(x,k=3),method="gam", col="navyblue",fill="blue",alpha=0.3)+
  theme_classic()+xlab("")+ylab("Prob. top predators")+scale_y_continuous(position = "left")+coord_cartesian( expand=c(0,0))
si_3=ggplot(alldata)+
  geom_point(aes(lB_marg,log(totalsp_marg),shape=alldata$Management,alpha=alldata$Management),col="goldenrod") +
  scale_shape_manual(values=c("Fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,log(totalsp_marg)),formula=y~s(x,k=3),method="gam", col="goldenrod",fill="yellow",alpha=0.3)+
  theme_classic()+xlab("")+ylab("log(Sp richness)")+scale_y_continuous(position = "left")+coord_cartesian( expand=c(0,0))
si_4=ggplot(alldata)+
  geom_point(aes(lB_marg,sqrt(Coral_marg),shape=alldata$Management,alpha=alldata$Management),col="darkmagenta") +
  scale_shape_manual(values=c("Fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,sqrt(Coral_marg)),formula=y~s(x,k=3),method="gam", col="darkmagenta",fill="magenta",alpha=0.3)+
  theme_classic()+xlab("")+ylab("sqrt(Coral cover)")+scale_y_continuous(position = "left",limits=c(0,10))+coord_cartesian( expand=c(0,0))
si_5=ggplot(alldata)+
  geom_point(aes(lB_marg,lh_marg,shape=alldata$Management,alpha=alldata$Management),col="violetred4") +
  scale_shape_manual(values=c("Fished"=16,"Restricted"=17))+
  scale_alpha_manual(values=c("Fished"=0.2,
                              "Restricted"=0.5))+guides(alpha=FALSE, shape=FALSE)+
  geom_smooth(aes(lB_marg,lh_marg),formula=y~s(x,k=3),method="gam", col="violetred4",fill="violetred",alpha=0.3)+
  theme_classic()+xlab("log(Biomass (t/km2))")+ylab("log(Scrapping (cm^2/min))")+scale_y_continuous(position = "left")+coord_cartesian( expand=c(0,0))

#Ecosystem metric distributions by management------------------------------
# we overlaid the observes distributions for remote reefs
a=ggplot(NULL) + 
  geom_density(aes(x=reefscale_data$meanSize_cm[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0]),lty=2)+
  geom_density(aes(x=alldata$Size_marg, fill=alldata$Management,alpha=alldata$Management),col="darkgreen",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="green","Restricted"="darkgreen"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+ scale_x_continuous(breaks = seq(10, 35, by = 5),limits=c(10,35))+
  labs(x="", y="") +
  coord_flip()+theme_classic()

b=ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$total_sp[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0])),lty=2)+
  geom_density(aes(x=log(alldata$totalsp_marg), fill=alldata$Management,alpha=alldata$Management),col="goldenrod",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="yellow","Restricted"="goldenrod"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(1, 7, by = 1))+xlim(c(1,7))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
c=ggplot(NULL) + 
  geom_density(aes(x=sqrt(reefscale_data$HardCoral[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0])),lty=2)+
  geom_density(aes(x=sqrt(alldata$Coral_marg), fill=alldata$Management,alpha=alldata$Management),col="magenta4",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="magenta","Restricted"="magenta4"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
d=ggplot(NULL) + 
  geom_density(aes(x=log(reefscale_data$scraping_potential[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0])),lty=2)+
  geom_density(aes(x=log(alldata$herb_marg), fill=alldata$Management,alpha=alldata$Management),col="violetred4",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="violetred2","Restricted"="violetred4"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 10, by = 2),limits=c(0,10))+
  labs(x="", y="") +
  coord_flip()+theme_classic()

e=ggplot(NULL) + geom_density(aes(x=reefscale_data$PA_tp[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0]),lty=2)+
  
  geom_density(aes(x=alldata$PA_tpmarg,fill=alldata$Management,alpha=alldata$Management),col="navyblue",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="blue","Restricted"="navyblue"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 1, by = 1))+
  labs(x="", y="") +
  coord_flip()+theme_classic()
f=ggplot(NULL) + geom_density(aes(x=log(reefscale_data$FamBiomass_tkm2[reefscale_data$remote_20h==1 & reefscale_data$inhabited==0] )),lty=2)+
  
  geom_density(aes(log(alldata$Biomass_marg),fill=alldata$Management,alpha=alldata$Management),col="darkred",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="red","Restricted"="darkred"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+#xlim(0,300)+
  labs(x="log(Biomass (t/km2))", y="") +theme_classic()+coord_flip()

windows()
ggarrange(si_1,a, si_2,e,si_3,b,si_4,c,si_5,d, nrow=5,ncol=2,widths=c(1.5,1))

#now we get the distributions covering the range in our suplus production curves
alldata2=alldata[alldata$Biomass_marg<301,]

a=ggplot(NULL) + 
  geom_density(aes(x=alldata2$Size_marg, fill=alldata2$Management,alpha=alldata2$Management),col="darkgreen",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="green","Restricted"="darkgreen"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+xlim(c(10,30))+
  guides(alpha=F,fill=F)+ # scale_x_continuous(breaks = seq(10, 40, by = 5))+
  labs(x="mean Length (cm)", y="") +
  coord_flip()+theme_classic()

b=ggplot(NULL) + 
  geom_density(aes(x=alldata2$totalsp_marg, fill=alldata2$Management,alpha=alldata2$Management),col="goldenrod",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="yellow","Restricted"="goldenrod"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  xlim(c(0,200))+
  labs(x="Total species richness", y="") +
  coord_flip()+theme_classic()
c=ggplot(NULL) + 
  geom_density(aes(x=alldata2$Coral_marg, fill=alldata2$Management,alpha=alldata2$Management),col="magenta4",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="magenta","Restricted"="magenta4"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+xlim(c(0,100))+
  guides(alpha=F,fill=F)+  #scale_x_continuous(breaks = seq(0, 10, by = 2))+
  labs(x="Hard Coral cover", y="") +
  coord_flip()+theme_classic()
d=ggplot(NULL) + 
  geom_density(aes(x=alldata2$herb_marg, fill=alldata2$Management,alpha=alldata2$Management),col="violetred4",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="violetred2","Restricted"="violetred4"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+xlim(c(0,750))+
  guides(alpha=F,fill=F)+  #scale_x_continuous(breaks = seq(0, 10, by = 2))+
  labs(x="Parrotfish scrapping potential (cm^2/min)", y="") +
  coord_flip()+theme_classic()

e=ggplot(NULL) + 
  geom_density(aes(x=alldata2$PA_tpmarg,fill=alldata2$Management,alpha=alldata2$Management),col="navyblue",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="blue","Restricted"="navyblue"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+  scale_x_continuous(breaks = seq(0, 1, by = 1))+
  labs(x="Prob. encountering top predators", y="") +
  coord_flip()+theme_classic()
f=ggplot(NULL) + 
  geom_density(aes(alldata2$Biomass_marg,fill=alldata2$Management,alpha=alldata2$Management),col="darkred",  
               na.rm=T) + scale_fill_manual(values=c("Fished"="red","Restricted"="darkred"))+
  scale_alpha_manual(values=c("Fished"=0.3,"Restricted"=0.3))+
  guides(alpha=F,fill=F)+#xlim(0,300)+
  labs(x="Biomass (t/km2)", y="") +theme_classic()+coord_flip()

windows()
ggarrange(a,b,c,d,e,f,nrow=2, ncol=3)

#Generalized additive models------------------------------------------------
#Ecosystem metrics in the manuscript

#FISH LENGTH (cm)
gam_size=gam(Size_marg~s(lB_marg,k=3),data=alldata)

#model fit
hist(resid(gam_size))  

#predicted length
size_pred=as.data.frame(predict(gam_size,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
size_pred$lB_marg=sort(alldata$lB_marg)
size_pred$upr=size_pred$fit+(2*size_pred$se.fit)
size_pred$low=size_pred$fit-(2*size_pred$se.fit)

#biomass arithmetic scale (for the trade-offs graph)
size_pred$Biomass=exp(size_pred$lB_marg)

#Length at the different reference points
currentB=median(alldata$Biomass_marg)
ScurrentB=size_pred$fit[size_pred$Biomass>currentB][1]
SBMSY=size_pred$fit[size_pred$Biomass>BMMSY_median][1]
SPGY=size_pred$fit[size_pred$Biomass>B_pgy][1]
SBO=size_pred$fit[size_pred$Biomass>B0_median][1]

#change with respect to unfished biomass conditions (unfished biomass 100%)
100-(ScurrentB*100)/SBO
100-(SPGY*100)/SBO
100-(SBMSY*100)/SBO
(100-(SBMSY*100)/SBO)-(100-(ScurrentB*100)/SBO)
(100-(SPGY*100)/SBO)-(100-(SBMSY*100)/SBO)
(100-(SBMSY*100)/SBO)-(100-(SPGY*100)/SBO)

#TOP PREDATOR PRESENCE/ABSENCE
gam_tp=gam(PA_tpmarg~s(lB_marg,k=3),data=alldata)
plot(fitted(gam_tp), resid(gam_tp))  
tp_pred=as.data.frame(predict(gam_tp,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
tp_pred$lB_marg=sort(alldata$lB_marg)
tp_pred$upr=tp_pred$fit+(2*tp_pred$se.fit)
tp_pred$low=tp_pred$fit-(2*tp_pred$se.fit)
tp_pred$Biomass=exp(tp_pred$lB_marg)
pBMSY=tp_pred$fit[tp_pred$Biomass>BMMSY_median][1]
pPGY=tp_pred$fit[tp_pred$Biomass>B_pgy][1]
pBO=tp_pred$fit[tp_pred$Biomass>B0_median][1]
pcurrentB=tp_pred$fit[tp_pred$Biomass>currentB][1]
100-(pcurrentB*100)/pBO
100-(pPGY*100)/pBO
100-(pBMSY*100)/pBO
(100-(pPGY*100)/pBO)-(100-(pBMSY*100)/pBO)
(100-(pcurrentB*100)/pBO)-(100-(pBMSY*100)/pBO)
(100-(pBMSY*100)/pBO)-(100-(pPGY*100)/pBO)

#TOTAL SPECIES RICHNESS
gam_tr=gam(log(totalsp_marg)~s(lB_marg,k=3),data=alldata)
hist(resid(gam_tr))  
tr_pred=as.data.frame(predict(gam_tr,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
tr_pred$lB_marg=sort(alldata$lB_marg)
tr_pred$fit=exp(tr_pred$fit)
tr_pred$upr=tr_pred$fit+(2*exp(tr_pred$se.fit))
tr_pred$low=tr_pred$fit-(2*exp(tr_pred$se.fit))
tr_pred$Biomass=exp(tr_pred$lB_marg)
rBMSY=tr_pred$fit[tr_pred$Biomass>BMMSY_median][1]
rPGY=tr_pred$fit[tr_pred$Biomass>B_pgy][1]
rBO=tr_pred$fit[tr_pred$Biomass>B0_median][1]
rcurrentB=tr_pred$fit[tr_pred$Biomass>currentB][1]
100-(rcurrentB*100)/rBO
100-(rPGY*100)/rBO
100-(rBMSY*100)/rBO
(100-(rPGY*100)/rBO)-(100-(rBMSY*100)/rBO)
(100-(rcurrentB*100)/rBO)-(100-(rBMSY*100)/rBO)
(100-(rBMSY*100)/rBO)-(100-(rPGY*100)/rBO)

#herb
gam_her=gam(log(herb_marg)~s(lB_marg,k=3),data=alldata)
hist(resid(gam_her))  
her_pred=as.data.frame(predict(gam_her,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
her_pred$lB_marg=sort(alldata$lB_marg)
her_pred$fit=exp(her_pred$fit)
her_pred$upr=her_pred$fit+(2*exp(her_pred$se.fit))
her_pred$low=her_pred$fit-(2*exp(her_pred$se.fit))
her_pred$Biomass=exp(her_pred$lB_marg)
hBMSY=her_pred$fit[her_pred$Biomass>BMMSY_median][1]
hPGY=her_pred$fit[her_pred$Biomass>B_pgy][1]
hBO=her_pred$fit[her_pred$Biomass>B0_median][1]
hcurrentB=her_pred$fit[her_pred$Biomass>currentB][1]
100-(hcurrentB*100)/hBO
100-(hPGY*100)/hBO
100-(hBMSY*100)/hBO
(100-(hBMSY*100)/hBO)-(100-(hPGY*100)/hBO)
(100-(hcurrentB*100)/hBO)-(100-(hBMSY*100)/hBO)
(100-(hBMSY*100)/hBO)-(100-(hPGY*100)/hBO)

#coral cover
gam_coral=gam(sqrt(Coral_marg)~s(lB_marg,k=3),data=alldata)
hist(resid(gam_coral))  
coral_pred=as.data.frame(predict(gam_coral,newdata=list(lB_marg= sort(alldata$lB_marg)),se.fit=2))
coral_pred$lB_marg=sort(alldata$lB_marg)
coral_pred$fit=(coral_pred$fit)^2
coral_pred$upr=coral_pred$fit+(2*(coral_pred$se.fit)^2)
coral_pred$low=coral_pred$fit-(2*(coral_pred$se.fit)^2)
coral_pred$Biomass=exp(coral_pred$lB_marg)
ccurrentB=coral_pred$fit[coral_pred$Biomass>currentB][1]
cBMSY=coral_pred$fit[coral_pred$Biomass>BMMSY_median][1]
cPGY=coral_pred$fit[coral_pred$Biomass>B_pgy][1]
cBO=coral_pred$fit[coral_pred$Biomass>B0_median][1]
100-(ccurrentB*100)/cBO
100-(cPGY*100)/cBO
100-(cBMSY*100)/cBO

#trade-offs figure matching with the distributions
windows()
a=ggplot(NULL)+geom_line(aes(y=popgrowth_null2$median[!is.na(popgrowth_null2$median)], x=popgrowth_null2$Biomass[!is.na(popgrowth_null2$median)]), col="black", lwd=2)+
  geom_ribbon(aes(x=popgrowth_null$Biomass,ymin=popgrowth_null$`5%`[!is.na(popgrowth_null$`5%`)],ymax=popgrowth_null$`95%`[!is.na(popgrowth_null$`95%`)]), fill="grey", alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ylim(c(0,3.5))+xlim(c(0,300))
b=ggplot(size_pred[size_pred$Biomass<300,])+
  geom_line(aes(y=fit, x=Biomass), col='darkgreen',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='green', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("mean Size")+  coord_cartesian(xlim=c(0,300),ylim=c(10,30),expand=c(0,0))+scale_y_continuous(position = "right")
c=ggplot(tp_pred[tp_pred$Biomass<300,])+
  geom_line(aes(y=fit, x=Biomass), col='navyblue',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='blue', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Prb. top predators")+  coord_cartesian(xlim=c(0,300),ylim=c(0,1),expand=c(0,0))+scale_y_continuous(position = "right")
d=ggplot(tr_pred[tr_pred$Biomass<300,])+
  geom_line(aes(y=fit, x=Biomass), col='goldenrod',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='yellow', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Total species richness")+  coord_cartesian(xlim=c(0,300),ylim=c(0,200),expand=c(0,0))+scale_y_continuous(position = "right")
e=ggplot(her_pred[her_pred$Biomass<300,])+
  geom_line(aes(y=fit, x=Biomass), col='violetred4',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='violetred3', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Parrotfish scraping potential (cm^2/min)")+  coord_cartesian(xlim=c(0,300),ylim=c(0,750),expand=c(0,0))+scale_y_continuous(position = "right")
f=ggplot(coral_pred[coral_pred$Biomass<300,])+
  geom_line(aes(y=fit, x=Biomass), col='darkmagenta',lwd=1.3)+
  geom_ribbon(aes(ymin=low, ymax=upr, x=Biomass), fill='magenta', alpha=0.3)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Hard Coral Cover")+  coord_cartesian(xlim=c(0,300),ylim=c(0,100),expand=c(0,0))+scale_y_continuous(position = "right")

ggarrange(a,b,c,d,e,f,nrow=2,ncol=3)


#######
#Potential jurisdiction-scale explanatory variables of observed status...............................................................................................

#check for correlations among socio-economic explanatory variables
windows()
pairs(~ log(meanTgrav_nh2+1)+log(Tourists_km2l)+log(Rfishers_km2r+1)+HDI+popgrow_prop+log(mpa_perc+1)+log(meanTTmarket_h), data=jurisdictionscale_data,lower.panel=panel.cor )

#check correlations with population size (which is used by the SAUP between anchor points): Total gravity is slighty correlated because it includes the population size within 500km2 of reefs.
#However, we keep it in because, theoretically, we know it is a good predictor of biomass (Cinner et al. 2018)
pairs(~ log(pop_n+1)+log(meanTgrav_nh2+1)+log(Tourists_km2l)+log(Rfishers_km2r+1)+HDI+popgrow_prop+log(mpa_perc+1)+log(meanTTmarket_h), data=jurisdictionscale_data,lower.panel=panel.cor,
      pch = 21, bg = "darkgrey",labels=c("Population size (log+1)","Total gravity(log+1)","Tourist density (log)","Fisher density (log+1)" ,"HDI","Population growth","% protected waters (log+1)","Travel time markets (log)"),cex.labels=1,font.labels=1)

#variance inflation factors below 2.2
vif(lm(log(fishingstatus)~ log(meanTgrav_nh2+1)+log(Tourists_km2l)+log(Rfishers_km2r+1)+HDI+popgrow_prop+log(mpa_perc+1)+log(meanTTmarket_h), data=jurisdictionscale_data))

#response variables (fishing status and biomass status)
#transform response variables
jurisdictionscale_data$lfishingstatus=log(jurisdictionscale_data$fishingstatus)
jurisdictionscale_data$lbiomassstatus=log(jurisdictionscale_data$medianBstatus)

#alternative binomial models for wheter jurisdictions were overfished/overfishing or not
#make dummy variables
jurisdictionscale_data$overfishing_dummy=ifelse(is.na(jurisdictionscale_data$overfishing),NA,ifelse(jurisdictionscale_data$overfishing=="Overfishing",1,0))
jurisdictionscale_data$overfished_dummy=ifelse(is.na(jurisdictionscale_data$overfished),NA,ifelse(jurisdictionscale_data$overfished=="Overfished",0,1))


#standardise explanatory variables and transform if neccesary
jurisdictionscale_data=jurisdictionscale_data %>%
  mutate(sfisherdens=standardise(log(Rfishers_km2r+1)),
         stouristdens=standardise(log(Tourists_km2l)),
         sHDI=standardise(HDI),
         spopgrowth=standardise(popgrow_prop),
         smeanttm=standardise(log(meanTTmarket_h)),
         smpa=standardise(log(mpa_perc+1)),
         stotalgravity=standardise(log(meanTgrav_nh2+1)))

#explanatory variables data
explanatorydata=jurisdictionscale_data[,c("area_name","stotalgravity","smeanttm","sfisherdens","stouristdens","sHDI", "spopgrowth","smpa" )]

#check missingness
apply(explanatorydata,2,pMiss)
mice::md.pattern(explanatorydata)
VIM::aggr(explanatorydata, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(explanatorydata), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#imputing missing data by predictive mean matching
impdata=mice(explanatorydata,m=5,maxit=50,meth='pmm',seed=500)
#examine imputed data
stripplot(impdata, pch = 20, cex = 1.2)

#get one of the runs to complete the explanatory data
completeimputeddata=mice::complete(impdata,3)
colnames(completeimputeddata)=c("area_name", "stotalgravity_imp","smeanttm_imp","sfisherdens_imp","stouristdens_imp","sHDI_imp", "spopgrowth_imp","smpa_imp" )

#combine imputations with jurisdiction level data
jurisdictionscale_data=merge(jurisdictionscale_data,completeimputeddata, by="area_name", all.x=T)
a=ggplot(jurisdictionscale_data,aes(x=prob_overfishing))+geom_histogram()+theme_classic()
b=ggplot(jurisdictionscale_data,aes(x=prob_overfished))+geom_histogram()+theme_classic()
ggarrange(a,b)

#now we test the 1st, 2nd and 3rd order polynomial of individual fits with the response variables
#and select the best fit relationship (based on AIC) to build the global model

#FISHING STATUS
Flinear1=lm(lfishingstatus~sfisherdens, data=jurisdictionscale_data)
F2poly1=lm(lfishingstatus~sfisherdens+I(sfisherdens^2), data=jurisdictionscale_data)
F3poly1=lm(lfishingstatus~sfisherdens+I(sfisherdens^2)+I(sfisherdens^3), data=jurisdictionscale_data)
aicff=AIC(Flinear1,F2poly1,F3poly1)
col2=aicff$AIC

Flinear2=lm(lfishingstatus~stotalgravity, data=jurisdictionscale_data)
F2poly2=lm(lfishingstatus~stotalgravity+I(stotalgravity^2), data=jurisdictionscale_data)
F3poly2=lm(lfishingstatus~stotalgravity+I(stotalgravity^2)+I(stotalgravity^3), data=jurisdictionscale_data)
aicff2=AIC(Flinear2,F2poly2,F3poly2)
col3=aicff2$AIC

Flinear3=lm(lfishingstatus~sHDI, data=jurisdictionscale_data)
F2poly3=lm(lfishingstatus~sHDI+I(sHDI^2), data=jurisdictionscale_data)
F3poly3=lm(lfishingstatus~sHDI+I(sHDI^2)+I(sHDI^3), data=jurisdictionscale_data)
aicff3=AIC(Flinear3,F2poly3,F3poly3)
col4=aicff3$AIC

Flinear4=lm(lfishingstatus~stouristdens, data=jurisdictionscale_data)
F2poly4=lm(lfishingstatus~stouristdens+I(stouristdens^2), data=jurisdictionscale_data)
F3poly4=lm(lfishingstatus~stouristdens+I(stouristdens^2)+I(stouristdens^3), data=jurisdictionscale_data)
aicff4=AIC(Flinear4,F2poly4,F3poly4)
col5=aicff4$AIC

Flinear5=lm(lfishingstatus~spopgrowth, data=jurisdictionscale_data)
F2poly5=lm(lfishingstatus~spopgrowth+I(spopgrowth^2), data=jurisdictionscale_data)
F3poly5=lm(lfishingstatus~spopgrowth+I(spopgrowth^2)+I(spopgrowth^3), data=jurisdictionscale_data)
aicff5=AIC(Flinear5,F2poly5,F3poly5)
col6=aicff5$AIC


Flinear6=lm(lfishingstatus~smpa, data=jurisdictionscale_data)
F2poly6=lm(lfishingstatus~smpa+I(smpa^2), data=jurisdictionscale_data)
F3poly6=lm(lfishingstatus~smpa+I(smpa^2)+I(smpa^3), data=jurisdictionscale_data)
aicff6=AIC(Flinear6,F2poly6,F3poly6)
col7=aicff6$AIC

Flinear7=lm(lfishingstatus~smeanttm, data=jurisdictionscale_data)
F2poly7=lm(lfishingstatus~smeanttm+I(smeanttm^2), data=jurisdictionscale_data)
F3poly7=lm(lfishingstatus~smeanttm+I(smeanttm^2)+I(smeanttm^3), data=jurisdictionscale_data)
aicff7=AIC(Flinear7,F2poly7,F3poly7)
col8=aicff7$AIC


#BIOMASS STATUS
Blinear1=lm(lbiomassstatus~sfisherdens, data=jurisdictionscale_data)
B2poly1=lm(lbiomassstatus~sfisherdens+I(sfisherdens^2), data=jurisdictionscale_data)
B3poly1=lm(lbiomassstatus~sfisherdens+I(sfisherdens^2)+I(sfisherdens^3), data=jurisdictionscale_data)
aicbb=AIC(Blinear1,B2poly1,B3poly1)
colb2=aicbb$AIC

Blinear2=lm(lbiomassstatus~stotalgravity, data=jurisdictionscale_data)
B2poly2=lm(lbiomassstatus~stotalgravity+I(stotalgravity^2), data=jurisdictionscale_data)
B3poly2=lm(lbiomassstatus~stotalgravity+I(stotalgravity^2)+I(stotalgravity^3), data=jurisdictionscale_data)
aicbb2=AIC(Blinear2,B2poly2,B3poly2)
colb3=aicbb2$AIC


Blinear3=lm(lbiomassstatus~sHDI, data=jurisdictionscale_data)
B2poly3=lm(lbiomassstatus~sHDI+I(sHDI^2), data=jurisdictionscale_data)
B3poly3=lm(lbiomassstatus~sHDI+I(sHDI^2)+I(sHDI^3), data=jurisdictionscale_data)
aicbb3=AIC(Blinear3,B2poly3,B3poly3)
colb4=aicbb3$AIC

Blinear4=lm(lbiomassstatus~stouristdens, data=jurisdictionscale_data)
B2poly4=lm(lbiomassstatus~stouristdens+I(stouristdens^2), data=jurisdictionscale_data)
B3poly4=lm(lbiomassstatus~stouristdens+I(stouristdens^2)+I(stouristdens^3), data=jurisdictionscale_data)
aicbb4=AIC(Blinear4,B2poly4,B3poly4)
colb5=aicbb4$AIC

Blinear5=lm(lbiomassstatus~spopgrowth, data=jurisdictionscale_data)
B2poly5=lm(lbiomassstatus~spopgrowth+I(spopgrowth^2), data=jurisdictionscale_data)
B3poly5=lm(lbiomassstatus~spopgrowth+I(spopgrowth^2)+I(spopgrowth^3), data=jurisdictionscale_data)
aicbb5=AIC(Blinear5,B2poly5,B3poly5)
colb6=aicbb5$AIC

Blinear6=lm(lbiomassstatus~smpa, data=jurisdictionscale_data)
B2poly6=lm(lbiomassstatus~smpa+I(smpa^2), data=jurisdictionscale_data)
B3poly6=lm(lbiomassstatus~smpa+I(smpa^2)+I(smpa^3), data=jurisdictionscale_data)
aicbb6=AIC(Blinear6,B2poly6,B3poly6)
colb7=aicbb6$AIC

Blinear7=lm(lbiomassstatus~smeanttm, data=jurisdictionscale_data)
B2poly7=lm(lbiomassstatus~smeanttm+I(smeanttm^2), data=jurisdictionscale_data)
B3poly7=lm(lbiomassstatus~smeanttm+I(smeanttm^2)+I(smeanttm^3), data=jurisdictionscale_data)
aicbb7=AIC(Blinear7,B2poly7,B3poly7)
colb8=aicbb7$AIC



#AIC table of biomass and fishing status
col1=c("Linear", "Polynomial (2)","Polynomial (3)")
FAICdata=cbind(col1,col2,col3,col4,col5,col6,col7,col8)
BAICdata=cbind(col1,colb2,colb3,colb4,colb5,colb6,colb7,colb8)

AICDATA=rbind(FAICdata,BAICdata)
colnames(AICDATA)=c("Model","Fisher density","Mean gravity","HDI","Tourist density", "Population growth","%territoial waters protected", "Mean travel time markets")
AICDATA=as.data.frame(AICDATA)
AICDATA$response=rep(c("log(C/MMSY)","log(B/BMMSY)"),each=3)
#write.csv(AICDATA, "best_individualAIC_socioeconomic.csv")

#MODEL FOR FISHING STATUS
jurisdictionscale_data2=jurisdictionscale_data[!is.na(jurisdictionscale_data$lfishingstatus), ]
jurisdictionscale_data2=droplevels(jurisdictionscale_data2)

##global based on best fit individual
global_fstatus=lm(lfishingstatus~sfisherdens_imp+I(sfisherdens_imp^2)+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+I(smpa_imp^2)+I(smpa_imp^3)+smeanttm_imp+I(smeanttm_imp^2), data=jurisdictionscale_data2)
coefplot(global_fstatus)

#simplified model eliminating non significant polinomials
fstatus_model=lm(lfishingstatus~stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=jurisdictionscale_data2)
AIC(global_fstatus,fstatus_model)

#check for influential values
leveragePlots(fstatus_model)
plot(fstatus_model,5)#leverage 
plot(fstatus_model,4)#cooks distance

#eliminate the influential value
jurisdictionscale_data2$area_name[row.names(jurisdictionscale_data2)=="77"]
jurisdictionscale_data4=jurisdictionscale_data2[!jurisdictionscale_data2$area_name=="PRIA",]

#fit model
global_fstatus=lm(lfishingstatus~sfisherdens_imp+I(sfisherdens_imp^2)+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+I(smpa_imp^2)+I(smpa_imp^3)+smeanttm_imp+I(smeanttm_imp^2), data=jurisdictionscale_data4)
coefplot(global_fstatus)
fstatus_model=lm(lfishingstatus~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=jurisdictionscale_data4)
AIC(global_fstatus,fstatus_model)

#difference of effect sizes with and without imputation
fstatus_model_ni=lm(lfishingstatus~stotalgravity+sfisherdens+stouristdens+sHDI+I(sHDI^2)+spopgrowth+smpa+smeanttm, data=jurisdictionscale_data4)
diff_imp_fs=as.data.frame(cbind(c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets"),coef(fstatus_model)[-1],coef(fstatus_model_ni)[-1]))
colnames(diff_imp_fs)=c("Variable","Stnd.effect size imputed","Stnd.effect size non-imputed")
row.names(diff_imp_fs)=NULL
#write.csv(diff_imp_fs,"impvsnonimp_fstatusmodel_socioeconomic.csv",row.names=F)

#model fit
mf_fs1=ggplot(NULL)+geom_histogram(aes(x=resid(fstatus_model)))+xlab("Residuals fishing status model")+theme_classic()
mf_fs2=ggplot(NULL)+geom_point(aes(y=resid(fstatus_model),x=fitted(fstatus_model)))+ylab("Residuals")+xlab("Fitted log(C/MMSY)")+theme_classic()
ggarrange(mf_fs1,mf_fs2, labels=c("a","b"))

#coefficient plot
coefplotfstatusdat=as.data.frame(confint(fstatus_model, level=0.9))
coefplotfstatusdat$mean=fstatus_model$coefficients
coefplotfstatusdat=coefplotfstatusdat[-1,]
coefplotfstatusdat$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplotfstatusdat$sign=ifelse(coefplotfstatusdat$`5 %`<0 & coefplotfstatusdat$`95 %`<0, "negative",ifelse(coefplotfstatusdat$`5 %`>0 & coefplotfstatusdat$`95 %`>0, "positive", "no effect"))

coefplotf=
  ggplot(coefplotfstatusdat,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("Std. effect size")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()


#relationships
MyDataT<-expand.grid(sfisherdens_imp=seq(min(jurisdictionscale_data2$sfisherdens_imp,na.rm=T),max(jurisdictionscale_data2$sfisherdens_imp,na.rm=T),length=length(jurisdictionscale_data2$sfisherdens_imp)),
                     spopgrowth_imp=0,
                     sHDI_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
colnames(jurisdictionscale_data2)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp+
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp)

ed1=ggplot(MyDataT)+ geom_ribbon(aes(x=sfisherdens_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(sfisherdens_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=sfisherdens_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. Fisher density")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)+ylim(c(-4,8))

MyDataT<-expand.grid(sHDI_imp=seq(min(jurisdictionscale_data$sHDI_imp),max(jurisdictionscale_data$sHDI_imp),length=length(jurisdictionscale_data$sHDI_imp)),
                     spopgrowth_imp=0,
                     sfisherdens_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp+
                                                                                           coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp)



ed2=ggplot(MyDataT)+ geom_ribbon(aes(x=sHDI_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(sHDI_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=sHDI_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. HDI")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)


MyDataT<-expand.grid(spopgrowth_imp=seq(min(jurisdictionscale_data$spopgrowth_imp),max(jurisdictionscale_data$spopgrowth_imp),length=length(jurisdictionscale_data$spopgrowth_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp+
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp)


ed3=ggplot(MyDataT)+ geom_ribbon(aes(x=spopgrowth_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(spopgrowth_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=spopgrowth_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. Population growth")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)


MyDataT<-expand.grid(stotalgravity_imp=seq(min(jurisdictionscale_data$stotalgravity_imp),max(jurisdictionscale_data$stotalgravity_imp),length=length(jurisdictionscale_data$stotalgravity_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp+
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp)




ed4=ggplot(MyDataT)+ geom_ribbon(aes(x=stotalgravity_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(stotalgravity_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=stotalgravity_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. Total gravity")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(stouristdens_imp=seq(min(jurisdictionscale_data$stouristdens_imp),max(jurisdictionscale_data$stouristdens_imp),length=length(jurisdictionscale_data$stouristdens_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp +
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp)



ed5=ggplot(MyDataT)+ geom_ribbon(aes(x=stouristdens_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(stouristdens_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=stouristdens_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. Tourism density")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(smeanttm_imp=seq(min(jurisdictionscale_data$smeanttm_imp),max(jurisdictionscale_data$smeanttm_imp),length=length(jurisdictionscale_data$smeanttm_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     smpa_imp=0, 
                     stouristdens_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp +
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[8]*jurisdictionscale_data2$smpa_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp)



ed6=ggplot(MyDataT)+ geom_ribbon(aes(x=smeanttm_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(smeanttm_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=smeanttm_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. Market travel time")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)
MyDataT<-expand.grid(smpa_imp=seq(min(jurisdictionscale_data$smpa_imp),max(jurisdictionscale_data$smpa_imp),length=length(jurisdictionscale_data$smpa_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(fstatus_model)
MyDataT$pseT <- sqrt(diag(X %*% vcov(fstatus_model) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data2$predicted_lfishingstatus=jurisdictionscale_data2$lfishingstatus-(coef(fstatus_model)[2]*jurisdictionscale_data2$stotalgravity_imp+
                                                                                           coef(fstatus_model)[3]*jurisdictionscale_data2$sfisherdens_imp +
                                                                                           coef(fstatus_model)[5]*jurisdictionscale_data2$sHDI_imp+
                                                                                           coef(fstatus_model)[6]*jurisdictionscale_data2$sHDI_imp^2+
                                                                                           coef(fstatus_model)[7]*jurisdictionscale_data2$spopgrowth_imp+
                                                                                           coef(fstatus_model)[9]*jurisdictionscale_data2$smeanttm_imp+
                                                                                           coef(fstatus_model)[4]*jurisdictionscale_data2$stouristdens_imp)



ed7=ggplot(MyDataT)+ geom_ribbon(aes(x=smpa_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(smpa_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data2, aes(x=smpa_imp, y=predicted_lfishingstatus), col="grey",alpha=0.5)+
  xlab("Stnd. % waters protected")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

windows()
fighingrel=ggarrange(ed1,ed2,ed3,ed4,ed5,ed6,ed7,nrow=4,ncol=2, labels=c("a","b","c","d","e","f","g"))
ggarrange(fighingrel,coefplotf, nrow=1,ncol=2, widths=c(3,2), labels=c("","h"))



#alternative: use binomial with overfishing/not overfishing
fstatus_model2=glm(overfishing_dummy~sfisherdens_imp+stotalgravity_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp , data=jurisdictionscale_data4,family=binomial(link = "logit"))


#model fit
mf_fs_d1=ggplot(NULL)+geom_histogram(aes(x=resid(fstatus_model2)))+xlab("Residuals overfishing model")+theme_classic()
mf_fs_d2=ggplot(NULL)+geom_point(aes(y=resid(fstatus_model2),x=fitted(fstatus_model2)))+ylab("Residuals")+xlab("Fitted ")+theme_classic()
ggarrange(mf_fs_d1,mf_fs_d2, labels=c("a","b"))

#coefficient plot
coefplotfstatusdat2=as.data.frame(confint(fstatus_model2,level=0.9))
coefplotfstatusdat2$mean=coef(fstatus_model2)
coefplotfstatusdat2=coefplotfstatusdat2[-1,]
coefplotfstatusdat2$variable=c("Fisher density","Total Gravity",  "Tourism density","HDI", "HDI^2","Population growth","% protected territorial waters","Travel time nearest markets")
coefplotfstatusdat2$sign=ifelse(coefplotfstatusdat2$`5 %`<0 & coefplotfstatusdat2$`95 %` <0, "negative",ifelse(coefplotfstatusdat2$`5 %` >0 & coefplotfstatusdat2$`95 %` >0, "positive", "no effect"))

coefplotf2=
  ggplot(coefplotfstatusdat2,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("log odds ratio")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()

ed_fig_fs=ggarrange(fighingrel,coefplotf,coefplotf2, nrow=1,ncol=3, widths=c(3,1,2), labels=c("","h  ","i"))
annotate_figure(ed_fig_fs,left = text_grob("log(C/MMSY)", rot = 90))

#BIOMASS STATUS
#global model for biomass status vs simplier model
jurisdictionscale_data3=jurisdictionscale_data[!is.na(jurisdictionscale_data$lbiomassstatus),]

global_bstatus=lm(lbiomassstatus~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp, data=jurisdictionscale_data3)
coefplot(global_bstatus)
bstatus_model=lm(lbiomassstatus~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp, data=jurisdictionscale_data3)
AIC(global_bstatus,bstatus_model)


#check for influential values through cook distance
plot(global_bstatus,5)#leverage 
plot(global_bstatus,4)#cooks distance
cutoff <- 4/((nrow(jurisdictionscale_data3)-length(global_bstatus$coefficients)-2)) 
windows()
plot(global_bstatus, which=4, cook.levels=cutoff) #less than 1


#non imputed version
global_bstatus_ni=lm(lbiomassstatus~sfisherdens+stotalgravity+I(stotalgravity^2)+I(stotalgravity^3)+stouristdens+sHDI+spopgrowth+smpa+smeanttm, data=jurisdictionscale_data3)
#difference of effect sizes with and without imputation
diff_imp_bs=as.data.frame(cbind(c("Fisher density","Total Gravity","Total Gravity^2","Total Gravity^3",  "Tourism density","HDI", "Population growth","% protected territorial waters","Travel time nearest markets"),coef(global_bstatus)[-1],coef(global_bstatus_ni)[-1]))
colnames(diff_imp_bs)=c("Variable","Stnd.effect size imputed","Stnd.effect size non-imputed")
row.names(diff_imp_bs)=NULL
#write.csv(diff_imp_bs,"impvsnonimp_bstatusmodel_socioeconomic.csv",row.names=F)


#model fit
mf_bs1=ggplot(NULL)+geom_histogram(aes(x=resid(global_bstatus)))+theme_classic()+xlab ("Residuals biomass status model")
mf_bs2=ggplot(NULL)+geom_point(aes(x=fitted(global_bstatus),y=resid(global_bstatus)))+theme_classic()+ylab ("Residuals")+xlab("Fitted log(B/BMMSY)")
ggarrange(mf_fs1,mf_fs2,mf_bs1,mf_bs2, nrow=2, ncol=2, labels=c("a","b","c","d"))


#coefficient plot
coefplotbstatusdat=as.data.frame(confint(global_bstatus, level=0.9))
coefplotbstatusdat$mean=global_bstatus$coefficients
coefplotbstatusdat=coefplotbstatusdat[-1,]
coefplotbstatusdat$variable=c("Fisher density","Total Gravity","Total Gravity^2","Total Gravity^3","Tourism density","HDI","Population growth", "% protected territorial waters","Travel time nearest markets")
coefplotbstatusdat$sign=ifelse(coefplotbstatusdat$`5 %`<0 & coefplotbstatusdat$`95 %`<0, "negative",ifelse(coefplotbstatusdat$`5 %`>0 & coefplotbstatusdat$`95 %`>0, "positive", "no effect"))

coefplotb=
  ggplot(coefplotbstatusdat,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("Std. effect size")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(axis.text.y =element_blank() ,axis.text.x = element_text(angle = 0, hjust = 0),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()


#relationships
MyDataT<-expand.grid(sfisherdens_imp=seq(min(jurisdictionscale_data3$sfisherdens_imp,na.rm=T),max(jurisdictionscale_data3$sfisherdens_imp,na.rm=T),length=length(jurisdictionscale_data3$sfisherdens_imp)),
                     spopgrowth_imp=0,
                     sHDI_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
colnames(jurisdictionscale_data3)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                           coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                           coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp+
                                                                                           coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                           coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)
plot(jurisdictionscale_data3$sfisherdens_imp ,jurisdictionscale_data3$predicted_lbiomassstatus)

ed1=ggplot(MyDataT)+ geom_ribbon(aes(x=sfisherdens_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(sfisherdens_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=sfisherdens_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. Fisher density")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(sHDI_imp=seq(min(jurisdictionscale_data3$sHDI_imp),max(jurisdictionscale_data3$sHDI_imp),length=length(jurisdictionscale_data3$sHDI_imp)),
                     spopgrowth_imp=0,
                     sfisherdens_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                           coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                           coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp+
                                                                                           coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                           coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)

ed2=ggplot(MyDataT)+ geom_ribbon(aes(x=sHDI_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(sHDI_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=sHDI_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. HDI")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(spopgrowth_imp=seq(min(jurisdictionscale_data3$spopgrowth_imp),max(jurisdictionscale_data3$spopgrowth_imp),length=length(jurisdictionscale_data3$spopgrowth_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                           coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                           coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                           coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                           coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)

ed3=ggplot(MyDataT)+ geom_ribbon(aes(x=spopgrowth_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(spopgrowth_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=spopgrowth_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. Population growth")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)


MyDataT<-expand.grid(stotalgravity_imp=seq(min(jurisdictionscale_data3$stotalgravity_imp),max(jurisdictionscale_data3$stotalgravity_imp),length=length(jurisdictionscale_data3$stotalgravity_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stouristdens_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-( coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                            coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                            coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                            coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp+
                                                                                            coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                            coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)





ed4=ggplot(MyDataT)+ geom_ribbon(aes(x=stotalgravity_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(stotalgravity_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=stotalgravity_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. Total gravity")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)


MyDataT<-expand.grid(stouristdens_imp=seq(min(jurisdictionscale_data3$stouristdens_imp),max(jurisdictionscale_data3$stouristdens_imp),length=length(jurisdictionscale_data3$stouristdens_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     smpa_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp+
                                                                                           coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                           coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                           coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                           coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)


ed5=ggplot(MyDataT)+ geom_ribbon(aes(x=stouristdens_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(stouristdens_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=stouristdens_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. Tourism density")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(smeanttm_imp=seq(min(jurisdictionscale_data3$smeanttm_imp),max(jurisdictionscale_data3$smeanttm_imp),length=length(jurisdictionscale_data3$smeanttm_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     smpa_imp=0, 
                     stouristdens_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                           coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                           coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                           coef(global_bstatus)[9]*jurisdictionscale_data3$smpa_imp+
                                                                                           coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp)

ed6=ggplot(MyDataT)+ geom_ribbon(aes(x=smeanttm_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(smeanttm_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=smeanttm_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. Market travel time")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

MyDataT<-expand.grid(smpa_imp=seq(min(jurisdictionscale_data3$smpa_imp),max(jurisdictionscale_data3$smpa_imp),length=length(jurisdictionscale_data3$smpa_imp)),
                     sHDI_imp=0,
                     sfisherdens_imp=0,
                     spopgrowth_imp=0,
                     stotalgravity_imp=0,
                     stouristdens_imp=0, 
                     smeanttm_imp=0)
X <- model.matrix(~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp,data=MyDataT)

MyDataT$etaT <- X %*% coef(global_bstatus)
MyDataT$pseT <- sqrt(diag(X %*% vcov(global_bstatus) %*% t(X)))
MyDataT$ploT <- MyDataT$etaT - 1.96 * MyDataT$pseT
MyDataT$phiT <- MyDataT$etaT + 1.96 * MyDataT$pseT

jurisdictionscale_data3$predicted_lbiomassstatus=jurisdictionscale_data3$lbiomassstatus-(coef(global_bstatus)[3]*jurisdictionscale_data3$stotalgravity_imp+
                                                                                           coef(global_bstatus)[4]*jurisdictionscale_data3$stotalgravity_imp^2+
                                                                                           coef(global_bstatus)[5]*jurisdictionscale_data3$stotalgravity_imp^3+
                                                                                           coef(global_bstatus)[6]*jurisdictionscale_data3$stouristdens_imp+
                                                                                           coef(global_bstatus)[2]*jurisdictionscale_data3$sfisherdens_imp+
                                                                                           coef(global_bstatus)[7]*jurisdictionscale_data3$sHDI_imp+
                                                                                           coef(global_bstatus)[8]*jurisdictionscale_data3$spopgrowth_imp+
                                                                                           coef(global_bstatus)[10]*jurisdictionscale_data3$smeanttm_imp)


ed7=ggplot(MyDataT)+ geom_ribbon(aes(x=smpa_imp,ymax=phiT,ymin=ploT),fill="dodgerblue3",alpha=0.2)+
  geom_line(aes(smpa_imp,etaT),colour="dodgerblue3",size=1)+
  geom_point(data=jurisdictionscale_data3, aes(x=smpa_imp, y=predicted_lbiomassstatus),col="grey",alpha=0.5)+
  xlab("Stnd. % waters protected")+
  ylab("")+theme_classic()+geom_hline(yintercept=0, lty=2)

windows()
biomassrel=ggarrange(ed1,ed2,ed3,ed4,ed5,ed6,ed7,nrow=4,ncol=2, labels=c("a","b","c","d","e","f","g"))
ggarrange(biomassrel,coefplotb, nrow=1,ncol=2, widths=c(3,1.5), labels=c("","h"))


#now binomial mode #dummy
global_bstatus2=glm(overfished_dummy~sfisherdens_imp+stotalgravity_imp+I(stotalgravity_imp^2)+I(stotalgravity_imp^3)+stouristdens_imp+sHDI_imp+spopgrowth_imp+smpa_imp+smeanttm_imp, data=jurisdictionscale_data3,family=binomial(link = "logit"))

coefplotbstatusdat2=as.data.frame(confint(global_bstatus2, level=0.9))
coefplotbstatusdat2$mean=coef(global_bstatus2)
coefplotbstatusdat2=coefplotbstatusdat2[-1,]
coefplotbstatusdat2$variable=c("Fisher density","Total Gravity","Total Gravity^2","Total Gravity^3","Tourism density","HDI","Population growth", "% protected territorial waters","Travel time nearest markets")
coefplotbstatusdat2$sign=ifelse(coefplotbstatusdat2$`5 %`<0 & coefplotbstatusdat2$`95 %` <0, "negative",ifelse(coefplotbstatusdat2$`5 %` >0 & coefplotbstatusdat2$`95 %` >0, "positive", "no effect"))

coefplotb2=
  ggplot(coefplotbstatusdat2,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("log odds ratio")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()

ed_fig_bs=ggarrange(biomassrel,coefplotb,coefplotb2, nrow=1,ncol=3, widths=c(3,1,2), labels=c("","h","i"))
annotate_figure(ed_fig_bs,left = text_grob("log(B/BMMSY)", rot = 90))


#SHIFTING BASELINES---------------------------------------------------------------------------------------------------------

#using different travel time cut off to define remoteness
remote_complete_10h=data_benchmark[data_benchmark$remote_10h==1,]

#using 10 hours instead of 20h
stanDat_sb_10h <- list(p=0,b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                       N = nrow(reserves_complete),
                       br=log(remote_complete_10h$FamBiomass_tkm2),
                       RE=nrow(remote_complete_10h),
                       newage=newdata$Closure.age,
                       Bio=Bio, B=B)
Fit_null_sb_10h <- stan(file = "Null_exports_schaeffer.stan", data = stanDat_sb_10h , iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_sb_10h, pars = c("r","bmin","B0","lp__","sigma_e"))
Fit_sb_10h <- summary(Fit_null_sb_10h,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_sb_10h=as.data.frame(Fit_sb_10h$summary)
list_of_draws_sb_10h <- as.data.frame(Fit_null_sb_10h)
newpred_10h=output_Fit_sb_10h[200:269,]
colnames(newpred_10h)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
newpred_10h$Closure.age=newdata$Closure.age
popgrowth_10h=output_Fit_sb_10h[391:741,]
colnames(popgrowth_10h)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
popgrowth_10h$mean_estimate=ifelse(popgrowth_10h$mean_estimate<0,0,popgrowth_10h$mean_estimate)
popgrowth_10h$median=ifelse(popgrowth_10h$median<0,NA,popgrowth_10h$median)
popgrowth_10h$Biomass=Bio
popgrowth_10h$`5%`=ifelse(popgrowth_10h$`5%`<0,0,popgrowth_10h$`5%`)
popgrowth_10h$`95%`=ifelse(popgrowth_10h$`95%`<0,0,popgrowth_10h$`95%`)
popgrowth_10h$`25%`=ifelse(popgrowth_10h$`25%`<0,0,popgrowth_10h$`25%`)
popgrowth_10h$`75%`=ifelse(popgrowth_10h$`75%`<0,0,popgrowth_10h$`75%`)

#fixing b0 to different lower values
stanDat_null3 <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                      N = nrow(reserves_complete),
                      B0=100,RE=nrow(remote_complete),
                      newage=newdata$Closure.age,
                      Bio=Bio, B=B)
Fit_null3 <- stan(file = "fixed_lowerB0.stan", data = stanDat_null3,iter=10000,warmup=9000, chains = 4, control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
stanDat_null4 <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                      N = nrow(reserves_complete),
                      B0=50,RE=nrow(remote_complete),
                      newage=newdata$Closure.age,
                      Bio=Bio, B=B)
Fit_null4 <- stan(file = "fixed_lowerB0.stan", data = stanDat_null4, iter=10000,warmup=9000,chains = 4, control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null4, pars = c("r","bmin","lp__","sigma_e"))
stan_hist(Fit_null4, pars = c("MMSY"))
Fit_reserves_summary3 <- summary(Fit_null3,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_reserves3=as.data.frame(Fit_reserves_summary3$summary)
Fit_reserves_summary4 <- summary(Fit_null4,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_reserves4=as.data.frame(Fit_reserves_summary4$summary)
list_of_draws_reserves3 <- as.data.frame(Fit_null3) 
list_of_draws_reserves4 <- as.data.frame(Fit_null4) 

#plot surplus production curve along a gradient of biomass
popgrowth_reserves3=output_Fit_reserves3[146:496,]
colnames(popgrowth_reserves3)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median","75%","95%","n_eff","Rhat")
popgrowth_reserves3$mean_estimate=ifelse(popgrowth_reserves3$mean_estimate<0,-NA,popgrowth_reserves3$mean_estimate)
popgrowth_reserves3$median=ifelse(popgrowth_reserves3$median<0,-NA,popgrowth_reserves3$median)
popgrowth_reserves3$Biomass=Bio
popgrowth_reserves3$`5%`=ifelse(popgrowth_reserves3$`5%`<0,0,popgrowth_reserves3$`5%`)
popgrowth_reserves3$`95%`=ifelse(popgrowth_reserves3$`95%`<0,0,popgrowth_reserves3$`95%`)
popgrowth_reserves3$`25%`=ifelse(popgrowth_reserves3$`25%`<0,0,popgrowth_reserves3$`25%`)
popgrowth_reserves3$`75%`=ifelse(popgrowth_reserves3$`75%`<0,0,popgrowth_reserves3$`75%`)
newpred_reserves3=output_Fit_reserves3[76:145,]
colnames(newpred_reserves3)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
newpred_reserves3$Closure.age=newdata$Closure.age

popgrowth_reserves4=output_Fit_reserves4[146:496,]
colnames(popgrowth_reserves4)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median","75%","95%","n_eff","Rhat")
popgrowth_reserves4$mean_estimate=ifelse(popgrowth_reserves4$mean_estimate<0,-NA,popgrowth_reserves4$mean_estimate)
popgrowth_reserves4$median=ifelse(popgrowth_reserves4$median<0,-NA,popgrowth_reserves4$median)
popgrowth_reserves4$Biomass=Bio
popgrowth_reserves4$`5%`=ifelse(popgrowth_reserves4$`5%`<0,0,popgrowth_reserves4$`5%`)
popgrowth_reserves4$`95%`=ifelse(popgrowth_reserves4$`95%`<0,0,popgrowth_reserves4$`95%`)
popgrowth_reserves4$`25%`=ifelse(popgrowth_reserves4$`25%`<0,0,popgrowth_reserves4$`25%`)
popgrowth_reserves4$`75%`=ifelse(popgrowth_reserves4$`75%`<0,0,popgrowth_reserves4$`75%`)
newpred_reserves4=output_Fit_reserves4[76:145,]
colnames(newpred_reserves4)=c("mean_estimate","se_mean","sd_estimate","5%","25%","median", "75%","95%","n_eff","Rhat")
newpred_reserves4$Closure.age=newdata$Closure.age

#plot shifting baselines
null_grow_sb=ggplot(NULL)+geom_line(aes(y=popgrowth_null2$median[!is.na(popgrowth_null2$median)], x=popgrowth_null2$Biomass[!is.na(popgrowth_null2$median)]), col="black", lwd=2)+
  geom_line(aes(y=popgrowth_10h$median[!is.na(popgrowth_10h$median)], x=popgrowth_10h$Biomass[!is.na(popgrowth_10h$median)]), col="darkmagenta", lwd=2)+
  geom_line(aes(y=popgrowth_reserves4$median[!is.na(popgrowth_reserves4$median)], x=popgrowth_reserves4$Biomass[!is.na(popgrowth_reserves4$median)]), col="red", lwd=2)+
  geom_line(aes(y=popgrowth_reserves3$median[!is.na(popgrowth_reserves3$median)], x=popgrowth_reserves3$Biomass[!is.na(popgrowth_reserves3$median)]), col="darkred", lwd=2)+
  theme_classic()+xlab("Biomass (t/km2)")+ylab("Potential sustainable yield (t/km2/y)")+ labs(x = expression ("Biomass ("~t/km^2*")"),y = expression ("Potential sustainable yield ("~t/km^2/y*")"))+ylim(c(0,6))

null_rdist_sb=ggplot(NULL)+geom_density(aes(x=list_of_draws_null$r), col="black",fill="black",alpha=0.5)+
  geom_density(aes(x=list_of_draws_sb_10h$r), col="black",fill="darkmagenta",alpha=0.5)+
  geom_density(aes(x=list_of_draws_reserves3$r), col="black",fill="darkred",alpha=0.5)+
  geom_density(aes(x=list_of_draws_reserves4$r), col="black",fill="red",alpha=0.5)+
  theme_classic()+xlim(c(0,0.5))+ylab("Posterior density")+xlab("Biomass intrinsic growth rate")

ggarrange(null_rdist_sb,null_grow_sb, nrow=1, ncol=2, labels=c("a","b"))

#B0-r-MMSY trade-off
r_var=c(median(list_of_draws_null$r),median(list_of_draws_sb_10h$r), median(list_of_draws_reserves3$r), median(list_of_draws_reserves4$r))
b0_var=c(median(list_of_draws_null$B0),median(list_of_draws_sb_10h$B0), 100, 50)
MMSY_var=c(median(list_of_draws_null$MMSY),median(list_of_draws_sb_10h$MMSY), median(list_of_draws_reserves3$MMSY), median(list_of_draws_reserves4$MMSY))
tradeoff_data=as.data.frame(cbind(r_var,b0_var,MMSY_var))


#do it continuously
B0_range=seq(50,250,5)
MMSY_range=rep(NA, length(B0_range))
MMSY_ci=matrix(NA, nrow=length(B0_range), ncol=4)
r_range=rep(NA, length(B0_range))
r_ci=matrix(NA, nrow=length(B0_range), ncol=4)
for (i in 1:length(B0_range)){
  stanDat_null_range <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                             N = nrow(reserves_complete),
                             B0=B0_range[i],RE=nrow(remote_complete),
                             newage=newdata$Closure.age,
                             Bio=Bio, B=B)
  Fit_null_range <- stan(file = "fixed_lowerB0.stan", data = stanDat_null_range, iter=5000,chains = 4, control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
  
  list_of_draws_range <- as.data.frame(Fit_null_range) 
  MMSY_range[i]=median(list_of_draws_range$MMSY)
  MMSY_ci[i,]=quantile(list_of_draws_range$MMSY,  probs=c(0.05,0.25,0.75, 0.95))
  r_range[i]=median(list_of_draws_range$r)
  r_ci[i,]=quantile(list_of_draws_range$r,  probs=c(0.05,0.25,0.75, 0.95))
}
windows()
par(mar=c(5, 4, 4, 6) + 0.1)
plot(log(B0_range), MMSY_range, pch=16, axes=F, ylim=c(1,6), xlab="", ylab="", 
     type="l",col="black")

polygon(c(rev(log(B0_range)), log(B0_range)),c(rev(MMSY_ci[,2]),MMSY_ci[,3]), col=adjustcolor("grey",alpha=0.5), border=NA)
par(new=TRUE)
plot(log(B0_range), MMSY_range, pch=16, axes=F, ylim=c(1,6), xlab="", ylab="", 
     type="o",col="black")

axis(2, ylim=c(1,6),col="black",las=1)  ## las=1 makes horizontal labels
mtext("MMSY",side=2,line=3)
box()

## Allow a second plot on the same graph
par(new=TRUE)
plot(log(B0_range), r_range, pch=15,  xlab="", ylab="", ylim=c(0,0.5), 
     axes=F, type="l", col="navyblue")

## Plot the second plot and put axis scale on right
polygon(c(rev(log(B0_range)), log(B0_range)),c(rev(r_ci[,2]),r_ci[,3]), col=adjustcolor("navyblue",alpha=0.3), border=NA)
par(new=TRUE)
plot(log(B0_range), r_range, pch=15,  xlab="", ylab="", ylim=c(0,0.5), 
     axes=FALSE, type="o", col="navyblue")

## a little farther out (line=4) to make room for labels
mtext("Intrinsic growth rate",side=4,col="navyblue",line=3) 
axis(4, ylim=c(0,0.5), col="navyblue",col.axis="navyblue",las=1)

## Draw the time axis
#axis(1,at=pretty(range(log(B0_range),10)), xlim=c(3,10))
axis(1,at=c(4,4.5,5,5.5), labels=c(55,90,148,245))
mtext(expression(paste( plain("Unfished biomass (t/km") ^ plain("2"), plain(" )") )),side=1,col="black",line=2.5)  




###########
#Testing representativeness of sampled reefs...............................................................................................
#Are our sampled reefs representative of their jurisdiction?
#Testing potential bias towards more/less accesible reefs by looking at gravity

#calculate median gravity by country 
Tbycountry=ddply(alldata,.(Larger),summarize,median_ttmarket=median(TT_market_h,na.rm=T), median_ttpop=median(TT_pop_h,na.rm=T), mean_ttmarket=mean(TT_market_h,na.rm=T), mean_ttpop=mean(TT_pop_h,na.rm=T), mean_grav=mean(Grav_tot, na.rm=T), median_grav=median(Grav_tot, na.rm=T), median_Bmarg=median(Biomass_marg), mean_Bmarg=mean(Biomass_marg))
tt=merge(jurisdictionscale_data,Tbycountry,  by="Larger", all.y=T)
windows()
ggplot(tt)+geom_point(aes(y=log(mean_grav+1),x=log(meanTgrav_nh2+1)))+geom_abline(intercept=0, slope=1, lty=2)+geom_text_repel(aes(y=log(mean_grav+1),x=log(meanTgrav_nh2+1),label=tt$Larger))+ylab("log(mean gravity of sampled reefs+1)")+xlab("log(mean gravity of entire jurisdiction+1)")+geom_smooth(aes(y=log(mean_grav+1),x=log(meanTgrav_nh2+1)),method="lm",  level=0.95,col="navyblue",  alpha=0.5)+theme_classic()

#check if intercept overlaps 0 and slope overlaps 1
model=stan_glm(log(mean_grav+1)~log(meanTgrav_nh2+1), data=tt, family=gaussian)
pairs(model)


################################################################################################################################
#SENSITIVITY ANALYSES

#.........................................................................................................
#1: proof that including remote locations improved the certainty of estimates

#reserves only with graham-schaefer: Prove that using remote places informs the b0 and doesnt allow for high likely implausible values
stanDat_reserves <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                         N = nrow(reserves_complete),
                         newage=newdata$Closure.age,
                         Bio=Bio, B=B)
Fit_null_reserves <- stan(file = "Null_reserves_schaeffer.stan", data = stanDat_reserves, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))


windows()
pairs(Fit_null_reserves, pars = c("r","bmin","B0","lp__"))
list_of_draws_null_reserves <- as.data.frame(Fit_null_reserves) 

#quick check: B0 from reserve only submodel is higher that prior B0
#is it sensitive to prior b0?
Fit_null_reserves_priorsens <- stan(file = "Null_reserves_schaeffer_priorsensitivity.stan", data = stanDat_reserves, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
windows()
pairs(Fit_null_reserves_priorsens, pars = c("r","bmin","B0","lp__"))

#two step approach (suggested by reviewer): first fit a model to remote only and then use the posterior of that as a prior for the reserve only model
#note that this should give almost identical results to joint model given probability theory

#remote only model for B0
stanDat_remote<-list(br=log(remote_complete$FamBiomass_tkm2),
                     RE=nrow(remote_complete))
Fit_remoteonly <- stan(file = "Null_remote.stan", data = stanDat_remote,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))

#extract B0 to use as a prior for reserves only model
list_of_draws_null_remote <- as.data.frame(Fit_remoteonly) #rows are iterations and columns are parameters
ggplot(list_of_draws_null_remote, aes(x=B0)) +
  geom_density( alpha=0.5)+
  theme_classic()+ylab("Posterior density")+ labs(x =  expression ("Unfished Biomass ("~t/km^2*")"))

#now we use the mean and sd of this posterior as a orior for the reserve only submodel
stanDat_reserves2 <- list(sigma_B0=sd(list_of_draws_null_remote$B0), mean_B0=mean(list_of_draws_null_remote$B0), b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                          N = nrow(reserves_complete),
                          newage=newdata$Closure.age,
                          Bio=Bio, B=B)
Fit_null_reserves2 <- stan(file = "Null_reservesonly_remoteprior.stan", data = stanDat_reserves2, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
list_of_draws_null_reserves2 <- as.data.frame(Fit_null_reserves2) #rows are iterations and columns are parameters


#compare estimates
median(list_of_draws_null$B0)
median(list_of_draws_null_reserves$B0)
mean(list_of_draws_null$B0)
mean(list_of_draws_null_reserves$B0)
max(list_of_draws_null_reserves$B0)
max(list_of_draws_null$B0)
comparisondata=as.data.frame(cbind(c(list_of_draws_null$MMSY,list_of_draws_null_reserves$MMSY,list_of_draws_null_reserves2$MMSY),c(list_of_draws_null$BMMSY,list_of_draws_null_reserves$BMMSY,list_of_draws_null_reserves2$BMMSY),c(list_of_draws_null$r,list_of_draws_null_reserves$r,list_of_draws_null_reserves2$r), c(list_of_draws_null$B0,list_of_draws_null_reserves$B0,list_of_draws_null_reserves2$B0),c(list_of_draws_null$bmin,list_of_draws_null_reserves$bmin,list_of_draws_null_reserves2$bmin),c(rep("reserves_remote", length.out=length(list_of_draws_null$r)),rep("reserves_only", length.out=length(list_of_draws_null_reserves$r)),rep("reserves_remote prior", length.out=length(list_of_draws_null_reserves2$r)))))
colnames(comparisondata)=c("MMSY","BMMSY","r","B0","bmin","model")
comparisondata$r=as.numeric(as.character(comparisondata$r))
comparisondata$B0=as.numeric(as.character(comparisondata$B0))
comparisondata$bmin=as.numeric(as.character(comparisondata$bmin))
comparisondata$MMSY=as.numeric(as.character(comparisondata$MMSY))
comparisondata$BMMSY=as.numeric(as.character(comparisondata$BMMSY))
comparisondata$model=relevel(comparisondata$model, ref="reserves_remote")
a=ggplot(comparisondata, aes(x=bmin,  group=model)) +
  geom_density( aes(fill=comparisondata$model),lwd=1, alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("Posterior density")+ labs(x =  expression ("Biomass reserve age 0 ("~t/km^2*")"))
b=ggplot(comparisondata, aes(x=B0,  group=model)) +
  geom_density( aes(fill=comparisondata$model), lwd=1, alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+xlab("B0 (t/km2)")+ylab("")+ labs(x =  expression ("Unfished Biomass ("~t/km^2*")"))
c=ggplot(comparisondata, aes(x=r,  group=model)) +guides(lty=F,fill=FALSE)+
  geom_density( aes(fill=comparisondata$model), lwd=1, alpha=0.5)+
  labs(fill="Model")+theme_classic()+xlab("intrinsic growth rate")+ylab("")
d=ggplot(comparisondata, aes(x=MMSY,  group=model)) +
  geom_density( aes(fill=comparisondata$model), lwd=1, alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))
e=ggplot(comparisondata, aes(x=BMMSY,  group=model)) +
  geom_density( aes(fill=comparisondata$model), lwd=1, alpha=0.5)  + labs(fill="Model")+
  theme_classic()+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"))
windows()
ggarrange(a,b,c,d,e, nrow=1,ncol=5,widths=c(1,1,1,1,1.8) )

#status results with reserve only submodel
MMSY_reserves=median(list_of_draws_null_reserves$MMSY)
BMMSY_reserves=median(list_of_draws_null_reserves$BMMSY)
B0_reserves=median(list_of_draws_null_reserves$B0)

#median surplus production curve for  reserves
Fit_null_reserves_summary=summary(Fit_null_reserves)
output_Fit_null_reserves=as.data.frame(Fit_null_reserves_summary$summary)
popgrowth_null_reserves=output_Fit_null_reserves[148:498,c("50%")]
popgrowth_null_reserves=ifelse(popgrowth_null_reserves<0,0.0001,popgrowth_null_reserves)
popgrowth_null_reserves=as.data.frame(popgrowth_null_reserves)
colnames(popgrowth_null_reserves)=c("reserves_surplus")
popgrowth_null_reserves$Biomass=Bio
jurisdictionscale_data$Biomass=round(jurisdictionscale_data$medianB_tkm2)
popgrowth_null_reserves2=merge(jurisdictionscale_data,popgrowth_null_reserves,by="Biomass",all.x=T)
popgrowth_null_reserves2$biomassstatus_curve=popgrowth_null_reserves2$medianB_tkm2/BMMSY_reserves
popgrowth_null_reserves2$collapsed=ifelse(is.na(popgrowth_null_reserves2$medianB_tkm2),NA,ifelse(popgrowth_null_reserves2$medianB_tkm2<(0.1*B0_reserves),1,0))

#jurisdictions fishing and fishery status: 
popgrowth_null_reserves2$fishingstatus=popgrowth_null_reserves2$catch_tkm2/MMSY_reserves
popgrowth_null_reserves2$overfishing=ifelse(is.na(popgrowth_null_reserves2$fishingstatus),NA,ifelse(popgrowth_null_reserves2$fishingstatus>1,"Overfishing","Not overfishing"))
popgrowth_null_reserves2$fishingstatus_curve=popgrowth_null_reserves2$catch_tkm2/popgrowth_null_reserves2$reserves_surplus
popgrowth_null_reserves2$overfishing_curve=ifelse(is.na(popgrowth_null_reserves2$fishingstatus_curve),NA,ifelse(popgrowth_null_reserves2$fishingstatus_curve>1,1,0))
popgrowth_null_reserves2$fisherystatus=as.character(ifelse(is.na(popgrowth_null_reserves2$fishingstatus_curve),NA, ifelse(popgrowth_null_reserves2$fishingstatus_curve<1 & popgrowth_null_reserves2$biomassstatus_curve>1, "Sustainable", ifelse(popgrowth_null_reserves2$fishingstatus_curve>1 & popgrowth_null_reserves2$biomassstatus_curve<1, "Unsustainable", ifelse(popgrowth_null_reserves2$fishingstatus_curve>1 & popgrowth_null_reserves2$biomassstatus_curve>1, "Warning", "Rebuilding")))))
#change quasisustainable ones
popgrowth_null_reserves2$fisherystatus=as.character(ifelse(is.na(popgrowth_null_reserves2$fishingstatus_curve),NA, ifelse(popgrowth_null_reserves2$fishingstatus_curve<1 & popgrowth_null_reserves2$biomassstatus_curve>1, "Sustainable", ifelse(popgrowth_null_reserves2$fishingstatus_curve>1 & popgrowth_null_reserves2$biomassstatus_curve<1, "Unsustainable", ifelse(popgrowth_null_reserves2$fishingstatus_curve>1 & popgrowth_null_reserves2$biomassstatus_curve>1, "Warning", "Rebuilding")))))
popgrowth_null_reserves2$fisherystatus=ifelse(is.na(popgrowth_null_reserves2$fisherystatus),NA,ifelse(popgrowth_null_reserves2$fisherystatus=="Warning" &popgrowth_null_reserves2$fishingstatus<1, "Sustainable",as.character(popgrowth_null_reserves2$fisherystatus)))

#save results
reserves_status_results=as.data.frame(c(length(popgrowth_null_reserves2$fishingstatus[!is.na(popgrowth_null_reserves2$fishingstatus)]),
                                        length(popgrowth_null_reserves2$overfishing[!is.na(popgrowth_null_reserves2$overfishing)&popgrowth_null_reserves2$overfishing=="Overfishing"])/length(popgrowth_null_reserves2$overfishing[!is.na(popgrowth_null_reserves2$overfishing)]),
                                        length(popgrowth_null_reserves2$collapsed[!is.na(popgrowth_null_reserves2$collapsed)]),
                                        length(popgrowth_null_reserves2$collapsed[!is.na(popgrowth_null_reserves2$collapsed)&popgrowth_null_reserves2$collapsed==1])/length(popgrowth_null_reserves2$collapsed[!is.na(popgrowth_null_reserves2$collapsed)]),
                                        length(popgrowth_null_reserves2$biomassstatus_curve[!is.na(popgrowth_null_reserves2$biomassstatus_curve)&popgrowth_null_reserves2$biomassstatus_curve<1])/length(popgrowth_null_reserves2$biomassstatus_curve[!is.na(popgrowth_null_reserves2$biomassstatus_curve)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)]),
                                        length(popgrowth_null_reserves2$overfishing_curve[!is.na(popgrowth_null_reserves2$overfishing_curve)&popgrowth_null_reserves2$overfishing_curve>0.9])/length(popgrowth_null_reserves2$overfishing_curve[!is.na(popgrowth_null_reserves2$overfishing_curve)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)&popgrowth_null_reserves2$fisherystatus=="Sustainable"])/length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)&popgrowth_null_reserves2$fisherystatus=="Unsustainable"])/length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)&popgrowth_null_reserves2$fisherystatus=="Warning"])/length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)&popgrowth_null_reserves2$fisherystatus=="Rebuilding"])/length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)]),
                                        length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)&!popgrowth_null_reserves2$fisherystatus=="Sustainable"])/length(popgrowth_null_reserves2$fisherystatus[!is.na(popgrowth_null_reserves2$fisherystatus)])))

colnames(reserves_status_results)=c("reserve_only")
row.names(reserves_status_results)=c("n_fishingstatus","overfishing_MMSY","n_biomassstatus","collapsed","overfished","n_fisherystatus","overfishing_curve", "sustainable","unsustainable","warning","rebuilding","conservation_concern")
#write.csv(reserves_status_results,"sensitivity_status_reserveonly_quasi.csv")


#.....................................................................................................
#2: Proof estimated bmin is representative of degraded systems

#model bounding bmin with high gravity fished locations: above the 0.75 quantile fo the gravity distribution
hist(log(reefscale_data$Grav_tot+1))
reefscale_data$lGrav_tot=log(reefscale_data$Grav_tot+1)
qgrav=quantile(reefscale_data$lGrav_tot, 0.75)
fishedbmin=reefscale_data[reefscale_data$Management=="Fished"& reefscale_data$lGrav_tot>qgrav,]
stanDat_bminb1 <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                       N = nrow(reserves_complete),
                       br=log(remote_complete$FamBiomass_tkm2),
                       RE=nrow(remote_complete),
                       newage=newdata$Closure.age,
                       Bio=Bio, B=B,
                       FI=nrow(fishedbmin),
                       bf=log(fishedbmin$FamBiomass_tkm2))
Fit_bminb <- stan(file = "Null_bminbounded.stan", data = stanDat_bminb1, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20), init=list(list(bmin=1.1),list(bmin=10),list(bmin=5), list(bmin=12)))
windows()
pairs(Fit_bminb , pars = c("r","bmin","B0","lp__"))
list_of_draws_null_boundedbmin <- as.data.frame(Fit_bminb) 

comparisondata2=as.data.frame(cbind(c(list_of_draws_null$MMSY,list_of_draws_null_boundedbmin$MMSY),c(list_of_draws_null$BMMSY,list_of_draws_null_boundedbmin$BMMSY),c(list_of_draws_null$r,list_of_draws_null_boundedbmin$r), c(list_of_draws_null$B0,list_of_draws_null_boundedbmin$B0),c(list_of_draws_null$bmin,list_of_draws_null_boundedbmin$bmin),c(rep("reserves_remote", length.out=length(list_of_draws_null$r)),rep("bmin_fished", length.out=length(list_of_draws_null_boundedbmin$r)))))
colnames(comparisondata2)=c("MMSY","BMMSY","r","B0","bmin","model")
comparisondata2$r=as.numeric(as.character(comparisondata2$r))
comparisondata2$B0=as.numeric(as.character(comparisondata2$B0))
comparisondata2$bmin=as.numeric(as.character(comparisondata2$bmin))
comparisondata2$MMSY=as.numeric(as.character(comparisondata2$MMSY))
comparisondata2$BMMSY=as.numeric(as.character(comparisondata2$BMMSY))

a=ggplot(comparisondata2, aes(x=bmin,  group=model)) +
  geom_density( aes(fill=comparisondata2$model), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("Posterior density")+ labs(x =  expression ("Biomass reserve age 0 ("~t/km^2*")"))
b=ggplot(comparisondata2, aes(x=B0,  group=model)) +
  geom_density( aes(fill=comparisondata2$model), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+xlab("B0 (t/km2)")+ylab("")+ labs(x =  expression ("Unfished Biomass ("~t/km^2*")"))
c=ggplot(comparisondata2, aes(x=r,  group=model)) +
  geom_density( aes(fill=comparisondata2$model), alpha=0.5)+guides(fill=FALSE)+
  labs(fill="Model")+theme_classic()+xlab("intrinsic growth rate")+ylab("")
d=ggplot(comparisondata2, aes(x=MMSY,  group=model)) +
  geom_density( aes(fill=comparisondata2$model), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))+xlim(c(0,10))
e=ggplot(comparisondata2, aes(x=BMMSY,  group=model)) +
  geom_density( aes(fill=comparisondata2$model), alpha=0.5)  +labs(fill="Model")+
  theme_classic()+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"))
ggarrange(a,b, c,d,e,widths=c(1,1,1,1,1.5), nrow=1, ncol=5)

#............................................................................................
#3.  Full model modelling bmin also as a function of covariates

Fit_full_bminalsocovariates <- stan(file = "Full_exports_schaeffer_noncentered_bmin.stan", data = stanDat_full,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
windows()
pairs(Fit_full_bminalsocovariates, pars = c("r","bmin","B0","lp__"))
Fit_bminalsocovariates_summary <- summary(Fit_full_bminalsocovariates,probs=c(0.05,0.25,0.5,0.75, 0.95))
list_of_draws_null_bmincovariates=as.data.frame(Fit_full_bminalsocovariates)
a=ggplot(NULL)+geom_density(aes(x=list_of_draws_null$MMSY), fill=NA,lty=2,lwd=2)+geom_density(aes(x=list_of_draws_null_bmincovariates$MMSY), fill="green",alpha=0.5)+theme_classic()+ylab("Posterior density")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))
b=ggplot(NULL)+geom_density(aes(x=list_of_draws_null$BMMSY), fill=NA,lty=2,lwd=2)+geom_density(aes(x=list_of_draws_null_bmincovariates$BMMSY), fill="green",alpha=0.5)+theme_classic()+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"))
ggarrange(a,b,labels=c("a","b"))

#model selection
#fit  model with cross validation
full_bmin_kcross <- stan_kfold(file="Full_exports_schaeffer_noncentered_bmin_kfold.stan",stanDat_full_l ,chains=4,cores=2)
loglik_full_bmin_kcross <- extract_log_lik_K(full_bmin_kcross,holdout_data)
kcross_full_bmin <- kfold(loglik_full_bmin_kcross) 
elpd_full_bmin=kcross_full_bmin$elpd_kfold
elpd_se_full_bmin=kcross_full_bmin$se_elpd_kfold
logpreddensity_full_bmin=kcross_full_bmin$pointwise

#compare the models (negative elpd_diff favors  the first model)
#poinwise difference in expected log predictive density: difference in their expected predictive accuracy
elpd_diff_null_full_bmin=sum(elpd_diffs(kcross_null,kcross_full_bmin))
se_elpd_diff_null_full_bmin=se_elpd_diff(elpd_diffs(kcross_null,kcross_full_bmin))
windows()
c=ggplot(NULL, aes(x=elpd_diffs(kcross_null,kcross_full_bmin)))+geom_density(fill="grey", alpha=0.5)+geom_vline(xintercept=0)+xlab("pointwise expected log predictive density difference")+
  geom_text(aes(x=-30,y=0.05,label=paste("sum elpd diff=",round(elpd_diff_null_full_bmin,2),"(SE:",round(se_elpd_diff_null_full_bmin,2),")")))+theme_classic()
ggarrange(c,a,b,labels=c("a","b","c"), nrow=1, ncol=3)


#............................................................................................

#5: Full model vs null model
#status results with reserve only submodel
MMSY_fullmodel=median(list_of_draws_full$MMSY)
BMMSY_fullmodel=median(list_of_draws_full$BMMSY)
B0_fullmodel=median(list_of_draws_full$B0)

#median surplus production curve for  fullmodel
Fit_full_summary=summary(Fit_full)
output_Fit_full=as.data.frame(Fit_full_summary$summary)
popgrowth_full=output_Fit_full[520:870,c("50%")]
popgrowth_full=ifelse(popgrowth_full<0,0.0001,popgrowth_full)
popgrowth_full=as.data.frame(popgrowth_full)
colnames(popgrowth_full)=c("fullmodel_surplus")
popgrowth_full$Biomass=Bio
popgrowth_full2=merge(jurisdictionscale_data,popgrowth_full,by="Biomass",all.x=T)
popgrowth_full2$fullmodel_surplus[!is.na(popgrowth_full2$medianB_tkm2)]
popgrowth_full2$biomassstatus_curve=popgrowth_full2$medianB_tkm2/BMMSY_fullmodel
popgrowth_full2$collapsed=ifelse(is.na(popgrowth_full2$medianB_tkm2),NA,ifelse(popgrowth_full2$medianB_tkm2<(0.1*B0_fullmodel),1,0))
#jurisdictions fishing and fishery status: 
popgrowth_full2$fishingstatus=popgrowth_full2$catch_tkm2/MMSY_fullmodel
popgrowth_full2$overfishing=ifelse(is.na(popgrowth_full2$fishingstatus),NA,ifelse(popgrowth_full2$fishingstatus>1,"Overfishing","Not overfishing"))
popgrowth_full2$fishingstatus_curve=popgrowth_full2$catch_tkm2/popgrowth_full2$fullmodel_surplus
popgrowth_full2$overfishing_curve=ifelse(is.na(popgrowth_full2$fishingstatus_curve),NA,ifelse(popgrowth_full2$fishingstatus_curve>1,1,0))
popgrowth_full2$fisherystatus=as.character(ifelse(is.na(popgrowth_full2$fishingstatus_curve),NA, ifelse(popgrowth_full2$fishingstatus_curve<1 & popgrowth_full2$biomassstatus_curve>1, "Sustainable", ifelse(popgrowth_full2$fishingstatus_curve>1 & popgrowth_full2$biomassstatus_curve<1, "Unsustainable", ifelse(popgrowth_full2$fishingstatus_curve>1 & popgrowth_full2$biomassstatus_curve>1, "Warning", "Rebuilding")))))
#add quasi-sustainable places
popgrowth_full2$fisherystatus=ifelse(is.na(popgrowth_full2$fisherystatus),NA,ifelse(popgrowth_full2$fisherystatus=="Warning" &popgrowth_full2$fishingstatus<1, "Sustainable",as.character(popgrowth_full2$fisherystatus)))

#results
fullmodel_status_results=as.data.frame(c(length(popgrowth_full2$fishingstatus[!is.na(popgrowth_full2$fishingstatus)]),
                                         length(popgrowth_full2$overfishing[!is.na(popgrowth_full2$overfishing)&popgrowth_full2$overfishing=="Overfishing"])/length(popgrowth_full2$overfishing[!is.na(popgrowth_full2$overfishing)]),
                                         length(popgrowth_full2$collapsed[!is.na(popgrowth_full2$collapsed)]),
                                         length(popgrowth_full2$collapsed[!is.na(popgrowth_full2$collapsed)&popgrowth_full2$collapsed==1])/length(popgrowth_full2$collapsed[!is.na(popgrowth_full2$collapsed)]),
                                         length(popgrowth_full2$biomassstatus_curve[!is.na(popgrowth_full2$biomassstatus_curve)&popgrowth_full2$biomassstatus_curve<1])/length(popgrowth_full2$biomassstatus_curve[!is.na(popgrowth_full2$biomassstatus_curve)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)]),
                                         length(popgrowth_full2$overfishing_curve[!is.na(popgrowth_full2$overfishing_curve)&popgrowth_full2$overfishing_curve>0.9])/length(popgrowth_full2$overfishing_curve[!is.na(popgrowth_full2$overfishing_curve)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)&popgrowth_full2$fisherystatus=="Sustainable"])/length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)&popgrowth_full2$fisherystatus=="Unsustainable"])/length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)&popgrowth_full2$fisherystatus=="Warning"])/length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)&popgrowth_full2$fisherystatus=="Rebuilding"])/length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)]),
                                         length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)&!popgrowth_full2$fisherystatus=="Sustainable"])/length(popgrowth_full2$fisherystatus[!is.na(popgrowth_full2$fisherystatus)])))

colnames(fullmodel_status_results)=c("full_model")
row.names(fullmodel_status_results)=c("n_fishingstatus","overfishing_MMSY","n_biomassstatus","collapsed","overfished","n_fisherystatus","overfishing_curve", "sustainable","unsustainable","warning","rebuilding","conservation_concern")
#write.csv(fullmodel_status_results,"sensitivity_status_fullmodelonly_quasi.csv")

#............................................................................................
#4: Sesitivity analysis: Choice of surplus production model

#Now we try the Pella-tomlinson surplus model
Fit_null_pt <- stan(file = "Null_pellatomlinson_free.stan", data = stanDat_null, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_pt, pars = c("r","bmin","B0","n","lp__"))
pairs(Fit_null_pt, pars = c("MMSY","BMMSY", "lp__"))
#this model does not converge: you can get virtually identical fits by covarying two of the model parameters in a particular way (i.e., a ridge in the likelihood function). Pointed out by Fletcher in the 1970s

#Consequently we try the Pella-tomlinson putting a contraint in n from 0.5 to 4. Prior: uniform (0.5,4)
Fit_null_pt2 <- stan(file = "Null_pellatomlinson_constraints_free.stan", data = stanDat_null, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_pt2, pars = c("r","bmin","B0","n","lp__"))
pairs(Fit_null_pt2, pars = c("MMSY","BMMSY", "lp__"))
# converges (but gives any potentialvalye of n

#Thus, we fit the Pella-tomlinson fixed with 3 and 4  for n (Quinn and Deriso)
stanDat_null_fixedn3 <- list(n=3,b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                             N = nrow(reserves_complete),
                             br=log(remote_complete$FamBiomass_tkm2),
                             RE=nrow(remote_complete),
                             newage=newdata$Closure.age,
                             Bio=Bio, B=B)

Fit_null_pt3 <- stan(file = "Null_pellatomlinson_fixedn.stan", data = stanDat_null_fixedn3, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_pt3, pars = c("r","bmin","B0","lp__"))
pairs(Fit_null_pt3, pars = c("MMSY","BMMSY", "lp__"))
stanDat_null_fixedn4 <- list(n=4,b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                             N = nrow(reserves_complete),
                             br=log(remote_complete$FamBiomass_tkm2),
                             RE=nrow(remote_complete),
                             newage=newdata$Closure.age,
                             Bio=Bio, B=B)
Fit_null_pt4 <- stan(file = "Null_pellatomlinson_fixedn.stan", data = stanDat_null_fixedn4, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_pt4, pars = c("r","bmin","B0","lp__"))
pairs(Fit_null_pt4, pars = c("MMSY","BMMSY", "lp__"))
#when we fix the n parameter the model does converge


#fox model (Gompertz growth function instead of logistic or equivalently, the Pella-Tomlinson with n=1)
Fit_null_fox <- stan(file = "Null_Gompertzfox.stan", data = stanDat_null, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
pairs(Fit_null_fox, pars = c("r","bmin","B0","lp__"))
pairs(Fit_null_fox, pars = c("MMSY","BMMSY", "lp__"))

#extract posteriors
list_of_draws_null_pt3 <- as.data.frame(Fit_null_pt3) 
list_of_draws_null_pt4 <- as.data.frame(Fit_null_pt4) 
list_of_draws_null_fox <- as.data.frame(Fit_null_fox) 

#difference in parameter distributions
MMSYdata=as.data.frame(cbind(c(list_of_draws_null$MMSY,list_of_draws_null_fox$MMSY,list_of_draws_null_pt3$MMSY,list_of_draws_null_pt4$MMSY), c(list_of_draws_null$BMMSY,list_of_draws_null_fox$BMMSY,list_of_draws_null_pt3$BMMSY,list_of_draws_null_pt4$BMMSY),c(list_of_draws_null$uMMSY,list_of_draws_null_fox$uMMSY,list_of_draws_null_pt3$uMMSY,list_of_draws_null_pt4$uMMSY),c(rep("Schaefer", length.out=length(list_of_draws_null$MMSY)),rep("Gompertz-Fox", length.out=length(list_of_draws_null_fox$MMSY)),rep("Pella-Tomlinson 3", length.out=length(list_of_draws_null_pt3$MMSY)),rep("Pella-Tomlinson 4", length.out=length(list_of_draws_null_pt4$MMSY)))))
colnames(MMSYdata)=c("MMSY","BMMSY","uMMSY","model")
MMSYdata$MMSY=as.numeric(as.character(MMSYdata$MMSY))
MMSYdata$BMMSY=as.numeric(as.character(MMSYdata$BMMSY))
MMSYdata$uMMSY=as.numeric(as.character(MMSYdata$uMMSY))

a=ggplot(MMSYdata, aes(x=MMSY,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+xlab("MMSY (t/km2/y)")+ylab("Posterior density")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))

b=ggplot(MMSYdata, aes(x=BMMSY,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+
  theme_classic()+xlab("BMMSY (t/km2)")+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"),fill="Model")
c=ggplot(MMSYdata, aes(x=uMMSY,  group=model)) +
  geom_density( aes(fill=MMSYdata$model), alpha=0.5)+
  labs(fill="Model")+theme_classic()+xlab("uMMSY (1/y)")+ylab("")
windows()
ggarrange(a,b, widths=c(1,1.5), nrow=1, ncol=2)


#fishing status and biomass status under different surplus models
#median MMSY fro each
MMSY_surplus=c( median(list_of_draws_null_fox$MMSY),median(list_of_draws_null$MMSY),median(list_of_draws_null_pt3$MMSY), median(list_of_draws_null_pt4$MMSY))
BMMSY_surplus=c( median(list_of_draws_null_fox$BMMSY),median(list_of_draws_null$BMMSY),median(list_of_draws_null_pt3$BMMSY), median(list_of_draws_null_pt4$BMMSY))

#surplus production curve for  each one
fishingstatus_surplus=matrix(NA, nrow=length(jurisdictionscale_data$catch_tkm2),ncol=length(MMSY_surplus))
biomassstatus_surplus=matrix(NA, nrow=length(jurisdictionscale_data$medianB_tkm2),ncol=length(BMMSY_surplus))
for (i in 1:length(MMSY_surplus)){
  fishingstatus_surplus[,i]=jurisdictionscale_data$catch_tkm2/MMSY_surplus[i]
  biomassstatus_surplus[,i]=jurisdictionscale_data$medianB_tkm2/BMMSY_surplus[i]
}

overfishing_surplus=ifelse(is.na(fishingstatus_surplus),NA, ifelse(fishingstatus_surplus>1,"Overfishing","Not overfishing"))
overfished_surplus=ifelse(is.na(biomassstatus_surplus),NA, ifelse(biomassstatus_surplus<1,"Overfished","Not overfished"))
surplus_status_sensitivity=matrix(NA,nrow=2,ncol=length(MMSY_surplus))
for (i in 1:length(MMSY_surplus)){
  surplus_status_sensitivity[,i]=c(length(overfishing_surplus[,i][!is.na(overfishing_surplus[,i])&overfishing_surplus[,i]=="Overfishing"])/length(overfishing_surplus[,i][!is.na(overfishing_surplus[,i])]),
                                   length(overfished_surplus[,i][!is.na(overfished_surplus[,i])&overfished_surplus[,i]=="Overfished"])/length(overfished_surplus[,i][!is.na(overfished_surplus[,i])]))
}
row.names(surplus_status_sensitivity)=c("Overfishing","Overfished")
colnames(surplus_status_sensitivity)=c("Gompertz-Fox", "Graham-Schaeffer","Pella-Tomlinson n=3","Pella-Tomlinson n=4")
#write.csv(surplus_status_sensitivity,"sensitivity_status_surplus.csv")
#............................................................................................

#5: Sesitivity analysis: CLOSED VS OPEN POPULATIONS IN RESERVES 
#different assumptions on exports
#(a) exports as a function of reserve age (i.e., smaller reserves export more)
stanDat_null_exprts_size <- list(b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,si=reserves_complete$sClosure.size,
                                 N = nrow(reserves_complete),
                                 br=log(remote_complete$FamBiomass_tkm2),
                                 RE=nrow(remote_complete),
                                 newage=newdata$Closure.age,
                                 Bio=Bio, B=B)

Fit_null_exports_size <- stan(file = "Null_exportsfunctionreservesize.stan", data = stanDat_null_exprts_size , iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
windows()
pairs(Fit_null_exports_size  , pars = c("r","bmin","B0","sigma_e","slope","lp__"))
list_of_draws_null_exports_size<- as.data.frame(Fit_null_exports_size) 

#(b) exports as a function of the difference between observed b and the estimated bmin (i.e., greater difference more exports due to density dependence)
Fit_null_exports_bbmin <- stan(file = "Null_exportsfunctionbminusbmin.stan", data = stanDat_null, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
windows()
pairs(Fit_null_exports_bbmin , pars = c("r","bmin","B0","sigma_e","slope","lp__"))
list_of_draws_null_exports_bbmin<- as.data.frame(Fit_null_exports_bbmin) 

comparisondata2=as.data.frame(cbind(c(list_of_draws_null$MMSY,list_of_draws_null_exports_bbmin$MMSY,list_of_draws_null_exports_size$MMSY),c(list_of_draws_null$BMMSY,list_of_draws_null_exports_bbmin$BMMSY,list_of_draws_null_exports_size$BMMSY),c(list_of_draws_null$r,list_of_draws_null_exports_bbmin$r,list_of_draws_null_exports_size$r), c(list_of_draws_null$B0,list_of_draws_null_exports_bbmin$B0,list_of_draws_null_exports_size$B0),c(list_of_draws_null$bmin,list_of_draws_null_exports_bbmin$bmin,list_of_draws_null_exports_size$bmin),c(rep("zero (original)", length.out=length(list_of_draws_null$r)),rep("f(B-bmin)", length.out=length(list_of_draws_null_exports_bbmin$r)),rep("f(Reserve size)", length.out=length(list_of_draws_null_exports_size$r)))))
colnames(comparisondata2)=c("MMSY","BMMSY","r","B0","bmin","exports")
comparisondata2$r=as.numeric(as.character(comparisondata2$r))
comparisondata2$B0=as.numeric(as.character(comparisondata2$B0))
comparisondata2$bmin=as.numeric(as.character(comparisondata2$bmin))
comparisondata2$MMSY=as.numeric(as.character(comparisondata2$MMSY))
comparisondata2$BMMSY=as.numeric(as.character(comparisondata2$BMMSY))

a=ggplot(comparisondata2, aes(x=bmin,  group=exports)) +
  geom_density( aes(fill=comparisondata2$exports), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("Posterior density")+ labs(x =  expression ("Biomass reserve age 0 ("~t/km^2*")"))
b=ggplot(comparisondata2, aes(x=B0,  group=exports)) +
  geom_density( aes(fill=comparisondata2$exports), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+xlab("B0 (t/km2)")+ylab("")+ labs(x =  expression ("Unfished Biomass ("~t/km^2*")"))
c=ggplot(comparisondata2, aes(x=r,  group=exports)) +guides(fill=FALSE)+
  geom_density( aes(fill=comparisondata2$exports), alpha=0.5)+
  labs(fill="exports")+theme_classic()+xlab("intrinsic growth rate")+ylab("")
d=ggplot(comparisondata2, aes(x=MMSY,  group=exports)) +
  geom_density( aes(fill=comparisondata2$exports), alpha=0.5)+guides(fill=FALSE)+
  theme_classic()+ylab("")+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))+xlim(c(0,10))
e=ggplot(comparisondata2, aes(x=BMMSY,  group=exports)) +
  geom_density( aes(fill=comparisondata2$exports), alpha=0.5)  +labs(fill="exports")+
  theme_classic()+ylab("")+ labs(x =  expression ("BMMSY ("~t/km^2*")"))
ggarrange(a,b,c,d,e, widths=c(1,1,1,1,1.5), nrow=1, ncol=5)


#(c): exports as a fixed proportion
exportvector=c(0,0.05,0.10,0.15)
unfishedbiomass=matrix(NA, nrow=4000, ncol=length(exportvector))
littler=matrix(NA, nrow=4000, ncol=length(exportvector))
bioage0=matrix(NA, nrow=4000, ncol=length(exportvector))
MMSY_exports=matrix(NA, nrow=4000, ncol=length(exportvector))
BMMSY_exports=matrix(NA, nrow=4000, ncol=length(exportvector))
median_model_estimates=matrix(NA, nrow=length(output_Fit_null$`50%`), ncol=length(exportvector))
for (i in 1: length(exportvector)){
  stanDat_null_exports <- list(p=exportvector[i], b = log(reserves_complete$FamBiomass_tkm2),ag=reserves_complete$Closure.age,
                               N = nrow(reserves_complete),
                               br=log(remote_complete$FamBiomass_tkm2),
                               RE=nrow(remote_complete),
                               newage=newdata$Closure.age,
                               Bio=Bio, B=B) 
  Fit_null_exports <- stan(file = "Null_exports_schaeffer.stan", data = stanDat_null_exports, iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.999,stepsize = 0.005,max_treedepth = 20))
  list_of_draws_null_exports <- as.data.frame(Fit_null_exports) 
  Fit_null_exports_summary <- summary(Fit_null_exports,probs=c(0.05,0.25,0.5,0.75, 0.95))
  output_Fit_null_exports=as.data.frame(Fit_null_exports_summary$summary)
  median_model_estimates[,i]=output_Fit_null_exports$`50%`
  unfishedbiomass[,i]=list_of_draws_null_exports$B0
  littler[,i]=list_of_draws_null_exports$r
  bioage0[,i]=list_of_draws_null_exports$bmin
  MMSY_exports[,i]=list_of_draws_null_exports$MMSY
  BMMSY_exports[,i]=list_of_draws_null_exports$BMMSY
}
unfishedbiomass=as.data.frame(unfishedbiomass)
colnames(unfishedbiomass)=c("0%","5%","10%","15%")
unfishedbiomass1=melt(unfishedbiomass)
a=ggplot(unfishedbiomass1, aes (value)) + guides(fill=F)+
  geom_density(aes(fill = variable), alpha=0.3)+theme_classic()+ labs(x =  expression ("Unfished Biomass ("~t/km^2*")"))+ylab("Posterior density")

littler=as.data.frame(littler)
colnames(littler)=c("0%","5%","10%","15%")
littler1=melt(littler)
b=ggplot(littler1, aes (value)) +guides(fill=F)+
  geom_density(aes(fill = variable), alpha=0.3)+theme_classic()+xlab("intrinsic growth rate")+ylab("")

bioage0=as.data.frame(bioage0)
colnames(bioage0)=c("0%","5%","10%","15%")
bioage01=melt(bioage0)
c=ggplot(bioage01, aes (value)) +guides(fill=F)+
  geom_density(aes(fill = variable), alpha=0.3)+theme_classic()+ labs(x =  expression ("Biomass reserve age 0 ("~t/km^2*")"))+ylab("")

MMSY_exports=as.data.frame(MMSY_exports)
colnames(MMSY_exports)=c("0%","5%","10%","15%")
MMSY_exports1=melt(MMSY_exports)
d=ggplot(MMSY_exports1, aes (value)) +guides(fill=F)+
  geom_density(aes(fill = variable), alpha=0.3)+theme_classic()+ labs(x =  expression ("MMSY ("~t/km^2/y*")"))+ylab("")

BMMSY_exports=as.data.frame(BMMSY_exports)
colnames(BMMSY_exports)=c("0%","5%","10%","15%")
BMMSY_exports1=melt(BMMSY_exports)
e=ggplot(BMMSY_exports1, aes (value)) +
  geom_density(aes(fill = variable), alpha=0.3)+theme_classic()+labs(fill="Fixed exports")+ labs(x =  expression ("BMMSY ("~t/km^2*")"))+ylab("")
windows()
ggarrange(a,b,c,d,e,nrow=1,ncol=5, widths=c(1,1,1,1,1.5))

#median MMSY fro each
MMSY_fixedexports=matrixStats::colMedians(as.matrix.data.frame(MMSY_exports))

#surplus production curve for  each one
jurisdictionscale_data2$Biomass=jurisdictionscale_data2$roundedmedianB
popgrowth_null_exports=median_model_estimates[309:659,]
popgrowth_null_exports=ifelse(popgrowth_null_exports<0,0.0001,popgrowth_null_exports)
popgrowth_null_exports=as.data.frame(popgrowth_null_exports)
colnames(popgrowth_null_exports)=c("0%e","5%e","10%e","15%e")
popgrowth_null_exports$Biomass=Bio
popgrowth_null_exports2=merge(jurisdictionscale_data2,popgrowth_null_exports,by="Biomass",all.x=T)
popgrowth_null_exports2$biomassstatus_curve=popgrowth_null_exports2$medianB_tkm2/BMMSY_median 
popgrowth_null_exports3=popgrowth_null_exports2[,c("0%e","5%e","10%e","15%e")]

#jurisdictions fishing and fishery status: 
export_sensitivity_fishingstatus=matrix(NA, nrow=length(popgrowth_null_exports2$catch_tkm2),ncol=length(MMSY_fixedexports))
export_sensitivity_overfishing=matrix(NA, nrow=length(popgrowth_null_exports2$catch_tkm2),ncol=length(MMSY_fixedexports))
export_sensitivity_fishingstatus_curve=matrix(NA, nrow=length(popgrowth_null_exports2$catch_tkm2),ncol=length(MMSY_fixedexports))
export_sensitivity_overfishing_curve=matrix(NA, nrow=length(popgrowth_null_exports2$catch_tkm2),ncol=length(MMSY_fixedexports))
export_sensitivity_fisherystatus=matrix(NA, nrow=length(popgrowth_null_exports2$catch_tkm2),ncol=length(MMSY_fixedexports))
export_status_results=matrix(NA,nrow=nrow(status_results)+1,ncol=length(MMSY_fixedexports))

for (i in 1:length(MMSY_fixedexports)){
  export_sensitivity_fishingstatus[,i]=popgrowth_null_exports2$catch_tkm2/MMSY_fixedexports[i]
  export_sensitivity_overfishing[,i]=ifelse(is.na(export_sensitivity_fishingstatus[,i]),NA,ifelse(export_sensitivity_fishingstatus[,i]>1,"Overfishing","Not overfishing"))
  print(length(export_sensitivity_overfishing[,i][!is.na(export_sensitivity_overfishing[,i])&export_sensitivity_overfishing[,i]=="Overfishing"])/length(export_sensitivity_overfishing[,i][!is.na(export_sensitivity_overfishing[,i])]))
  export_sensitivity_fishingstatus_curve[,i]=popgrowth_null_exports2$catch_tkm2/popgrowth_null_exports3[,i]
  export_sensitivity_overfishing_curve[,i]=ifelse(is.na(export_sensitivity_fishingstatus_curve[,i]),NA,ifelse(export_sensitivity_fishingstatus_curve[,i]>1,1,0))
  export_sensitivity_fisherystatus[,i]=as.character(ifelse(is.na(export_sensitivity_fishingstatus_curve[,i]),NA, ifelse(is.na(popgrowth_null_exports2$biomassstatus_curve),NA,ifelse(export_sensitivity_fishingstatus_curve[,i]<1 & popgrowth_null_exports2$biomassstatus_curve>1, "Sustainable", ifelse(export_sensitivity_fishingstatus_curve[,i]>1 & popgrowth_null_exports2$biomassstatus_curve<1, "Unsustainable", ifelse(export_sensitivity_fishingstatus_curve[,i]>1 & popgrowth_null_exports2$biomassstatus_curve>1, "Warning", "Rebuilding"))))))
  #add quasi-sustainable places
  export_sensitivity_fisherystatus[,i]=as.character(ifelse(is.na(export_sensitivity_fishingstatus_curve[,i]),NA, ifelse(is.na(popgrowth_null_exports2$biomassstatus_curve),NA,ifelse(export_sensitivity_fishingstatus[,i]<1 & popgrowth_null_exports2$biomassstatus_curve>1, "Sustainable", ifelse(export_sensitivity_fishingstatus_curve[,i]>1 & popgrowth_null_exports2$biomassstatus_curve<1, "Unsustainable", ifelse(export_sensitivity_fishingstatus[,i]>1 & popgrowth_null_exports2$biomassstatus_curve>1, "Warning", "Rebuilding"))))))
  export_status_results[,i]=c(length(jurisdictionscale_data$overfishing[!is.na(jurisdictionscale_data$catch_tkm2)]),
                              length(export_sensitivity_overfishing[,i][!is.na(export_sensitivity_overfishing[,i])&export_sensitivity_overfishing[,i]=="Overfishing"])/length(export_sensitivity_overfishing[,i][!is.na(export_sensitivity_overfishing[,i])]),
                              length(Bbycountry$collapsed),
                              length(Bbycountry$collapsed[Bbycountry$collapsed==1])/length(Bbycountry$collapsed),
                              length(Bbycountry$overfished[Bbycountry$overfished=="Overfished"])/length(Bbycountry$overfished),
                              length(export_sensitivity_fishingstatus_curve[,i][!is.na(export_sensitivity_fishingstatus_curve[,i])]),
                              length(export_sensitivity_overfishing_curve[,i][!is.na(export_sensitivity_overfishing_curve[,i])&export_sensitivity_overfishing_curve[,i]>0.9])/length(export_sensitivity_overfishing_curve[,i][!is.na(export_sensitivity_overfishing_curve[,i])]),
                              length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])&export_sensitivity_fisherystatus[,i]=="Sustainable"])/length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])]),
                              length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])&export_sensitivity_fisherystatus[,i]=="Unsustainable"])/length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])]),
                              length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])&export_sensitivity_fisherystatus[,i]=="Warning"])/length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])]),
                              length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])&export_sensitivity_fisherystatus[,i]=="Rebuilding"])/length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])]),
                              length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])&!export_sensitivity_fisherystatus[,i]=="Sustainable"])/length(export_sensitivity_fisherystatus[,i][!is.na(export_sensitivity_fisherystatus[,i])]))
}

export_sensitivity_fisherystatus=as.data.frame(export_sensitivity_fisherystatus)
export_status_results=as.data.frame( export_status_results)
colnames( export_status_results)=c("0%","5%","10%","15%")
row.names(export_status_results)=c("n_fishingstatus","overfishing_MMSY","n_biomassstatus","collapsed","overfished","n_fisherystatus","overfishing_curve", "sustainable","unsustainable","warning","rebuilding","conservation_concern")
#write.csv(export_status_results,"sensitivity_status_exports_quasi.csv")


##..............................................................................................
#6: Choice of catch data
catchdataoptions=jurisdictionscale_data[,c("mean_spatial_totalcatch_reefs_t","mean_totalcatch_t","mean_tcatch_notindus_t","mean_tcatch_onlyreported_t", "mean_tcatch_notindus_onlyrep_t")]
#tnnes per km2 reef
catch_km2_options=catchdataoptions/jurisdictionscale_data$CoralReefArea_km2

#fishing status
fishingstatus_catch=catch_km2_options/MMSY_median
pairs(~., data=log(fishingstatus_catch))
overfishing_catch=ifelse(is.na(fishingstatus_catch),NA,ifelse(fishingstatus_catch>1,"Overfishing","Not overfishing"))

#fishing sttaus based on surplus curve
catch_km2_options2=catch_km2_options
catch_km2_options2$area_name=jurisdictionscale_data$area_name
catch_km2_options2=merge(catch_km2_options2,jurisdictionscale_data2[,c("area_name","Biomass","medianB_tkm2","stotalgravity_imp","smeanttm_imp","sfisherdens_imp","stouristdens_imp","sHDI_imp","spopgrowth_imp","smpa_imp")],by="area_name",all.x=T)
catch_km2_options2=merge(catch_km2_options2,popgrowth_null[,c("Biomass","median")],by="Biomass",all.x=T)
catch_km2_options2$avyield=ifelse(is.na(catch_km2_options2$median),NA, catch_km2_options2$median)
catch_km2_options2$median=NULL
catch_km2_options2$biomassstatus=catch_km2_options2$medianB_tkm2/BMMSY_median
catch_km2_options3=catch_km2_options2[,3:7]
fishingstatus_curve_catch=catch_km2_options3/catch_km2_options2$avyield
overfishing_curve_catch=ifelse(is.na(fishingstatus_curve_catch),NA,ifelse(fishingstatus_curve_catch>1,1,0))
fishingstatus_curve_catch2=fishingstatus_curve_catch
fishingstatus_curve_catch2$biomassstatus=catch_km2_options2$biomassstatus
fisherystatus_catch=ifelse(is.na(fishingstatus_curve_catch),NA, ifelse(fishingstatus_curve_catch<1 & fishingstatus_curve_catch2$biomassstatus>1, "Sustainable", ifelse(fishingstatus_curve_catch>1 & fishingstatus_curve_catch2$biomassstatus<1, "Unsustainable", ifelse(fishingstatus_curve_catch>1 & fishingstatus_curve_catch2$biomassstatus>1, "Warning", "Rebuilding"))))
#add quasi
fishingstatus_catch2=catch_km2_options3/MMSY_median
fisherystatus_catch=ifelse(is.na(fisherystatus_catch),NA,ifelse(fisherystatus_catch=="Warning" &fishingstatus_catch2<1, "Sustainable",as.character(fisherystatus_catch)))

#sumamry results
catch_status_sensitivity=matrix(NA,nrow=nrow(status_results)+1,ncol=ncol(overfishing_catch))

for (i in 1:ncol(overfishing_catch)){
  catch_status_sensitivity[,i]=c(length(fishingstatus_catch[,i][!is.na(fishingstatus_catch[,i])]),
                                 length(overfishing_catch[,i][!is.na(overfishing_catch[,i])&overfishing_catch[,i]=="Overfishing"])/length(overfishing_catch[,i][!is.na(overfishing_catch[,i])]),
                                 length(Bbycountry$collapsed),
                                 length(Bbycountry$collapsed[Bbycountry$collapsed==1])/length(Bbycountry$collapsed),
                                 length(Bbycountry$overfished[Bbycountry$overfished=="Overfished"])/length(Bbycountry$overfished),
                                 length(fishingstatus_curve_catch[,i][!is.na(fishingstatus_curve_catch[,i])]),
                                 length(overfishing_curve_catch[,i][!is.na(overfishing_curve_catch[,i])&overfishing_curve_catch[,i]>0.9])/length(overfishing_curve_catch[,i][!is.na(overfishing_curve_catch[,i])]),
                                 length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Sustainable"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]),
                                 length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Unsustainable"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]),
                                 length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Warning"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]),
                                 length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&fisherystatus_catch[,i]=="Rebuilding"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]),
                                 length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])&!fisherystatus_catch[,i]=="Sustainable"])/length(fisherystatus_catch[,i][!is.na(fisherystatus_catch[,i])]))
}

catch_status_sensitivity=as.data.frame( catch_status_sensitivity)
colnames(catch_km2_options)
colnames( catch_status_sensitivity)=c("spatial","total reef catch","non-industrial reef catch","reported reef catch","reported non-industrial reef catch")
row.names(catch_status_sensitivity)=c("n_fishingstatus","overfishing_MMSY","n_biomassstatus","collapsed","overfished","n_fisherystatus","overfishing_curve", "sustainable","unsustainable","warning","rebuilding","conservation_concern")
#write.csv(catch_status_sensitivity,"sensitivity_status_catchdata_quasi.csv")

#covariate model robustness (log(fishing status))
lfishingstatus_catch=log(fishingstatus_catch+0.00000000001)
lfishingstatus_catch$area_name=jurisdictionscale_data$area_name
lfishingstatus_catch=merge(lfishingstatus_catch,jurisdictionscale_data4[,c("area_name","stotalgravity_imp","smeanttm_imp","sfisherdens_imp","stouristdens_imp","sHDI_imp","spopgrowth_imp","smpa_imp")],by="area_name",all.x=T)


fstatus_model_original=lm(mean_spatial_totalcatch_reefs_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=lfishingstatus_catch)
fstatus_model_totalcatch=lm(mean_totalcatch_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=lfishingstatus_catch)
fstatus_model_nonindustrial=lm(mean_tcatch_notindus_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=lfishingstatus_catch)
fstatus_model_reported=lm(mean_tcatch_onlyreported_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=lfishingstatus_catch)
fstatus_model_reportednonindustrial=lm(mean_tcatch_notindus_onlyrep_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=lfishingstatus_catch)

#coefplots
coefplot_original=as.data.frame(confint(fstatus_model_original, level=0.9))
coefplot_original$mean=fstatus_model_original$coefficients
coefplot_original=coefplot_original[-1,]
coefplot_original$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_original$sign=ifelse(coefplot_original$`5 %`<0 & coefplot_original$`95 %`<0, "negative",ifelse(coefplot_original$`5 %`>0 & coefplot_original$`95 %`>0, "positive", "no effect"))

coefplotf_original=
  ggplot(coefplot_original,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Spatial")
coefplot_totalcatch=as.data.frame(confint(fstatus_model_totalcatch, level=0.9))
coefplot_totalcatch$mean=fstatus_model_totalcatch$coefficients
coefplot_totalcatch=coefplot_totalcatch[-1,]
coefplot_totalcatch$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_totalcatch$sign=ifelse(coefplot_totalcatch$`5 %`<0 & coefplot_totalcatch$`95 %`<0, "negative",ifelse(coefplot_totalcatch$`5 %`>0 & coefplot_totalcatch$`95 %`>0, "positive", "no effect"))

coefplotf_totalcatch=
  ggplot(coefplot_totalcatch,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Total")

coefplot_nonindustrial=as.data.frame(confint(fstatus_model_nonindustrial, level=0.9))
coefplot_nonindustrial$mean=fstatus_model_nonindustrial$coefficients
coefplot_nonindustrial=coefplot_nonindustrial[-1,]
coefplot_nonindustrial$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_nonindustrial$sign=ifelse(coefplot_nonindustrial$`5 %`<0 & coefplot_nonindustrial$`95 %`<0, "negative",ifelse(coefplot_nonindustrial$`5 %`>0 & coefplot_nonindustrial$`95 %`>0, "positive", "no effect"))

coefplotf_nonindustrial=
  ggplot(coefplot_nonindustrial,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Non-industrial")

coefplot_reported=as.data.frame(confint(fstatus_model_reported, level=0.9))
coefplot_reported$mean=fstatus_model_reported$coefficients
coefplot_reported=coefplot_reported[-1,]
coefplot_reported$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_reported$sign=ifelse(coefplot_reported$`5 %`<0 & coefplot_reported$`95 %`<0, "negative",ifelse(coefplot_reported$`5 %`>0 & coefplot_reported$`95 %`>0, "positive", "no effect"))

coefplotf_reported=
  ggplot(coefplot_reported,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Reported")

coefplot_reportednonindustrial=as.data.frame(confint(fstatus_model_reportednonindustrial, level=0.9))
coefplot_reportednonindustrial$mean=fstatus_model_reportednonindustrial$coefficients
coefplot_reportednonindustrial=coefplot_reportednonindustrial[-1,]
coefplot_reportednonindustrial$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_reportednonindustrial$sign=ifelse(coefplot_reportednonindustrial$`5 %`<0 & coefplot_reportednonindustrial$`95 %`<0, "negative",ifelse(coefplot_reportednonindustrial$`5 %`>0 & coefplot_reportednonindustrial$`95 %`>0, "positive", "no effect"))

coefplotf_reportednonindustrial=
  ggplot(coefplot_reportednonindustrial,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Reported non-industrial")
windows()
catch_sensitivity_coefplots=ggarrange(coefplotf_original,coefplotf_totalcatch,coefplotf_nonindustrial,coefplotf_reported,coefplotf_reportednonindustrial, nrow=1,ncol=5, widths=c(1,1,1,1,1.7), labels=c("a","b","c","d","e"))
annotate_figure(catch_sensitivity_coefplots,bottom = "Stnd. effect size")

#alternative using binomial models

overfishing_catch2=ifelse(is.na(fishingstatus_catch),NA,ifelse(fishingstatus_catch>1,1,0))
overfishing_catch2=as.data.frame(overfishing_catch2)
overfishing_catch2$area_name=jurisdictionscale_data$area_name
overfishing_catch2=merge(overfishing_catch2,jurisdictionscale_data4[,c("area_name","stotalgravity_imp","smeanttm_imp","sfisherdens_imp","stouristdens_imp","sHDI_imp","spopgrowth_imp","smpa_imp")],by="area_name",all.x=T)

overfishing_model_original=glm(mean_spatial_totalcatch_reefs_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=overfishing_catch2, family=binomial(link="logit"))
overfishing_model_totalcatch=glm(mean_totalcatch_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=overfishing_catch2, family=binomial(link="logit"))
overfishing_model_nonindustrial=glm(mean_tcatch_notindus_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=overfishing_catch2, family=binomial(link="logit"))
overfishing_model_reported=glm(mean_tcatch_onlyreported_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=overfishing_catch2, family=binomial(link="logit"))
overfishing_model_reportednonindustrial=glm(mean_tcatch_notindus_onlyrep_t~stotalgravity_imp+sfisherdens_imp+stouristdens_imp+sHDI_imp+I(sHDI_imp^2)+spopgrowth_imp+smpa_imp+smeanttm_imp, data=overfishing_catch2, family=binomial(link="logit"))

#coefplots
coefplot_original=as.data.frame(confint(overfishing_model_original, level=0.9))
coefplot_original$mean=overfishing_model_original$coefficients
coefplot_original=coefplot_original[-1,]
coefplot_original$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_original$sign=ifelse(coefplot_original$`5 %`<0 & coefplot_original$`95 %`<0, "negative",ifelse(coefplot_original$`5 %`>0 & coefplot_original$`95 %`>0, "positive", "no effect"))

coefplotf_original=
  ggplot(coefplot_original,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Spatial")
coefplot_totalcatch=as.data.frame(confint(overfishing_model_totalcatch, level=0.9))
coefplot_totalcatch$mean=overfishing_model_totalcatch$coefficients
coefplot_totalcatch=coefplot_totalcatch[-1,]
coefplot_totalcatch$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_totalcatch$sign=ifelse(coefplot_totalcatch$`5 %`<0 & coefplot_totalcatch$`95 %`<0, "negative",ifelse(coefplot_totalcatch$`5 %`>0 & coefplot_totalcatch$`95 %`>0, "positive", "no effect"))

coefplotf_totalcatch=
  ggplot(coefplot_totalcatch,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Total")

coefplot_nonindustrial=as.data.frame(confint(overfishing_model_nonindustrial, level=0.9))
coefplot_nonindustrial$mean=overfishing_model_nonindustrial$coefficients
coefplot_nonindustrial=coefplot_nonindustrial[-1,]
coefplot_nonindustrial$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_nonindustrial$sign=ifelse(coefplot_nonindustrial$`5 %`<0 & coefplot_nonindustrial$`95 %`<0, "negative",ifelse(coefplot_nonindustrial$`5 %`>0 & coefplot_nonindustrial$`95 %`>0, "positive", "no effect"))

coefplotf_nonindustrial=
  ggplot(coefplot_nonindustrial,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Non-industrial")

coefplot_reported=as.data.frame(confint(overfishing_model_reported, level=0.9))
coefplot_reported$mean=overfishing_model_reported$coefficients
coefplot_reported=coefplot_reported[-1,]
coefplot_reported$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_reported$sign=ifelse(coefplot_reported$`5 %`<0 & coefplot_reported$`95 %`<0, "negative",ifelse(coefplot_reported$`5 %`>0 & coefplot_reported$`95 %`>0, "positive", "no effect"))

coefplotf_reported=
  ggplot(coefplot_reported,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),axis.text.y =element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Reported")

coefplot_reportednonindustrial=as.data.frame(confint(overfishing_model_reportednonindustrial, level=0.9))
coefplot_reportednonindustrial$mean=overfishing_model_reportednonindustrial$coefficients
coefplot_reportednonindustrial=coefplot_reportednonindustrial[-1,]
coefplot_reportednonindustrial$variable=c("Total Gravity", "Fisher density", "Tourism density","HDI","HDI^2", "Population growth", "% protected territorial waters","Travel time nearest markets")
coefplot_reportednonindustrial$sign=ifelse(coefplot_reportednonindustrial$`5 %`<0 & coefplot_reportednonindustrial$`95 %`<0, "negative",ifelse(coefplot_reportednonindustrial$`5 %`>0 & coefplot_reportednonindustrial$`95 %`>0, "positive", "no effect"))

coefplotf_reportednonindustrial=
  ggplot(coefplot_reportednonindustrial,aes(x=variable,y=mean,ymin=`5 %`,ymax=`95 %`))+
  geom_pointrange(aes(colour=sign),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "top") +scale_y_continuous(position = "left")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0),panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()+ggtitle("   Reported non-industrial")
windows()
catch_sensitivity_coefplots=ggarrange(coefplotf_original,coefplotf_totalcatch,coefplotf_nonindustrial,coefplotf_reported,coefplotf_reportednonindustrial, nrow=1,ncol=5, widths=c(1,1,1,1,1.7), labels=c("a","b","c","d","e"))
annotate_figure(catch_sensitivity_coefplots,bottom = "log odds ratio")


#############################################################
#save.image(file='Zamborain-Masonetal2020_Revised.RData')
#load('Zamborain-Masonetal2020_Revised.RData')


