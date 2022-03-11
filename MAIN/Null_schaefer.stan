data {
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
 real<lower=0,upper=1> p;// proportion of biomass exported 
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
generated quantities {
 real <lower=0> BMMSY; //BMMSY for averge and most common environemental conditions
 real <lower=0> MMSY;//MMSY for average and most common environmental conditions
 vector[res+rem+fis] log_lik; //log-likelihood (for loo)
 for (n in 1:res) {
 log_lik[n] = normal_lpdf(b[n]| mu[n], sigma_e);
 }
for (n in (res+1):(res+rem)) {
 log_lik[n] = normal_lpdf(b2[n-res]| mu2[n-res], sigma_r);
}
for (n in (res+rem+1):(res+rem+fis)) {
 log_lik[n] = normal_lpdf(b3[n-(res+rem)]| mu3[n-(res+rem)], sigma_f);
}
//average and common conditions
 BMMSY=B0/2;
 MMSY=((r*B0)/4); 
}



