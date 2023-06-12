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
 int<lower=1> R3; //number of data regions fished(groups)
 int<lower=1, upper=R3> pr3[fis]; //region id 
  real p; //fixed export
 }
 parameters {
   vector[12] beta; //slopes for sampling/environmental variables 
   real<lower=0> sigma_e; //error sd for biomass reserves
   real<lower=0> sigma_r; //error sd for biomass remote
   real<lower=0> sigma_f; //error sd for biomass fished
   vector[R3] u_tilde3; // random effects for fished regions
   real<lower=0> sigma_u; //region sd
   real log_r;//intrinsic growth rate
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
   vector[R3] u3; // random effects for fished regions
   real<lower=0> r;
   real<lower=0> bmin; //biomass at reserve age 0
   real<lower=0> B0; //unfished biomass 
   r=exp(log_r); 
   bmin=exp(log_bmin);
   B0=exp(log_B0);
   // non centered version for random effects
   for(i in 1:R3){
     u3[i] = sigma_u .* u_tilde3[i];//non-centered
   }
   //reserve component
   for (i in 1:res){ 
     k[i] = exp(log(B0) +beta[1] * op[i]+beta[2] * sst[i] +beta[3]*hc[i]+beta[4]*at[i]); 
     mu[i] =    log(k[i]*(bmin/k[i])^exp(r*ag[i]*(p - 1)))+beta[12]*d[i]+beta[5]*rh_c[i]+beta[6]*rh_f[i]+beta[7]*rh_b[i]+beta[8]*cm_pc[i]+beta[9]*sa[i]+ beta[10]*si[i];
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
   u_tilde3 ~ normal(0,1); //prior re reserves
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
 BMMSY=B0/2.718281828;
 MMSY=((r*B0)/2.718281828); 
}

