data {
 int<lower=1> res; //number of reserve  data points
 real b[res]; //Biomass response variable (y axes) (log-transformed)
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
 real gr[res]; //predictor grav
 int at[res];//predictor atoll
 int rh_b[res];//predictor backreef/lagoon habitat
 int rh_f[res];//predictor flat habitat
 int rh_c[res];//predictor crest habitat
 int cm_pc[res];//predictor census method point count
//fished
 real d3[fis]; //predictor depth 
 real sst3[fis]; //predictor sea surface temperature
 real sa3[fis]; //predictor sampling area
 real op3[fis]; //predictor ocean productivity 
 real gr3[fis]; //predictor ocean productivity 
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
}
parameters {
 vector[13] beta; //slopes for sampling/environmental variables 
 real<lower=0> sigma_e; //error sd for biomass reserves
 real<lower=0> sigma_f; //error sd for biomass fished
 vector[R3] u3; // random effects for fished regions
 real<lower=0> sigma_u; //region sd
 real log_r;//community biomass intrinsic growth rate
 real log_bmin; //biomass at reserve age 0
 real log_B0; //unfished biomass 
 real I_fished; //intercept for fished reefs
}
transformed parameters {
 vector[res] mu;//mean log-biomass  reserves
 vector[fis] mu3;//mean log-biomass  fished
 vector[res] k; //site-specific reserve carrying capacity given enevironmental factors
 real<lower=0> r;
 real<lower=0> bmin; 
 real<lower=0> B0; 
 r=exp(log_r); 
 bmin=exp(log_bmin); 
 B0=exp(log_B0);  
//reserve component
for (i in 1:res){ 
       k[i] = exp(log(B0) +beta[1] * op[i]+beta[2] * sst[i] +beta[3]*hc[i]+beta[4]*at[i]); 
       mu[i] = log(exp(log(k[i]*exp(log(bmin/k[i])*exp(-r*ag[i])))+beta[12]*d[i]+beta[5]*rh_c[i]+beta[6]*rh_f[i]+beta[7]*rh_b[i]+beta[8]*cm_pc[i]+beta[9]*sa[i]+ beta[10]*si[i]+ beta[13]*gr[i]));
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
 sigma_f ~ cauchy(0,1); //uninformative prior sd
 log_bmin ~ normal (log(10),1); //weakly informative prior reserve biomass at age 0
 log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
 I_fished~ normal(5,5); // prior for intercept in fished reefs
//likelihoods  
 for(n in 1:res){
            target += normal_lpdf(b[n] | mu[n], sigma_e);
    }
 for(n in 1:fis){
            target += normal_lpdf(b3[n] | mu3[n], sigma_f);
   }
}
