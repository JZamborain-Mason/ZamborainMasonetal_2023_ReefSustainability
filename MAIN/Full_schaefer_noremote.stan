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
 //to simulate from posterior
 real newage[res];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
 }
 parameters {
   vector[12] beta; //slopes for sampling/environmental variables 
   real<lower=0> sigma_e; //error sd for biomass reserves
   real<lower=0> sigma_f; //error sd for biomass fished
   vector[R3] u3; // random effects for fished regions
   real<lower=0> sigma_u; //region sd
   real log_r;//intrinsic growth rate
   real log_bmin; //biomass at reserve age 0
   real log_B0; //unfished biomass 
   real I_fished; //intercept for fished reefs
   real<lower=0,upper=1> p;// proportion of log biomass exported
 }
 transformed parameters {
   vector[res] mu;//mean log-biomass  reserves
   vector[fis] mu3;//mean log-biomass  fished
   vector[res] k; //site-specific reserve carrying capacity given enevironmental factors
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
   sigma_f ~ cauchy(0,1); //uninformative prior sd
   log_bmin ~ normal (log(40),1); //weakly informative prior reserve biomass at age 0
   log_B0 ~ normal(log(120),1);//weakly informative prior for unfished biomass
   I_fished~ normal(5,5); // prior for intercept in fished reefs
   p~ uniform (0,1); // uninformative prior for export proportion
   //likelihoods  
   for(n in 1:res){
     target += normal_lpdf(b[n] | mu[n], sigma_e);
   }
    for(n in 1:fis){
     target += normal_lpdf(b3[n] | mu3[n], sigma_f);
   }
 }
   generated quantities {
    real <lower=0> BMMSY; //BMMSY for averge and most common environemental conditions
   real <lower=0> MMSY;//MMSY for average and most common environmental conditions
   real <lower=0> predreserveB[res]; //expected reserve recovery for average and most common conditions
   real  sustyield[B]; //expected surplus production for average and most common conditions
   real site_B0[res+fis]; //site-specific B0 based on environmental conditions
   real site_BMMSY[res+fis];//site-specific BMMSY based on environmental conditions
   real site_MMSY[res+fis];//site-specific MMSY based on environmental conditions
   real site_Bmarg[res+fis];// marginalize biomass by sampling for status
   real site_Bmarg_all[res+fis];// marginalize biomass by sampling and environmental conditions (for ecosystem metrics)
   real site_status[res+fis]; //site-specific biomass status
   vector[res+fis] log_lik; //log-likelihood (for loo)
   for (n in 1:res) {
     log_lik[n] = normal_lpdf(b[n]| mu[n], sigma_e);
     site_Bmarg[n] = exp(b[n]-(beta[12]*d[n]+beta[5]*rh_c[n]+beta[6]*rh_f[n]+beta[7]*rh_b[n]+beta[8]*cm_pc[n]+beta[9]*sa[n]));
     site_Bmarg_all[n] = exp(b[n]-(beta[1] * op[n] +beta[2] * sst[n]+beta[3]*hc[n]+beta[4]*at[n]+beta[12]*d[n]+beta[5]*rh_c[n]+beta[6]*rh_f[n]+beta[7]*rh_b[n]+beta[8]*cm_pc[n]+beta[9]*sa[n]));
     site_B0[n] = exp(log(B0)+(beta[1] * op[n] +beta[2] * sst[n]+beta[3]*hc[n]+beta[4]*at[n]));
     site_BMMSY[n] = site_B0[n]/2;
     site_MMSY[n] = (site_B0[n]*r)/4;
     site_status[n] = site_Bmarg[n]/ site_BMMSY[n];
   }
     for (n in (res+1):(res+fis)) {
     log_lik[n] = normal_lpdf(b3[n-(res)]| mu3[n-(res)], sigma_f);
     site_Bmarg[n] = exp(b3[n-(res)]-(beta[12]*d3[n-(res)]+beta[8]*cm_pc3[n-(res)]+beta[9]*sa3[n-(res)]+beta[11]*cm_ds3[n-(res)]+beta[5]*rh_c3[n-(res)]+beta[6]*rh_f3[n-(res)]+beta[7]*rh_b3[n-(res)]));
     site_Bmarg_all[n] = exp(b3[n-(res)]-(beta[1] * op3[n-(res)] +beta[2] * sst3[n-(res)] +beta[3]*hc3[n-(res)]+beta[4]*at3[n-(res)]+beta[12]*d3[n-(res)]+beta[8]*cm_pc3[n-(res)]+beta[9]*sa3[n-(res)]+beta[11]*cm_ds3[n-(res)]+beta[5]*rh_c3[n-(res)]+beta[6]*rh_f3[n-(res)]+beta[7]*rh_b3[n-(res)]));
     site_B0[n] = exp(log(B0)+(beta[1] * op3[n-(res)] +beta[2] * sst3[n-(res)]+beta[3]*hc3[n-(res)]+beta[4]*at3[n-(res)]));
     site_BMMSY[n] = site_B0[n]/2; 
     site_MMSY[n] = (site_B0[n]*r)/4;
     site_status[n] = site_Bmarg[n]/ site_BMMSY[n];
   }
   //average and common conditions
   BMMSY=B0/2;
   MMSY=((r*B0)/4); 
   //predicted biomass in reserves following logistic growth
   for (i in 1:res) {
     predreserveB[i] =  (1-p)*exp(log(B0/(1+(( B0-bmin)/bmin)*exp(-r*newage[i]))));
   }
   //surplus production curve along a gradient of biomass for average/common categories
   for (i in 1:B){
     sustyield[i]=r* Bio[i] * (1-(Bio[i]/B0)); //common condition
   }
 }

