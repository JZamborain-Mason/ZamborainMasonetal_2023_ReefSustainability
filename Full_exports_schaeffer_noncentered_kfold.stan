data {
 real<lower=0> p; //proportion of net exports
 int<lower=1> N; //number of reserve data points 
 int<lower=1> RE;//number of remote data points
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 real br[RE]; //Biomass response variable (y axes) in remote reefs (log-transformed)
 //explanatory variables
 int<lower=0> ag[N]; //age reserve
 real hc[N]; //predictor hard coral
 real si[N]; //predictor reserve size
 real d[N]; //predictor depth
 real op[N]; //predictor ocean productivity
 //explanantory variables remote
 real hc2[RE]; //predictor hard coral
 real d2[RE]; //predictor depth
 real op2[RE]; //predictor ocean productivity
 int<lower=1> R; //number of data regions for reserves(groups)
 int<lower=1, upper=R> pr[N]; //region id for reserves
 int<lower=1> R2; //number of data regions for remote(groups)
 int<lower=1, upper=R2> pr2[RE]; //region id for remote
 real newage[N];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
  int<lower=0,upper=1> holdout[N]; 
  //index whether the observation should be held out (1) or used (0)
 int<lower=0,upper=1> holdout2[RE]; 
  //index whether the observation should be held out (1) or used (0)
}
parameters {
 vector[4] beta; //slopes of explanatory variables 
 real<lower=0> sigma_e; //error sd
 real<lower=0> sigma_r; //error sd
 vector[R] u_tilde;
 real<lower=0> sigma_u; //region sd
 vector[R2] u2_tilde; 
 real<lower=0> sigma_u2; //region sd
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass (carrying capacity)
 }
transformed parameters {
 vector[N] mu;
 vector[N] k;
 vector[RE] mu2;
 vector[RE] k2;
 vector[R] u; //region intercepts
 vector[R2] u2;
  //reserve submodel
  for(i in 1:R){
    u[i] = sigma_u .* u_tilde[i];
  }
 for (i in 1:N) {
  k[i] = log(B0) +beta[1] * op[i] +beta[2]*d[i]+beta[3]*hc[i];
  mu[i] = log(exp(k[i])/(1+((exp(k[i])-bmin)/bmin)*exp(-(r-p)*ag[i])))+ beta[4]*si[i]+u[pr[i]];
 }
 //remote submodel
  for(i in 1:R2){
    u2[i] = sigma_u2 .* u2_tilde[i];
  }
 for (i in 1:RE) {
  k2[i]=log(B0)+beta[1]*op2[i]+ beta[2]*d2[i]+beta[3]*hc2[i];
  mu2[i] = log(exp(k2[i]))+u2[pr2[i]];
 }
}
model {
 //priors
 beta[1] ~ normal (0,10); //prior slope
 beta[2] ~ normal (0,10); //prior slope
 beta[3] ~ normal (0,10); //prior slope
 beta[4] ~ normal (0,10); //prior slope
 sigma_u2 ~ cauchy (0,1); //prior sd for varying intercept
 u2_tilde ~ normal (0, 1); //prior noncenteredversion
 sigma_u ~ cauchy (0,1); //prior sd for varying intercept (strong based on Gelman)
 u_tilde ~ normal(0, 1); //prior noncenteredversion
 r ~ normal (0.2,1); //weekly informative prior biomass intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 sigma_r ~ cauchy(0,1); //uninformative prior sd
 bmin ~ normal (0,200); //weakly informative prior lbiomass at age 0
 B0 ~ normal(120,200);//weakly informative prior for unfished biomass
  //likelihood holding out some data
    for(n in 1:N){
        if(holdout[n] == 0){
            target += normal_lpdf(b[n] | mu[n], sigma_e);
        }
    }
  for(n in 1:RE){
        if(holdout2[n] == 0){
            target += normal_lpdf(br[n] | mu2[n], sigma_r);
        }
    }
}
generated quantities {
 real <lower=0> BMMSY;
 real <lower=0> MMSY;
 real <lower=0> predreserveB[N];
 real <lower=0> predremoteB[RE];
 real  sustyield[B];
 vector[N+RE] log_lik;
 for (n in 1:N) {
 log_lik[n] = normal_lpdf(b[n]| mu, sigma_e);
 }
 for (n in (N+1):(N+RE)) {
 log_lik[n] = normal_lpdf(br[(n-N)]| mu2, sigma_r);
 }
//BMMSY and MMSY
 BMMSY=B0/2;
 MMSY=((r*B0)/4); 
 //predicted biomass in reserves following logistic growth
 for (i in 1:N) {
  predreserveB[i] = B0/(1+((B0-bmin)/bmin)*exp(-r*newage[i]));
 }
 //predicted biomass in remote reefs (i.e., unfished biomass estimate)
 for (i in 1:RE) {
  predremoteB[i] = B0;
 }
 //surplus production curve along a gradient of biomass
 for (i in 1:B){
  sustyield[i]=r* Bio[i] * (1-(Bio[i]/B0)); 
 }
}

