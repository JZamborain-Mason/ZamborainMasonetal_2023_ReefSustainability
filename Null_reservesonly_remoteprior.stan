data {
//use as prior for B0 the posterior from remote locations
 real sigma_B0; 
 real mean_B0;
 int<lower=1> N; //number of reserve data points 
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 int<lower=0> ag[N]; //age reserve 
 real newage[N];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
 }
parameters {
 real<lower=0> sigma_e; //error sd
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass (carrying capacity)
 }
transformed parameters {
 vector[N] mu;
 //reserve submodel
 for (i in 1:N) {
  mu[i] = log(B0/(1+((B0-bmin)/bmin)*exp(-r*ag[i])));
 }
}
model {
 //priors
 r ~ normal (0.2,1); //weekly informative prior biomass intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 bmin ~ normal (0,200); //weakly informative prior lbiomass at age 0
 B0 ~ normal(mean_B0,sigma_B0);//weakly informative prior for unfished biomass
 // likelihood
 b~ normal(mu, sigma_e);
 }
generated quantities {
 real <lower=0> BMMSY;
 real <lower=0> MMSY;
 real <lower=0> predreserveB[N];
 real  sustyield[B];
 //BMMSY and MMSY
 BMMSY=B0/2;
 MMSY=((r*B0)/4); 
 //predicted biomass in reserves following logistic growth
 for (i in 1:N) {
  predreserveB[i] = B0/(1+((B0-bmin)/bmin)*exp(-r*newage[i]));
 }
  //surplus production curve along a gradient of biomass
 for (i in 1:B){
  sustyield[i]=r* Bio[i] * (1-(Bio[i]/B0)); 
 }
}