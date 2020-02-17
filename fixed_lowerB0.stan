data {
 int<lower=1> N; //number of reserve data points (observations)
 real b[N]; //Biomass response variable (y axes) sqrt
 //explanatory variables (all slopes and standrad belt transect)
 int<lower=0> ag[N]; //age reserve
 real B0;//fixed unfished biomass
 real newage[N];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
}
parameters {
 real<lower=0> sigma_e; //error sd
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 }
transformed parameters {
 vector[N] mu;
 for (i in 1:N) {
  mu[i] = log(B0/(1+((B0-bmin)/bmin)*exp(-r*ag[i])));
 }
}
model {
 //priors
  r ~ normal (0.2,1); //uninformative prior intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 bmin ~ normal (0,200); //weakly informative prior lbiomass at age 0
 // likelihood
 b~ normal(mu, sigma_e);
}
generated quantities {
 real <lower=0> MMSY;
 real <lower=0> BMMSY;
 real <lower=0> predreserveB[N];
 real  sustyield[B];
 //BMMSY and MMSY
 BMMSY=B0/2;
 MMSY=((r*B0)/4); 
 //predicted biomass in reserves
 for (i in 1:N) {
  predreserveB[i] = B0/(1+((B0-bmin)/bmin)*exp(-r*newage[i]));
 }
 for (i in 1:B){
  sustyield[i]=r* Bio[i] * (1-(Bio[i]/B0)); 
 }
}