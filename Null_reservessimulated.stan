data {
 int<lower=1> N; //number of reserve data points 
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 int<lower=0> ag[N]; //age reserve 
 real newage[N];//added the seq reserve age to predict
 }
parameters {
 real<lower=0> sigma;
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass (carrying capacity)
 }
transformed parameters {
 vector[N] mu;
 //reserve submodel
 for (i in 1:N) {
  mu[i] = B0/(1+((B0-bmin)/bmin)*exp(-r*ag[i]));
 }
}
model {
//priors
r~normal(0.2,1);
B0~normal (3000,1000);
bmin~normal (0,1000);
sigma~cauchy(1,100);
// likelihood
 b~ normal(mu, sigma);
}
generated quantities {
 real <lower=0> predreserveB[N];
  //predicted biomass in reserves following logistic growth
 for (i in 1:N) {
  predreserveB[i] = B0/(1+((B0-bmin)/bmin)*exp(-r*newage[i]));
 }
}