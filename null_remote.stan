data {
 int<lower=1> RE;//number of remote data points
 real br[RE]; //Biomass response variable (y axes) in remote reefs (log-transformed)
 }
parameters {
 real<lower=0> sigma_r; //error sd
 real<lower=0> B0; //unfished biomass (carrying capacity)
 }
transformed parameters {
 vector[RE] mu2;
  //remote submodel
 for (i in 1:RE) {
  mu2[i] = log(B0);
 }
}
model {
 //priors
 sigma_r ~ cauchy(0,100); //uninformative prior sd
 B0 ~ normal(120,300);//weakly informative prior for unfished biomass
 // likelihood
 br~ normal (mu2,sigma_r);
}
