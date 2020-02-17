data {
 int<lower=1> N; //number of reserve data points 
 int<lower=1> RE;//number of remote data points
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 real br[RE]; //Biomass response variable (y axes) in remote reefs (log-transformed)
 int<lower=0> ag[N]; //age reserve 
 real si[N]; //age reserve
 real newage[N];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
}
parameters {
 real<lower=0> sigma_e; //error sd
 real<lower=0> sigma_r; //error sd
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass (carrying capacity)
 real<upper=0> slope; //slope of the relation between spillover and reserve size. as reserves get bigger they export less
}
transformed parameters {
 vector[N] mu;
 vector[N] p;//proportion of exports (spillover)
 vector[RE] mu2;
 //reserve submodel
 for (i in 1:N) {
 p[i]=slope*si[i];
  mu[i] = log(B0/(1+((B0-bmin)/bmin)*exp(-(r-p[i])*ag[i])));
 }
 //remote submodel
 for (i in 1:RE) {
  mu2[i] = log(B0);
 }
}
model {
 //priors
 slope~ normal(-0.002,1);
 r ~ normal (1,3); //weekly informative prior biomass intrinsic growth rate
 sigma_e ~ cauchy(0,100); //uninformative prior sd
 sigma_r ~ cauchy(0,100); //uninformative prior sd
 bmin ~ normal (0,300); //weakly informative prior lbiomass at age 0
 B0 ~ normal(120,300);//weakly informative prior for unfished biomass
 // likelihood
 b~ normal(mu, sigma_e);
 br~ normal (mu2,sigma_r);
}
generated quantities {
 real <lower=0> BMMSY;
 real <lower=0> MMSY;
 real <lower=0> uMMSY;
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
 uMMSY=r/2;
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