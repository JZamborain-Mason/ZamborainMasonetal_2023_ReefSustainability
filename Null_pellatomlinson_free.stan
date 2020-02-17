data {
 int<lower=1> N; //number of reserve data points 
 int<lower=1> RE;//number of remote data points
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 real br[RE]; //Biomass response variable (y axes) in remote reefs (log-transformed)
 int<lower=0> ag[N]; //age reserve 
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
 real n;//parameter that modifies at which biomass the surplus production peaks 
 }
transformed parameters {
 real y;//parameter from fletcher to re-parametrize the pella-tomlinson
 real m;//equivalent to MMSY
 vector[N] mu;
 vector[RE] mu2;
 y=(n^(n/(n-1)))/(n-1);
 m=(r*B0)*(1/(1+(n/2)))^((1/(n/2))+1);
 //reserve submodel
 for (i in 1:N) {
  mu[i] = log((B0^(1-n)+(bmin^(1-n)-B0^(1-n))*exp(((y*m)/B0)*(1-n)*ag[i]))^(1/(1-n)));
 }
 //remote submodel
 for (i in 1:RE) {
  mu2[i] = log(B0);
 }
}
model {
 //priors
 r ~ normal (0.2,1); //weekly informative prior biomass intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 sigma_r ~ cauchy(0,1); //uninformative prior sd
 bmin ~ normal (0,200); //weakly informative prior lbiomass at age 0
 B0 ~ normal(120,200);//weakly informative prior for unfished biomass
 n~normal(2,2);
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
 for (i in 1:N) {
 log_lik[i] = normal_lpdf(b[i]| mu, sigma_e);
 }
 for (i in (N+1):(N+RE)) {
 log_lik[i] = normal_lpdf(br[(i-N)]| mu2, sigma_r);
 }
 //BMMSY and MMSY
 uMMSY=r*(1/(1+(n/2)));
 BMMSY=n^(1/(1-n))*B0;
  MMSY=(r*B0)*(1/(1+(n/2)))^((1/(n/2))+1); //from merging FAO and QUinn and deriso we get that the FAO p is n/2 of the quinn and deriso
 //predicted biomass in reserves following logistic growth
 for (i in 1:N) {
  predreserveB[i] =(B0^(1-n)+(bmin^(1-n)-B0^(1-n))*exp(((y*((r*B0)/4))/B0)*(1-n)*newage[i]))^(1/(1-n));
 }
 //predicted biomass in remote reefs (i.e., unfished biomass estimate)
 for (i in 1:RE) {
  predremoteB[i] = B0;
 }
 //surplus production curve along a gradient of biomass
 for (i in 1:B){
  sustyield[i]=(y*MMSY* (Bio[i]/B0))- (y*MMSY* (Bio[i]/B0)^n); 
 }
}
