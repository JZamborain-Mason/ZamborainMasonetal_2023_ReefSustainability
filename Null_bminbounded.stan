data {
 int<lower=1> N; //number of reserve data points 
 int<lower=1> RE;//number of remote data points
 int<lower=1> FI;//number of fished data points contributing to bmin
 real b[N]; //Biomass response variable (y axes) in reserves(log-transformed)
 real br[RE]; //Biomass response variable (y axes) in remote reefs (log-transformed)
 real bf[FI]; //Biomass response variable (y axes) in FISHED reefs (log-transformed)
 int<lower=0> ag[N]; //age reserve 
 real newage[N];//added the seq reserve age to predict
 int<lower=1> B; //number of BIOMASS data points to predict popgrowth
 real<lower=0> Bio[B];//given Biomass to calculate sustainable yield at a biomass gradient
 }
parameters {
 real<lower=0> sigma_e; //error sd
 real<lower=0> sigma_r; //error sd
 real<lower=0> sigma_f; //error sd
 real<lower=0> r;//intrinsic growth rate
 real<lower=0> bmin; //biomass at reserve age 0
 real<lower=0> B0; //unfished biomass (carrying capacity)
 }
transformed parameters {
 vector[N] mu;
 vector[RE] mu2;
 vector[FI] mu3;
  //reserve submodel
 for (i in 1:N) {
  mu[i] = log(B0/(1+((B0-bmin)/bmin)*exp(-r*ag[i])));
 }
 //remote submodel
 for (i in 1:RE) {
  mu2[i] = log(B0);
 }
  //bmin fished submodel
 for (i in 1:FI) {
  mu3[i] = log(bmin);
 }
}
model {
 //priors
 r ~ normal (0.2,1); //weekly informative prior biomass intrinsic growth rate
 sigma_e ~ cauchy(0,1); //uninformative prior sd
 sigma_r ~ cauchy(0,1); //uninformative prior sd
 sigma_f ~ cauchy(0,1); //uninformative prior sd
 bmin ~ normal (0,200); //weakly informative prior lbiomass at age 0
 B0 ~ normal(120,200);//weakly informative prior for unfished biomass
 // likelihood
 b~ normal(mu, sigma_e);
 br~ normal (mu2,sigma_r);
 bf~ normal (mu3,sigma_f);
 }
generated quantities {
 real <lower=0> BMMSY;
 real <lower=0> MMSY;
 real <lower=0> predreserveB[N];
 real <lower=0> predremoteB[RE];
 real  sustyield[B];
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