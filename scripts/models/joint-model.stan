data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 vector[res] b;  // response variable (log-transformed)
 vector[fis] b2;  // response variable (log-transformed)
//explanatory variables for each component
 int mpa[res]; //age reserve (only for reserves)
 real ag[res]; //age reserve (only for reserves)
 real si[res]; //predictor reserve size (only for reserves)
 real env[res]; //predictor environment for reserves
 real env2[fis]; //predictor environment for not reserves
 int rest[fis]; //predictors of fishing restrictions for not reserves
}
parameters {
 vector[5] beta; //effect sizes
 real alpha; //intercept for reserves
 real<lower=0> sigma_r; //error sd for biomass reserves
 real<lower=0> sigma_f; //error sd for biomass fished
}
transformed parameters {
 vector[res] mu;//mean biomass  reserves
 vector[fis] mu2;//mean biomass fished
 
//reserve component
for (i in 1:res){ 
  mu[i] =  alpha + beta[1]*mpa[i] + beta[3]*ag[i] + beta[4]*si[i] + beta[5]*env[i];
}
//fished component
for (i in 1:fis){ 
  mu2[i] = alpha + beta[2]*rest[i] + beta[5]*env2[i];
  }
}
model {
 //priors
 beta[1] ~ normal (0,10); //prior slope
 beta[2] ~ normal (0,10); //prior slope
 beta[3] ~ normal (0,10); //prior slope
 beta[4] ~ normal (0,10); //prior slope
 beta[5] ~ normal (0,10); //prior slope
 sigma_r ~ cauchy(0,2.5); //uninformative prior sd
 sigma_f ~ cauchy(0,2.5); //uninformative prior sd
 alpha ~ normal(5,10);
 //likelihoods  
 for(n in 1:res){
      b[n] ~ normal(mu[n],sigma_r);
}
 for(n in 1:fis){
      b2[n] ~ normal(mu2[n],sigma_f);
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
 for (n in 1:res) {
 log_lik[n] = normal_lpdf(b[n]| mu[n], sigma_r);
}
for (n in 1:fis) {
 log_lik[n+res] = normal_lpdf(b2[n]| mu2[n], sigma_f);
}
}
