data {
  int<lower=0> n; // num obs
  int<lower=0> p; // number of covariates (equipped with Cauchy prior)
  int<lower=0,upper=1> y[n]; // outcome
  matrix[n,p] x_centered; //covariates (no intercept)
  real<lower=0> beta_expit_shape;
  int<lower = 0, upper = 2> alpha_prior_type;//0 = student-t; 1 = logistic; 2 = exponential power
  real alpha_mean;
  real<lower = 0> alpha_scale;
  real<lower = 0> alpha_power;
  real<lower = 0> alpha_dof;
  int<lower=0,upper=1> only_prior;//if 1, ignore the model and data
}
parameters {
  real alpha;
  vector<lower=0,upper=1>[p] beta_expit;
}
transformed parameters {
  real alpha_copy; 
  vector[p] beta;//shrunk regression coefficients
  if(alpha_scale > 0) {
    alpha_copy = alpha;
  } else {
    alpha_copy = alpha_mean;  
  } 
  beta = logit(beta_expit);
}
model {
  beta_expit ~ beta(beta_expit_shape,beta_expit_shape);
  if(alpha_prior_type == 0 && alpha_scale > 0) {
    alpha ~ student_t(alpha_dof, alpha_mean, alpha_scale);
  } else if(alpha_prior_type == 1 && alpha_scale > 0) {
    alpha ~ logistic(alpha_mean, alpha_scale);
  } else if(alpha_prior_type == 2 && alpha_scale > 0) {
    target += -fabs((alpha-alpha_mean)/(sqrt(2.0)*alpha_scale))^(alpha_power);
  } else {
    alpha ~ normal(0.0, 1.0);
  }
  if(only_prior == 0) 
    y ~ bernoulli_logit(alpha_copy + x_centered * beta);
}

