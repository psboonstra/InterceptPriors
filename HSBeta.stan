
data {
  int<lower = 0> n; // num obs
  int<lower = 0> p; // number of covariates (equipped with Cauchy prior)
  int<lower = 0,upper = 1> y[n]; // outcome
  matrix[n,p] x_standardized; //covariates (no intercept)
  real<lower = 0> beta_global_scale;
  real<lower = 0> beta_local_dof;
  real<lower = 0> beta_global_dof;
  real<lower = 0> slab_precision;
  int<lower = 0, upper = 2> alpha_prior_type;//0 = student-t; 1 = logistic; 2 = exponential power
  real alpha_mean;
  real<lower = 0> alpha_scale;
  real<lower = 0> alpha_power;
  real<lower = 0> alpha_dof;
  int<lower = 0, upper = 1> only_prior;//if 1, ignore the model and data
}
parameters {
  real alpha;
  vector[p] beta_raw;
  real<lower = 0> r1_global;//numerator for global shrinkage factor
  real<lower = 0> r2_global;//denominator for global shrinkage factor
  vector<lower = 0>[p] r1_local;//numerator for local shrinkage factor
  vector<lower = 0>[p] r2_local;//denominator for local shrinkage factor
}
transformed parameters {
  vector[p] beta;//shrunk regression coefficients
  real<lower=0> tau;//global shrinkage factor
  vector<lower=0>[p] lambda;//local shrinkage factor
  vector<lower=0, upper = sqrt(1/slab_precision)>[p] beta_total_scale;
  real alpha_copy; 
  if(alpha_scale > 0) {
    alpha_copy = alpha;
  } else {
    alpha_copy = alpha_mean;  
  } 
  tau = r1_global * sqrt(r2_global);
  lambda = r1_local .* sqrt(r2_local);
  beta_total_scale = 1 ./ sqrt(slab_precision + (1 ./ (beta_global_scale^2 * tau^2 * square(lambda))));
  beta =  (beta_total_scale .* beta_raw);
}
model {
  real half_beta_local_dof;
  real half_beta_global_dof;
  half_beta_local_dof = 0.5 * beta_local_dof;
  half_beta_global_dof = 0.5 * beta_global_dof;
  beta_raw ~ normal(0.0, 1.0);
  // half-t for lambdas
  r1_local ~ normal(0.0, 1.0);
  r2_local ~ inv_gamma(half_beta_local_dof, half_beta_local_dof);
  // half-t for tau
  r1_global ~ normal(0.0, 1.0);
  r2_global ~ inv_gamma(half_beta_global_dof, half_beta_global_dof);
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
    y ~ bernoulli_logit(alpha_copy + x_standardized * beta);
}

