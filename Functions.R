
#DESCRIPTION: self-explanatory
expit = function(x) { 1/(1+exp(-x));}
logit = function(x) { log(x/(1-x));}

#DESCRIPTION: Error-handling function
#Borrowed from the R package simsalapar 
#https://www.rdocumentation.org/packages/simsalapar/versions/1.0-9/topics/tryCatch.W.E
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W,w)
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


#DESCRIPTION: Program for fitting a (logistic) GLM equipped with a hierarchical shrinkage prior on beta
#and any of the priors on alpha considered in Boonstra, et al.. It has two intended uses: 
#compile stan scripts or analyze data. First, if the user provides nothing but a valid 
#'stan_path', then the stan script is compiled. Second, the user provides both a compiled 
#stanfit object as well as values for y, x_standardized, and any other desired 
#arguments to actually fit a regression. 
#
#
#ARGUMENTS:
#
#stan_fit: an R object of class stanfit, which allows the function to run without recompiling the stan code.
#
#stan_path: (character) a path pointing to a .stan file, which indicates the stan code to compile and run. If
#both stan_fit and stan_path are provided, stan_fit takes precedence. 
#
#y (vector) outcomes corresponding to the type of glm desired. This should match whatever datatype is expected 
#by the stan program.
#
#x_standardized (matrix) matrix of numeric values with number of rows equal to the length of y and number of columns
#equal to p+q. It is assumed without verification that each column is standardized to whatever scale the prior 
#expects - in Boonstra and Barbaro, all predictors are marginally generated to have mean zero and unit variance, so no 
#standardization is conducted. In practice, all data should be standardized to have a common scale before model fitting. 
#If regression coefficients on the natural scale are desired, they be easily obtained through unstandardizing. 
#
#beta_global_scale (pos. real) constants indicating the prior scale of the horseshoe. Corresponds to
#to 'psi_n' in the notation of Boonstra, et al.
#
#local_dof, global_dof (pos. integer) numbers indicating the degrees of freedom for lambda_j and tau, respectively. Boonstra, 
#et al. never considered local_dof != 1 or global_dof != 1. 
#
#slab_precision (pos. real) the slab-part of the regularized horseshoe, this is equivalent to (1/d)^2 in the notation of
#Boonstra and Barbaro 
#
#alpha_prior_args (list)
#
#only_prior (logical) should all data be ignored, sampling only from the prior?
#
#ntries (pos. integer) the stan function will run up to this many times, stopping either when the number of 
#*divergent transitions* is zero or when ntries has been reached. The reported fit will be that with the fewest number of 
#divergent iterations. 

glm_HSBeta = function(stan_fit = NA, 
                      stan_path,
                      y = c(0,1),
                      x_standardized = matrix(0,length(y),6), 
                      beta_global_scale = 1, 
                      local_dof = 1, 
                      global_dof = 1, 
                      slab_precision = (1/15)^2, 
                      alpha_prior_args = list(type = 0, 
                                              mean = 0, 
                                              scale = 10, 
                                              power = 0,
                                              dof = 3),
                      only_prior = F, 
                      mc_warmup = 50, 
                      mc_iter_after_warmup = 50, 
                      mc_chains = 1, 
                      mc_thin = 1, 
                      mc_stepsize = 0.1, 
                      mc_adapt_delta = 0.9,
                      mc_max_treedepth = 15,
                      ntries = 1) {
  
  p = ncol(x_standardized);
  n = nrow(x_standardized);
  stopifnot(n == length(y));
  
  alpha_prior_type = alpha_prior_args$type;#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
  alpha_mean = alpha_prior_args$mean;
  alpha_scale = alpha_prior_args$scale;
  alpha_power = alpha_prior_args$power;#ignored unless alpha_prior_type = 2;
  alpha_dof = alpha_prior_args$dof;#ignored unless alpha_prior_type = 0;
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  
  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n = n,
                                                    p = p,
                                                    y = y,
                                                    x_standardized = x_standardized,
                                                    beta_global_scale = beta_global_scale,
                                                    beta_local_dof = local_dof,
                                                    beta_global_dof = global_dof,
                                                    slab_precision = slab_precision,
                                                    alpha_prior_type = alpha_prior_type,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                                    alpha_mean = alpha_mean,
                                                    alpha_scale = alpha_scale,
                                                    alpha_power = alpha_power,#ignored unless alpha_prior_type = 2;
                                                    alpha_dof = alpha_dof,#ignored unless alpha_prior_type = 0;
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below. 
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_alpha = as.numeric(foo$alpha);
      curr_beta = foo$beta;
      curr_theta = foo$beta_total_scale;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_alpha = curr_alpha,
         curr_beta = curr_beta,
         curr_theta = curr_theta);
  }
}

glm_LogisBeta = function(stan_fit = NA, 
                         stan_path,
                         y = c(0,1),
                         x_centered = matrix(0,length(y),6), 
                         beta_expit_shape = 1,
                         alpha_prior_args = list(type = 0, 
                                                 mean = 0, 
                                                 scale = 10, 
                                                 power = 0,
                                                 dof = 3),
                         only_prior = F, 
                         mc_warmup = 50, 
                         mc_iter_after_warmup = 50, 
                         mc_chains = 1, 
                         mc_thin = 1, 
                         mc_stepsize = 0.1, 
                         mc_adapt_delta = 0.9,
                         mc_max_treedepth = 15,
                         ntries = 1) {
  
  p = ncol(x_centered);
  n = nrow(x_centered);
  stopifnot(n == length(y));
  
  alpha_prior_type = alpha_prior_args$type;#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
  alpha_mean = alpha_prior_args$mean;
  alpha_scale = alpha_prior_args$scale;
  alpha_power = alpha_prior_args$power;#ignored unless alpha_prior_type = 2;
  alpha_dof = alpha_prior_args$dof;#ignored unless alpha_prior_type = 0;
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  
  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n = n,
                                                    p = p,
                                                    y = y,
                                                    x_centered = x_centered,
                                                    beta_expit_shape = beta_expit_shape,
                                                    alpha_prior_type = alpha_prior_type,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                                    alpha_mean = alpha_mean,
                                                    alpha_scale = alpha_scale,
                                                    alpha_power = alpha_power,#ignored unless alpha_prior_type = 2;
                                                    alpha_dof = alpha_dof,#ignored unless alpha_prior_type = 0;
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below. 
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_alpha = as.numeric(foo$alpha);
      curr_beta = foo$beta;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_alpha = curr_alpha,
         curr_beta = curr_beta);
  }
}



draw_data = function(n = 100,
                     n_new = 100,
                     true_alpha = 0,
                     true_betas = 0,
                     chol_latent_x = 1, 
                     mu_latent_x = 0,
                     x_binom = 1,
                     case_control_ratio = NA,
                     seed = .Machine$integer.max
) {
  set.seed(seed);
  p = length(true_betas);
  n_by_mu_latent_x = matrix(mu_latent_x,nrow=n,ncol=p,byrow=T);
  n_new_by_mu_latent_x = matrix(mu_latent_x,nrow=n_new,ncol=p,byrow=T);
  
  if(is.na(case_control_ratio)) {#Cohort Sampling
    x = matrix(rnorm(n*p),nrow=n)%*%chol_latent_x + n_by_mu_latent_x;
    x_new = matrix(rnorm(n_new*p),nrow=n_new)%*%chol_latent_x + n_new_by_mu_latent_x;
    
    if(length(x_binom)>0) {
      x[,x_binom] = 1*(x[,x_binom,drop = F]>0);
      x_new[,x_binom] = 1*(x_new[,x_binom,drop = F]>0);
    }
    y = rbinom(n,1,expit(drop(true_alpha+x%*%true_betas)));
    y_new = rbinom(n_new,1,expit(drop(true_alpha+x_new%*%true_betas)));
  } else {#Case-Control Sampling
    n_case = ceiling(n/(1+1/case_control_ratio));
    n_control = n - n_case;
    x_case = matrix(0,n_case,p);
    x_control = matrix(0,n_control,p);
    curr_case = curr_control = 0;
    while((curr_case<n_case)|(curr_control<n_control)) {
      foo = matrix(rnorm(n*p),nrow=n)%*%chol_latent_x + n_by_mu_latent_x;
      if(length(x_binom)>0) {
        foo[,x_binom] = 1*(foo[,x_binom,drop = F]>0);
      }
      foo2 = rbinom(n,1,expit(drop(true_alpha+foo%*%true_betas)));
      if(curr_case<n_case) {
        x_case[curr_case+(1:min(n_case-curr_case,sum(foo2))),] = foo[which(foo2==1)[1:min(n_case-curr_case,sum(foo2))],]
        curr_case = curr_case + min(n_case-curr_case,sum(foo2));
      }
      if(curr_control<n_control) {
        x_control[curr_control+(1:min(n_control-curr_control,sum(1-foo2))),] = foo[which(foo2==0)[1:min(n_control-curr_control,sum(1-foo2))],]
        curr_control = curr_control + min(n_control-curr_control,sum(1-foo2));
      }
    }
    x = rbind(x_case,x_control);
    y = c(rep(1,n_case),rep(0,n_control));
    
    x_new = matrix(rnorm(n_new*p),nrow=n_new)%*%chol_latent_x + n_new_by_mu_latent_x;
    if(length(x_binom)>0) {
      x_new[,x_binom] = 1*(x_new[,x_binom,drop = F]>0);
    }
    y_new = rbinom(n_new,1,expit(drop(true_alpha+x_new%*%true_betas)));
    
  }
  list(x = x, 
       y = y, 
       x_new = x_new,
       y_new = y_new);
}

#Main function to run a single simulation 

InterceptSim <- function(niter,#Number of simulations
                         n,#size of training per simulation
                         n_new,#size of validation per simulation
                         true_alpha,#generating intercept parameter
                         true_betas,#generating regression coefficients
                         x_binom,#vector of logicals to indicate which x's are binomial
                         case_control_ratio=NULL,#
                         mu_latent_x,
                         chol_latent_x,
                         priors_to_fit = NULL,
                         hiershrink_prior_num_relevant,
                         local_dof = 1,
                         global_dof = 1,
                         slab_precision = (1.0/15.0)^2,
                         exppow_prior_args,
                         stan_file_path, 
                         hsbeta_file_name = "HS_VariedIntPriors.stan",
                         logisbeta_file_name = "Logis_VariedIntPriors.stan",
                         calc_n_over = T,#should n overlap be calculated (set to F to decrease runtime)
                         calc_n_comp = T,#same for ncomplete
                         calc_n_piv = T,#same for npviot
                         algorithm2_maxk = 10,#arguments to algorithm 2
                         algorithm2_ndraws = 100,#arguments to algorithm 2
                         mc_chains, #number of chains
                         mc_iter_after_warmup,# number of total iterations, including warmup
                         mc_warmup,#number of warmup iterations
                         mc_thin,
                         mc_max_treedepth,
                         mc_adapt_delta,
                         mc_stepsize,
                         ntries_per_iter = 4,
                         random_seed,
                         array_id,
                         scenario,
                         fit_methods = T, 
                         #If both x and y are provided, a single dataset will be analyzed
                         data_analysis = F,#logical flag to indicate a data analysis instead of simulation
                         x_all = NULL, #design matrix data analysis
                         y_all = NULL,#outcomes
                         frac_training = NULL,#fraction of data to designate as validation
                         num_partitions = 100,#number of random partitions to implement
                         prespecified_training_subsets = NULL,#matrix, with each row indicating the observations to designate as training for that iteration
                         verbose = F#return lots of information from data analysis?
) {
  set.seed(random_seed);
  data_seeds = sample(.Machine$integer.max,niter);#Mersenne-Twister is seeded by 32 bit integers 
  
  if(data_analysis && (niter != num_partitions)) {
    stop("num_partitions must equal niter if data_analysis == T");
  }
  
  if(missing(logisbeta_file_name) || is.null(logisbeta_file_name)) {
    do_logisbeta_prior = F;
    warning("skipping logistic(1) prior on beta because no stan file was provided");
  } else {
    do_logisbeta_prior = T;
  }
  
  if(is.null(priors_to_fit)) {
    intercept_prior_names = c(
      "StudT10",
      "Normal10");
  } else {
    intercept_prior_names = priors_to_fit;
  }
  stopifnot(all(names(exppow_prior_args)==grep("EP",intercept_prior_names,value=T)));
  
  method_names = paste0(rep(c("HS_","Logis_"),each=length(intercept_prior_names)),rep(intercept_prior_names,times=2));
  
  
  if(data_analysis) {
    stopifnot(!is.null(x_all)&!is.null(y_all));
    stopifnot(is.null(prespecified_training_subsets) || ((nrow(prespecified_training_subsets) == niter) && (ncol(prespecified_training_subsets) <= nrow(x_all)))); 
    
    store_alpha  = 
      store_beta = 
      store_fitted_probs = vector("list",length(method_names));
    names(store_alpha) = 
      names(store_beta) = 
      names(store_fitted_probs) = method_names;
    if(!is.null(prespecified_training_subsets)) {
      if(!is.null(frac_training) && (frac_training * nrow(x_all) != ncol(prespecified_training_subsets) )) {warning("overriding given value of 'frac_training' because prespecified training subsets were provided")}
      frac_training = ncol(prespecified_training_subsets) / nrow(x_all) ;
    } 
    n = round(frac_training * nrow(x_all));
    p = ncol(x_all);
    
    #Placeholder values that are not used;
    true_alpha = 0;
    true_betas = numeric(p)
    x_binom = NULL;
    case_control_ratio = NULL;
    mu_latent_x = numeric(p);
    chol_latent_x = diag(1,p);
    if(nrow(x_all) > n) {#Check if frac_training == 1 (within rounding error);
      n_new = nrow(x_all) - n;
      niter = num_partitions;
    } else {
      n_new = n;
      niter = 1;
    }
  } 
  
  n = as.integer(n);
  p = length(true_betas);
  stopifnot(nrow(chol_latent_x)==p & ncol(chol_latent_x) == p);
  if(is.null(mu_latent_x)) {
    mu_latent_x = numeric(p);
  } else {
    stopifnot(length(mu_latent_x)==p)
  }
  #Global scale parameter to use for hierarchical shrinkage prior
  if(fit_methods) {
    stopifnot(hiershrink_prior_num_relevant < p);
    beta_global_scale = solve_for_hiershrink_scale(hiershrink_prior_num_relevant,
                                                   local_dof = local_dof,
                                                   global_dof = global_dof,
                                                   npar = p,
                                                   n = n,
                                                   sigma = 2)$tau0;
  } else {
    beta_global_scale = NA;
  }
  #
  
  colnames_x = paste("x",1:p,sep = "");
  mc_total_samps = mc_chains * mc_iter_after_warmup / mc_thin;
  
  matrix_true_betas = matrix(true_betas,
                             nrow = mc_total_samps,
                             ncol = p,
                             byrow=T);
  
  
  #Store average empirical covariance matrix of x_all
  avg_cov_obs_x = matrix(0,p,p);
  store_obs_mean_x = store_obs_sd_x = matrix(NA,nrow=niter,ncol=p);
  
  ##Store results for alpha
  store_mean_alpha  = matrix(NA,niter,length(method_names),dimnames = list(NULL,method_names));
  ##Store results for beta
  store_mean_beta = store_sd_beta = store_50CIcoverage_beta = store_50CIwidth_beta = vector("list",length(method_names));
  names(store_mean_beta) = names(store_sd_beta) = names(store_50CIcoverage_beta) = names(store_50CIwidth_beta) = method_names;
  for(i in method_names) {
    store_mean_beta[[i]] =  store_sd_beta[[i]] = store_50CIcoverage_beta[[i]] = store_50CIwidth_beta[[i]] =
      matrix(NA,niter,p,dimnames=list(NULL,true_betas));
  }
  ##store expected log posterior density (overall and for events only)
  elpd = elpd_events = elpd_nonevents = store_auc = bayes_rmse = matrix(NA,niter,length(method_names),dimnames=list(NULL,method_names));
  ##Population-prevalence
  population_prevalence = matrix(NA,nrow=niter,ncol=12,dimnames = list(NULL,c("avg_true_prob",
                                                                              "sd_true_prob",
                                                                              "avg_true_xbeta",
                                                                              "sd_true_xbeta",
                                                                              "avg_emp_prob",
                                                                              "extreme_true_prob",
                                                                              "min_true_prob",
                                                                              "lowerQ_true_prob",
                                                                              "middleQ_true_prob",
                                                                              "upperQ_true_prob",
                                                                              "max_true_prob",
                                                                              "optimal_auc"
  )));
  
  algorithm2_results = matrix(NA,niter,3,dimnames = list(NULL,c("n_over","n_comp","n_piv")));
  ##Store information on numerical challenging xs
  store_attributes_data = matrix(0,nrow=niter,ncol=5,dimnames = list(1:niter,c("const_y","const_x","copy_x","lt_full_rank","divergence")));
  max_divergences_by_prior =
    accepted_divergences_by_prior = 
    max_rhat_by_prior = matrix(NA,niter,length(method_names),dimnames = list(NULL,method_names));
  
  ##Store datasets that require very small step sizes
  interesting_datasets = list();
  
  #The value of the intercept such that, when sum(y)=0, the likelihood ratio of this value of the maximum likelihood estimate (i.e. infinity) is 1;
  #This is equal to logit(s_n) in the notation of the paper.
  logit_sn = -log(exp(1/2/n)-1);
  #Identify adaptive parameters in novel priors, which are functions of n by way of the likelihood_cutpoints
  if(any(grepl("EP",intercept_prior_names))) {
    for(k in 1:length(exppow_prior_args)) {
      curr_prior = names(exppow_prior_args)[k];
      if(curr_prior%in%intercept_prior_names) {
        if(exppow_prior_args[[k]]$adaptive_scale) {
          exppow_prior_args[[k]]$alpha_scale = solve_for_exppow_scale(logit_sn,0.01,exppow_prior_args[[k]]$alpha_power);
        }
      } 
    }
    rm(curr_prior);
  }
  if(any(grepl("LogisticAdaptive",intercept_prior_names))) {
    logistic_adaptive_alpha_scale = solve_for_logistic_scale(logit_sn,0.01);
  } else {
    logistic_adaptive_alpha_scale = NULL;
  }
  
  sd_actual_x = round(sqrt(diag(crossprod(chol_latent_x))),15);
  mu_actual_x = mu_latent_x;
  if(length(x_binom)>0) {
    mu_actual_x[x_binom] = pnorm((mu_latent_x/sd_actual_x)[x_binom]);
    sd_actual_x[x_binom] = sqrt(mu_actual_x[x_binom]*(1 - mu_actual_x[x_binom]));
  }
  
  i=1;  
  begin = Sys.time();
  stan_compiled=F;
  
  for(i in 1:niter) {
    curr_try = 1;
    dataset_stored = F;
    ## Data setup====##########################################################################
    if(!data_analysis) {
      #If data_analysis is FALSE, then this is a simulation study
      foo = draw_data(n = n,
                      n_new = n_new,
                      true_alpha = true_alpha,
                      true_betas = true_betas,
                      chol_latent_x = chol_latent_x, 
                      mu_latent_x = mu_latent_x,
                      x_binom = x_binom,
                      case_control_ratio = case_control_ratio,
                      seed = data_seeds[i]);
      x = foo$x;
      y = foo$y;
      x_new = foo$x_new;
      y_new = foo$y_new;
      rm(foo);
    } else if(!is.null(prespecified_training_subsets)){
      #Prespecified training/testing partitioning of the data
      curr_samp = prespecified_training_subsets[i,];
      if(!length(setdiff(curr_samp,1:n)) && !length(setdiff(1:n,curr_samp))) {
        x = x_new = x_all;
        y = y_new = y_all;
      } else {
        x = x_all[curr_samp,,drop = F];
        y = y_all[curr_samp];
        x_new = x_all[-curr_samp,,drop = F];
        y_new = y_all[-curr_samp];
      }
    } else if(nrow(x_all) > n) {
      #If training/testing partition is not prespecified but the total size of the data exceeds the value of n, 
      #then do random partitioning of the provided data
      curr_samp = sort(sample(nrow(x_all),n));
      x = x_all[curr_samp,,drop = F];
      y = y_all[curr_samp];
      x_new = x_all[-curr_samp,,drop = F];
      y_new = y_all[-curr_samp];
    } else {
      #Otherwise, use the provided data as both training and testing (which will of course be optimistic)
      x = x_new = x_all;
      y = y_new = y_all;
    }
    
    population_prevalence[i,"avg_emp_prob"] = mean(y);
    if(!data_analysis) {
      curr_true_probs = expit(drop(true_alpha+(foo<-x_new%*%true_betas)));
      population_prevalence[i,"avg_true_xbeta"] = mean(foo);
      population_prevalence[i,"sd_true_xbeta"] = sd(foo);
      population_prevalence[i,"avg_true_prob"] = mean(curr_true_probs);
      population_prevalence[i,"extreme_true_prob"] = mean(pmin(curr_true_probs,1-curr_true_probs));
      population_prevalence[i,"sd_true_prob"] = sd(curr_true_probs);
      population_prevalence[i,"lowerQ_true_prob"] = quantile(curr_true_probs,p = 0.25);
      population_prevalence[i,"middleQ_true_prob"] = quantile(curr_true_probs,p = 0.50);
      population_prevalence[i,"upperQ_true_prob"] = quantile(curr_true_probs,p = 0.75);
      population_prevalence[i,"min_true_prob"] = min(curr_true_probs);
      population_prevalence[i,"max_true_prob"] = max(curr_true_probs);
      if(diff(range(y_new)) > .Machine$double.eps^0.5) {
        population_prevalence[i,"optimal_auc"] = pROC::auc(pROC::roc(y_new, curr_true_probs));
      }
      rm(foo);
    }
    
    x = data.matrix(x);
    colnames(x) = colnames(x_new) = colnames_x;
    obs_mean_x = colMeans(x);
    matrix_obs_mean_x = matrix(obs_mean_x,nrow=mc_total_samps,ncol=p,byrow=T);
    obs_sd_x = sqrt(colMeans(x^2)-obs_mean_x^2);
    avg_cov_obs_x = avg_cov_obs_x + stats::cov(x_new)/niter; 
    if(any(obs_sd_x < .Machine$double.eps^0.5)) {#Constant valued columns
      obs_sd_x[curr_const <- which(obs_sd_x < .Machine$double.eps^0.5)] = Inf;
      store_attributes_data[i,"const_x"] = 1;
    }
    if((qr_x<-qr(crossprod(x)))$rank < p) {#LT Full Rank, copied columns
      store_attributes_data[i,"lt_full_rank"] = 1;
      if(sum(as.numeric(cor_x<-abs(cor(x)*lower.tri(avg_cov_obs_x)) > 1 - .Machine$double.eps^0.5),na.rm=T)) {
        store_attributes_data[i,"copy_x"] = 1;  
        obs_sd_x[unique(unlist(lapply(strsplit(names(unlist(apply(cor_x,2,which))),".",fixed=T),"[",2)))] = Inf;
      }
    }
    if(diff(range(y)) < .Machine$double.eps^0.5) {
      store_attributes_data[i,"const_y"] = 1;
      store_attributes_data[i,"divergence"] = NA;
      next;
    }
    store_obs_sd_x[i,] = obs_sd_x;
    store_obs_mean_x[i,] = obs_mean_x;
    matrix_obs_sd_x = matrix(obs_sd_x,nrow=mc_total_samps,ncol=p,byrow=T);
    x_standardized = scale(x,center=obs_mean_x,scale=obs_sd_x);
    x_centered = scale(x,center=obs_mean_x,scale=F);
    x_new_standardized = scale(x_new,center=obs_mean_x,scale=obs_sd_x);
    x_new_centered = scale(x_new,center=obs_mean_x,scale=F);
    matrix_y_new = matrix(y_new,nrow=mc_total_samps,ncol=n_new,byrow=T);
    
    ## Separation statistics ====##########################################################################
    if(calc_n_piv) {
      phil_approx = fast_calculate_sep_stat(x_standardized = x_standardized,
                                            y = y,
                                            which_stat = "piv",
                                            maxk = algorithm2_maxk,
                                            ndraws = algorithm2_ndraws);
      algorithm2_results[i,"n_piv"] = phil_approx$statistic;
      rm(phil_approx);
    }
    if(calc_n_comp) {
      phil_approx = fast_calculate_sep_stat(x_standardized = x_standardized,
                                            y = y,
                                            which_stat = "comp",
                                            maxk = algorithm2_maxk,
                                            ndraws = algorithm2_ndraws);
      algorithm2_results[i,"n_comp"] = min(phil_approx$statistic,algorithm2_results[i,"n_piv"],na.rm=T);
      rm(phil_approx);
    }
    if(calc_n_over) {
      phil_approx = fast_calculate_sep_stat(x_standardized = x_standardized,
                                            y = y,
                                            which_stat = "over",
                                            maxk = algorithm2_maxk,
                                            ndraws = algorithm2_ndraws);
      algorithm2_results[i,"n_over"] = min(phil_approx$statistic,algorithm2_results[i,"n_comp"],na.rm=T);
      rm(phil_approx);
    }
    
    
    if(fit_methods) {
      
      ## Compile templates ====##########################################################################
      if(!stan_compiled) {
        begin_compile = Sys.time();
        
        assign("hsbeta_template",glm_HSBeta(stan_path = paste0(stan_file_path,hsbeta_file_name)));
        assign("logisbeta_template",glm_LogisBeta(stan_path = paste0(stan_file_path,logisbeta_file_name)));
        
        end_compile = Sys.time();  
        stan_compiled = T;
      } 
      only_prior = F;
      
      ## t3(10)&HS ====##########################################################################
      curr_prior = "StudT10";
      if(curr_prior%in%intercept_prior_names) {
        #HS prior on Betas
        curr_method = paste0("HS_",curr_prior);
        alpha_prior_args = list(type = 0,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                mean = 0, 
                                scale = 10, 
                                power = 0,#ignored unless alpha_prior_type = 2;
                                dof = 3#ignored unless alpha_prior_type = 0;
        );
        foo = glm_HSBeta(stan_path = paste0(stan_file_path,hsbeta_file_name), 
                         stan_fit = hsbeta_template,
                         y = y, 
                         x_standardized = x_standardized, 
                         beta_global_scale = beta_global_scale,
                         local_dof = local_dof, 
                         global_dof = global_dof, 
                         slab_precision = slab_precision,
                         alpha_prior_args = alpha_prior_args,
                         only_prior = only_prior, 
                         mc_warmup = mc_warmup, 
                         mc_iter_after_warmup = mc_iter_after_warmup, 
                         mc_chains = mc_chains, 
                         mc_thin = mc_thin, 
                         mc_stepsize = mc_stepsize, 
                         mc_adapt_delta = mc_adapt_delta,
                         mc_max_treedepth = mc_max_treedepth,
                         ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        standardized_beta_samps = foo$curr_beta;
        standardized_alpha_samps = foo$curr_alpha;
        standardized_theta_samps = foo$curr_theta;
        
        store_mean_alpha[i,curr_method] = mean(standardized_alpha_samps);
        
        #Estimation
        mean_standardized_beta = colMeans(standardized_beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_standardized_beta/obs_sd_x;
        store_sd_beta[[curr_method]][i,] =  (apply(standardized_beta_samps,2,sd)/obs_sd_x);
        beta_low = round(apply(standardized_beta_samps,2,quantile,.25)/obs_sd_x,3);
        beta_high = round(apply(standardized_beta_samps,2,quantile,.75)/obs_sd_x,3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((standardized_beta_samps/obs_sd_x - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(standardized_alpha_samps + standardized_beta_samps%*%t(x_new_standardized));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new == 1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new == 0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_standardized_beta%*%t(x_new_standardized))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = standardized_alpha_samps;
          store_beta[[curr_method]] = (standardized_beta_samps/matrix_obs_sd_x);
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,standardized_beta_samps,standardized_alpha_samps,mean_standardized_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        if(!do_logisbeta_prior){rm(alpha_prior_args)};
      }
      
      ## t3(10)&Logistic(1) ====##########################################################################
      if(do_logisbeta_prior && curr_prior%in%intercept_prior_names) {
        #Logistic(1) prior on Betas
        curr_method = paste0("Logis_",curr_prior);
        foo = glm_LogisBeta(stan_path = paste0(stan_file_path,logisbeta_stan_filename), 
                            stan_fit = logisbeta_template,
                            y = y, 
                            x_centered = x_centered,
                            beta_expit_shape = 1,
                            alpha_prior_args = alpha_prior_args,
                            only_prior = only_prior, 
                            mc_warmup = mc_warmup, 
                            mc_iter_after_warmup = mc_iter_after_warmup, 
                            mc_chains = mc_chains, 
                            mc_thin = mc_thin, 
                            mc_stepsize = mc_stepsize, 
                            mc_adapt_delta = mc_adapt_delta,
                            mc_max_treedepth = mc_max_treedepth,
                            ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        #Centering the covariates only impacts the intercept
        #Scaling the covariates only impacts the regression coefficients (beta)
        #The Logistic(1) prior centers but does not scale the covariate
        beta_samps = foo$curr_beta;
        centered_alpha_samps = foo$curr_alpha;
        
        store_mean_alpha[i,curr_method] = mean(centered_alpha_samps);
        
        #Estimation
        mean_beta = colMeans(beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_beta;
        store_sd_beta[[curr_method]][i,] =  apply(beta_samps,2,sd);
        beta_low = round(apply(beta_samps,2,quantile,.25),3);
        beta_high = round(apply(beta_samps,2,quantile,.75),3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((beta_samps - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(centered_alpha_samps + beta_samps%*%t(x_new_centered));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_beta%*%t(x_new_centered))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = centered_alpha_samps;
          store_beta[[curr_method]] = beta_samps;
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,curr_prior,beta_samps,centered_alpha_samps,mean_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        rm(alpha_prior_args);
      }
      
      ## N(10)&HS ====##########################################################################
      curr_prior = "Normal10";
      if(curr_prior%in%intercept_prior_names) {
        #HS prior on Betas
        curr_method = paste0("HS_",curr_prior);
        alpha_prior_args = list(type = 2,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                mean = 0, 
                                scale = 10, 
                                power = 2,#ignored unless alpha_prior_type = 2;
                                dof = 3#ignored unless alpha_prior_type = 0;
        );
        foo = glm_HSBeta(stan_path = paste0(stan_file_path,hsbeta_file_name), 
                         stan_fit = hsbeta_template,
                         y = y, 
                         x_standardized = x_standardized, 
                         beta_global_scale = beta_global_scale,
                         local_dof = local_dof, 
                         global_dof = global_dof, 
                         slab_precision = slab_precision,
                         alpha_prior_args = alpha_prior_args,
                         only_prior = only_prior, 
                         mc_warmup = mc_warmup, 
                         mc_iter_after_warmup = mc_iter_after_warmup, 
                         mc_chains = mc_chains, 
                         mc_thin = mc_thin, 
                         mc_stepsize = mc_stepsize, 
                         mc_adapt_delta = mc_adapt_delta,
                         mc_max_treedepth = mc_max_treedepth,
                         ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        standardized_beta_samps = foo$curr_beta;
        standardized_alpha_samps = foo$curr_alpha;
        standardized_theta_samps = foo$curr_theta;
        
        store_mean_alpha[i,curr_method] = mean(standardized_alpha_samps);
        
        #Estimation
        mean_standardized_beta = colMeans(standardized_beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_standardized_beta/obs_sd_x;
        store_sd_beta[[curr_method]][i,] =  (apply(standardized_beta_samps,2,sd)/obs_sd_x);
        beta_low = round(apply(standardized_beta_samps,2,quantile,.25)/obs_sd_x,3);
        beta_high = round(apply(standardized_beta_samps,2,quantile,.75)/obs_sd_x,3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((standardized_beta_samps/obs_sd_x - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(standardized_alpha_samps + standardized_beta_samps%*%t(x_new_standardized));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_standardized_beta%*%t(x_new_standardized))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = standardized_alpha_samps;
          store_beta[[curr_method]] = (standardized_beta_samps/matrix_obs_sd_x);
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,standardized_beta_samps,standardized_alpha_samps,mean_standardized_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        if(!do_logisbeta_prior){rm(alpha_prior_args)};
      }
      
      ## N(10)&Logistic(1) ====##########################################################################
      if(do_logisbeta_prior && curr_prior%in%intercept_prior_names) {
        #Logistic(1) prior on Betas
        curr_method = paste0("Logis_",curr_prior);
        foo = glm_LogisBeta(stan_path = paste0(stan_file_path,logisbeta_stan_filename), 
                            stan_fit = logisbeta_template,
                            y = y, 
                            x_centered = x_centered,
                            beta_expit_shape = 1,
                            alpha_prior_args = alpha_prior_args,
                            only_prior = only_prior, 
                            mc_warmup = mc_warmup, 
                            mc_iter_after_warmup = mc_iter_after_warmup, 
                            mc_chains = mc_chains, 
                            mc_thin = mc_thin, 
                            mc_stepsize = mc_stepsize, 
                            mc_adapt_delta = mc_adapt_delta,
                            mc_max_treedepth = mc_max_treedepth,
                            ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        #Centering the covariates only impacts the intercept
        #Scaling the covariates only impacts the regression coefficients (beta)
        #The Logistic(1) prior centers but does not scale the covariate
        beta_samps = foo$curr_beta;
        centered_alpha_samps = foo$curr_alpha;
        
        store_mean_alpha[i,curr_method] = mean(centered_alpha_samps);
        
        #Estimation
        mean_beta = colMeans(beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_beta;
        store_sd_beta[[curr_method]][i,] =  apply(beta_samps,2,sd);
        beta_low = round(apply(beta_samps,2,quantile,.25),3);
        beta_high = round(apply(beta_samps,2,quantile,.75),3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((beta_samps - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(centered_alpha_samps + beta_samps%*%t(x_new_centered));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_beta%*%t(x_new_centered))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = centered_alpha_samps;
          store_beta[[curr_method]] = beta_samps;
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,curr_prior,beta_samps,centered_alpha_samps,mean_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        rm(alpha_prior_args);
      }
      
      ## EP&HS ====##########################################################################
      for(k in 1:length(exppow_prior_args)) {
        curr_prior = names(exppow_prior_args)[k];
        if(curr_prior%in%intercept_prior_names) {
          #HS prior on Betas
          curr_method = paste0("HS_",curr_prior);
          alpha_prior_args = list(type = 2,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                  mean = 0, 
                                  scale = exppow_prior_args[[k]]$alpha_scale,
                                  power = exppow_prior_args[[k]]$alpha_power,#ignored unless alpha_prior_type = 2;
                                  dof = 3#ignored unless alpha_prior_type = 0;
          );
          foo = glm_HSBeta(stan_path = paste0(stan_file_path,hsbeta_file_name), 
                           stan_fit = hsbeta_template,
                           y = y, 
                           x_standardized = x_standardized, 
                           beta_global_scale = beta_global_scale,
                           local_dof = local_dof, 
                           global_dof = global_dof, 
                           slab_precision = slab_precision,
                           alpha_prior_args = alpha_prior_args,
                           only_prior = only_prior, 
                           mc_warmup = mc_warmup, 
                           mc_iter_after_warmup = mc_iter_after_warmup, 
                           mc_chains = mc_chains, 
                           mc_thin = mc_thin, 
                           mc_stepsize = mc_stepsize, 
                           mc_adapt_delta = mc_adapt_delta,
                           mc_max_treedepth = mc_max_treedepth,
                           ntries = ntries_per_iter);
          
          store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
          max_divergences_by_prior[i,curr_method] = foo$max_divergences;
          accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
          max_rhat_by_prior[i,curr_method] = foo$max_rhat;
          
          standardized_beta_samps = foo$curr_beta;
          standardized_alpha_samps = foo$curr_alpha;
          standardized_theta_samps = foo$curr_theta;
          
          store_mean_alpha[i,curr_method] = mean(standardized_alpha_samps);
          
          #Estimation
          mean_standardized_beta = colMeans(standardized_beta_samps);
          store_mean_beta[[curr_method]][i,] = mean_standardized_beta/obs_sd_x;
          store_sd_beta[[curr_method]][i,] =  (apply(standardized_beta_samps,2,sd)/obs_sd_x);
          beta_low = round(apply(standardized_beta_samps,2,quantile,.25)/obs_sd_x,3);
          beta_high = round(apply(standardized_beta_samps,2,quantile,.75)/obs_sd_x,3);
          store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
          store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
          bayes_rmse[i,curr_method] = sqrt(mean(rowSums((standardized_beta_samps/obs_sd_x - matrix_true_betas)^2)));
          
          #Predictive density
          x_new_prob = expit(standardized_alpha_samps + standardized_beta_samps%*%t(x_new_standardized));
          x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
          elpd[i,curr_method] = sum(x_new_lpd);
          elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
          elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
          
          #AUC
          if(diff(range(y_new)) > .Machine$double.eps^0.5) {
            store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_standardized_beta%*%t(x_new_standardized))))
          }
          
          if(verbose && data_analysis && (niter == 1)) {
            store_alpha[[curr_method]] = standardized_alpha_samps;
            store_beta[[curr_method]] = (standardized_beta_samps/matrix_obs_sd_x);
            store_fitted_probs[[curr_method]] = x_new_prob;
          }
          rm(foo,curr_method,standardized_beta_samps,standardized_alpha_samps,mean_standardized_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
          rm(alpha_prior_args);
        }
      }
      
      ## EP&Logistic(1) ====##########################################################################
      for(k in 1:length(exppow_prior_args)) {
        curr_prior = names(exppow_prior_args)[k];
        if(do_logisbeta_prior && curr_prior%in%intercept_prior_names) {
          #Logistic(1) prior on Betas
          curr_method = paste0("Logis_",curr_prior);
          alpha_prior_args = list(type = 2,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                  mean = 0, 
                                  scale = exppow_prior_args[[k]]$alpha_scale,
                                  power = exppow_prior_args[[k]]$alpha_power,#ignored unless alpha_prior_type = 2;
                                  dof = 3#ignored unless alpha_prior_type = 0;
          );
          foo = glm_LogisBeta(stan_path = paste0(stan_file_path,logisbeta_stan_filename), 
                              stan_fit = logisbeta_template,
                              y = y, 
                              x_centered = x_centered,
                              beta_expit_shape = 1,
                              alpha_prior_args = alpha_prior_args,
                              only_prior = only_prior, 
                              mc_warmup = mc_warmup, 
                              mc_iter_after_warmup = mc_iter_after_warmup, 
                              mc_chains = mc_chains, 
                              mc_thin = mc_thin, 
                              mc_stepsize = mc_stepsize, 
                              mc_adapt_delta = mc_adapt_delta,
                              mc_max_treedepth = mc_max_treedepth,
                              ntries = ntries_per_iter);
          
          store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
          max_divergences_by_prior[i,curr_method] = foo$max_divergences;
          accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
          max_rhat_by_prior[i,curr_method] = foo$max_rhat;
          
          #Centering the covariates only impacts the intercept
          #Scaling the covariates only impacts the regression coefficients (beta)
          #The Logistic(1) prior centers but does not scale the covariate
          beta_samps = foo$curr_beta;
          centered_alpha_samps = foo$curr_alpha;
          
          store_mean_alpha[i,curr_method] = mean(centered_alpha_samps);
          
          #Estimation
          mean_beta = colMeans(beta_samps);
          store_mean_beta[[curr_method]][i,] = mean_beta;
          store_sd_beta[[curr_method]][i,] =  apply(beta_samps,2,sd);
          beta_low = round(apply(beta_samps,2,quantile,.25),3);
          beta_high = round(apply(beta_samps,2,quantile,.75),3);
          store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
          store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
          bayes_rmse[i,curr_method] = sqrt(mean(rowSums((beta_samps - matrix_true_betas)^2)));
          
          #Predictive density
          x_new_prob = expit(centered_alpha_samps + beta_samps%*%t(x_new_centered));
          x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
          elpd[i,curr_method] = sum(x_new_lpd);
          elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
          elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
          
          #AUC
          if(diff(range(y_new)) > .Machine$double.eps^0.5) {
            store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_beta%*%t(x_new_centered))))
          }
          
          if(verbose && data_analysis && (niter == 1)) {
            store_alpha[[curr_method]] = centered_alpha_samps;
            store_beta[[curr_method]] = beta_samps;
            store_fitted_probs[[curr_method]] = x_new_prob;
          }
          rm(foo,curr_method,curr_prior,beta_samps,centered_alpha_samps,mean_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
          rm(alpha_prior_args);
        }
      }
      
      ## Logis&HS ====##########################################################################
      curr_prior = "LogisticAdaptive";
      if(curr_prior%in%intercept_prior_names) {
        #HS prior on Betas
        curr_method = paste0("HS_",curr_prior);
        alpha_prior_args = list(type = 1,#0 = student-t; 1 = logistic; 2 = exponential power, inc. normal
                                mean = 0, 
                                scale = logistic_adaptive_alpha_scale,
                                power = 1,#ignored unless alpha_prior_type = 2;
                                dof = 3#ignored unless alpha_prior_type = 0;
        );
        foo = glm_HSBeta(stan_path = paste0(stan_file_path,hsbeta_file_name), 
                         stan_fit = hsbeta_template,
                         y = y, 
                         x_standardized = x_standardized, 
                         beta_global_scale = beta_global_scale,
                         local_dof = local_dof, 
                         global_dof = global_dof, 
                         slab_precision = slab_precision,
                         alpha_prior_args = alpha_prior_args,
                         only_prior = only_prior, 
                         mc_warmup = mc_warmup, 
                         mc_iter_after_warmup = mc_iter_after_warmup, 
                         mc_chains = mc_chains, 
                         mc_thin = mc_thin, 
                         mc_stepsize = mc_stepsize, 
                         mc_adapt_delta = mc_adapt_delta,
                         mc_max_treedepth = mc_max_treedepth,
                         ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        standardized_beta_samps = foo$curr_beta;
        standardized_alpha_samps = foo$curr_alpha;
        standardized_theta_samps = foo$curr_theta;
        
        store_mean_alpha[i,curr_method] = mean(standardized_alpha_samps);
        
        #Estimation
        mean_standardized_beta = colMeans(standardized_beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_standardized_beta/obs_sd_x;
        store_sd_beta[[curr_method]][i,] =  (apply(standardized_beta_samps,2,sd)/obs_sd_x);
        beta_low = round(apply(standardized_beta_samps,2,quantile,.25)/obs_sd_x,3);
        beta_high = round(apply(standardized_beta_samps,2,quantile,.75)/obs_sd_x,3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((standardized_beta_samps/obs_sd_x - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(standardized_alpha_samps + standardized_beta_samps%*%t(x_new_standardized));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_standardized_beta%*%t(x_new_standardized))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = standardized_alpha_samps;
          store_beta[[curr_method]] = (standardized_beta_samps/matrix_obs_sd_x);
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,standardized_beta_samps,standardized_alpha_samps,mean_standardized_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        if(!do_logisbeta_prior){rm(alpha_prior_args)};
      }
      
      ## Logis&Logistic(1) ====##########################################################################
      if(do_logisbeta_prior && curr_prior%in%intercept_prior_names) {
        #Logistic(1) prior on Betas
        curr_method = paste0("Logis_",curr_prior);
        foo = glm_LogisBeta(stan_path = paste0(stan_file_path,logisbeta_file_name), 
                            stan_fit = logisbeta_template,
                            y = y, 
                            x_centered = x_centered,
                            beta_expit_shape = 1,
                            alpha_prior_args = alpha_prior_args,
                            only_prior = only_prior, 
                            mc_warmup = mc_warmup, 
                            mc_iter_after_warmup = mc_iter_after_warmup, 
                            mc_chains = mc_chains, 
                            mc_thin = mc_thin, 
                            mc_stepsize = mc_stepsize, 
                            mc_adapt_delta = mc_adapt_delta,
                            mc_max_treedepth = mc_max_treedepth,
                            ntries = ntries_per_iter);
        
        store_attributes_data[i,"divergence"] = pmax(store_attributes_data[i,"divergence"],foo$max_divergences);
        max_divergences_by_prior[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_prior[i,curr_method] = foo$accepted_divergences;
        max_rhat_by_prior[i,curr_method] = foo$max_rhat;
        
        #Centering the covariates only impacts the intercept
        #Scaling the covariates only impacts the regression coefficients (beta)
        #The Logistic(1) prior centers but does not scale the covariate
        beta_samps = foo$curr_beta;
        centered_alpha_samps = foo$curr_alpha;
        
        store_mean_alpha[i,curr_method] = mean(centered_alpha_samps);
        
        #Estimation
        mean_beta = colMeans(beta_samps);
        store_mean_beta[[curr_method]][i,] = mean_beta;
        store_sd_beta[[curr_method]][i,] =  apply(beta_samps,2,sd);
        beta_low = round(apply(beta_samps,2,quantile,.25),3);
        beta_high = round(apply(beta_samps,2,quantile,.75),3);
        store_50CIcoverage_beta[[curr_method]][i,] = (beta_low < true_betas) & (beta_high > true_betas);
        store_50CIwidth_beta[[curr_method]][i,] = beta_high - beta_low;
        bayes_rmse[i,curr_method] = sqrt(mean(rowSums((beta_samps - matrix_true_betas)^2)));
        
        #Predictive density
        x_new_prob = expit(centered_alpha_samps + beta_samps%*%t(x_new_centered));
        x_new_lpd = log(colMeans((x_new_prob ^ matrix_y_new)*((1-x_new_prob) ^ (1-matrix_y_new))));
        elpd[i,curr_method] = sum(x_new_lpd);
        elpd_events[i,curr_method] = sum(x_new_lpd[which(y_new==1)]);
        elpd_nonevents[i,curr_method] = sum(x_new_lpd[which(y_new==0)]);
        
        #AUC
        if(diff(range(y_new)) > .Machine$double.eps^0.5) {
          store_auc[i,curr_method] = pROC::auc(pROC::roc(y_new, drop(mean_beta%*%t(x_new_centered))))
        }
        
        if(verbose && data_analysis && (niter == 1)) {
          store_alpha[[curr_method]] = centered_alpha_samps;
          store_beta[[curr_method]] = beta_samps;
          store_fitted_probs[[curr_method]] = x_new_prob;
        }
        rm(foo,curr_method,curr_prior,beta_samps,centered_alpha_samps,mean_beta,beta_low,beta_high,x_new_prob,x_new_lpd);
        rm(alpha_prior_args);
      }
      
      rm(x,x_centered,x_standardized,x_new,x_new_centered,x_new_standardized,
         y,y_new,matrix_y_new,
         matrix_obs_sd_x,matrix_obs_mean_x,obs_sd_x,obs_mean_x);
      if(!data_analysis) {rm(curr_true_probs);}
      gc();
    }
  }
  
  standardized_true_betas = sd_actual_x * true_betas;
  standardized_true_alpha = true_alpha + sum(true_betas*mu_actual_x);
  
  ##Average error over all iterations except those with unresolved divergences
  AllErrBeta =
    AllErrBeta_NotStndzd = matrix(NA,nrow=niter,ncol=length(method_names),dimnames = list(NULL,method_names));
  matrix_true_betas = matrix(true_betas,nrow=niter,ncol=p,byrow=T);
  for(k in method_names) {
    if(p>1) {
      #beta_specific_err is a matrix, with rows corresponding to different simulations
      #each row is comprised of the diagonal elements of (beta_est-beta)%*%t(beta_est-beta)%*%avg_cov_obs_x
      #the rowsum is of the elements is the mse of beta for that simulation (for that method)
      beta_specific_err = t(apply(store_mean_beta[[k]]-matrix_true_betas,1,"%*%",avg_cov_obs_x))*(store_mean_beta[[k]]-matrix_true_betas);
      AllErrBeta[,k] = sqrt(rowSums(beta_specific_err));
      AllErrBeta_NotStndzd[,k] = sqrt(rowSums((store_mean_beta[[k]]-matrix_true_betas)^2));
    } else {
      beta_specific_err = (store_mean_beta[[k]]-matrix_true_betas)^2*drop(avg_cov_obs_x)
      AllErrBeta[,k] = sqrt(beta_specific_err);
      AllErrBeta_NotStndzd[,k] = sqrt((store_mean_beta[[k]]-matrix_true_betas)^2);
    }
  }
  exclude_list = list(c("const_y","divergence"),c("const_y","const_x","copy_x","divergence"));
  MedianErrBeta = 
    MedianEffBeta = 
    MedianELPD = 
    MedianELPD_events = 
    MedianELPD_nonevents = 
    vector("list",length(exclude_list));
  names(MedianErrBeta) = 
    names(MedianEffBeta) = 
    names(MedianELPD) = 
    names(MedianELPD_events) = 
    names(MedianELPD_nonevents) = 
    list("ExcludeDivergentConstY","ExcludeAllProblems");
  for(j in 1:length(exclude_list)) {
    ##Average error over iterations without duplicated or constant columns
    curr_index = which(rowSums(store_attributes_data[,exclude_list[[j]],drop = F], na.rm = T) == 0);
    if(length(curr_index)) {
      matrix_true_betas = matrix(true_betas,nrow=length(curr_index),ncol=p,byrow=T);
      MedianErrBeta[[j]] = 
        MedianEffBeta[[j]] = matrix(NA,nrow=length(method_names),ncol=3,dimnames=list(method_names,c("all","|b|>0","|b|=0")));
      MedianELPD[[j]] = 
        MedianELPD_events[[j]] =  
        MedianELPD_nonevents[[j]] = matrix(NA,nrow=1,ncol=length(method_names),dimnames = list(NULL,method_names));
      for(k in method_names) {
        MedianErrBeta[[j]][k,"all"] = median(AllErrBeta[curr_index,k]);
        MedianEffBeta[[j]][k,"all"] = median(rowSums(store_50CIwidth_beta[[k]][curr_index,,drop = F]));
        MedianELPD[[j]][1,k] = median(elpd[curr_index,k]);
        MedianELPD_events[[j]][1,k] = median(elpd_events[curr_index,k]);
        MedianELPD_nonevents[[j]][1,k] = median(elpd_nonevents[curr_index,k]);
        
        if(any(abs(standardized_true_betas)>0)) {
          foo = which(abs(standardized_true_betas)>0);
          MedianErrBeta[[j]][k,"|b|>0"] = median(sqrt(rowSums((store_mean_beta[[k]][curr_index,,drop = F]-matrix_true_betas)[,foo,drop = F]^2*(sd_actual_x[foo]^2))));
          MedianEffBeta[[j]][k,"|b|>0"] = median(rowSums(store_50CIwidth_beta[[k]][curr_index,foo,drop = F]));
          rm(foo);
        }
        if(any(abs(standardized_true_betas)==0)) {
          foo = which(abs(standardized_true_betas)==0);
          MedianErrBeta[[j]][k,"|b|=0"] = median(sqrt(rowSums((store_mean_beta[[k]][curr_index,,drop = F]-matrix_true_betas)[,foo,drop = F]^2*(sd_actual_x[foo]^2))));
          MedianEffBeta[[j]][k,"|b|=0"] = median(rowSums(store_50CIwidth_beta[[k]][curr_index,foo,drop = F]));
          rm(foo);
        }
      }
    } 
  }
  
  generating_params = list(n = n,
                           n_new = n_new,
                           true_betas = true_betas,
                           standardized_true_betas = standardized_true_betas,
                           true_alpha = true_alpha,
                           standardized_true_alpha = standardized_true_alpha,
                           x_binom = x_binom,
                           case_control_ratio = case_control_ratio,
                           mu_latent_x = mu_latent_x, 
                           mu_actual_x = mu_actual_x, 
                           chol_latent_x = chol_latent_x, 
                           avg_cov_obs_x = avg_cov_obs_x,
                           store_obs_mean_x = store_obs_mean_x,
                           store_obs_sd_x = store_obs_sd_x,
                           population_prevalence = population_prevalence);
  model_params = list(hiershrink_prior_num_relevant = hiershrink_prior_num_relevant,
                      local_dof = local_dof,
                      global_dof = global_dof,
                      beta_global_scale = beta_global_scale,
                      exppow_prior_args = exppow_prior_args,
                      logistic_adaptive_alpha_scale = logistic_adaptive_alpha_scale);
  data_params = list(x_all = x_all, 
                     y_all = y_all,
                     frac_training = frac_training,
                     num_partitions = num_partitions,
                     prespecified_training_subsets = prespecified_training_subsets,
                     store_obs_mean_x = store_obs_mean_x,
                     store_obs_sd_x = store_obs_sd_x)
  sim_params = list(array_id = array_id,
                    scenario = scenario,
                    niter = niter,
                    random_seed = random_seed,
                    data_seeds = data_seeds);
  mc_params = list(mc_warmup = mc_warmup,
                   mc_iter_after_warmup = mc_iter_after_warmup,
                   mc_chains = mc_chains,
                   mc_thin = mc_thin,
                   mc_adapt_delta = mc_adapt_delta,
                   mc_stepsize = mc_stepsize,
                   mc_max_treedepth = mc_max_treedepth,
                   mc_ntries = ntries_per_iter);
  algorithm2_params = list(algorithm2_ndraws = algorithm2_ndraws,
                           algorithm2_maxk = algorithm2_maxk);
  
  if(data_analysis) {
    return(list(runtime = difftime(Sys.time(),begin,units="hours"),
                data_params = data_params,
                model_params = model_params,
                mc_params = mc_params,
                hiershrink_prior_num_relevant = generating_params$hiershrink_prior_num_relevant,
                algorithm2_params = algorithm2_params,
                store_mean_alpha = store_mean_alpha,
                store_mean_beta = store_mean_beta,
                store_alpha = store_alpha,
                store_beta = store_beta,
                store_fitted_probs = store_fitted_probs,
                store_attributes_data = store_attributes_data,
                accepted_divergences_by_prior = accepted_divergences_by_prior,
                max_divergences_by_prior = max_divergences_by_prior,
                algorithm2_results = algorithm2_results,
                MedianEffBeta = MedianEffBeta,
                AllELPD = elpd,
                AllELPD_events = elpd_events,
                AllELPD_nonevents = elpd_nonevents,
                MedianELPD = MedianELPD,
                MedianELPD_events = MedianELPD_events,
                MedianELPD_nonevents = MedianELPD_nonevents,
                store_auc = store_auc
    ));
  } else {
    return(list(runtime = difftime(Sys.time(),begin,units="hours"),
                model_params = model_params,
                generating_params = generating_params,
                sim_params = sim_params,
                mc_params = mc_params,
                algorithm2_params = algorithm2_params,
                store_attributes_data = store_attributes_data,
                accepted_divergences_by_prior = accepted_divergences_by_prior,
                max_divergences_by_prior = max_divergences_by_prior,
                interesting_datasets = interesting_datasets,
                store_mean_alpha = store_mean_alpha,
                store_mean_beta = store_mean_beta,
                store_sd_beta = store_sd_beta,
                store_50CIcoverage_beta = store_50CIcoverage_beta,
                store_50CIwidth_beta = store_50CIwidth_beta,
                algorithm2_results = algorithm2_results,
                AllErrBeta = AllErrBeta,
                AllErrBeta_NotStndzd = AllErrBeta_NotStndzd,
                MedianErrBeta = MedianErrBeta,
                MedianEffBeta = MedianEffBeta,
                AllELPD = elpd,
                AllELPD_events = elpd_events,
                AllELPD_nonevents = elpd_nonevents,
                MedianELPD = MedianELPD,
                MedianELPD_events = MedianELPD_events,
                MedianELPD_nonevents = MedianELPD_nonevents,
                store_auc = store_auc,
                bayes_rmse = bayes_rmse
    ));              
  }  
  
}

fast_calculate_sep_stat = function(x_standardized, y, which_stat = "comp", maxk = 10, ndraws = 1e3, seed = .Machine$integer.max) {
  set.seed(seed);
  if(which_stat == "comp") {
    #calculates number of observations needed to remove to induce complete separation
    overlap_allowed = F; 
    no_intercept = F; 
  } else if(which_stat == "piv") {
    #calculates number of observations needed to remove to induce pivotal separation
    overlap_allowed = F; 
    no_intercept = T; 
  } else if(which_stat == "over") {
    #calculates number of observations needed to remove to induce quasi-complete separation
    overlap_allowed = T; 
    no_intercept = F; 
  } else {
    stop("'which_stat' must be one of 'comp', 'piv', or 'over'");
  }
  
  n = length(y);
  p = ifelse(class(x_standardized) == "numeric", 1, ncol(x_standardized));
  separating_plane = numeric(p+1);
  min_ub = n;
  lowest_ub_found = F;
  solitary_separators = NULL;
  
  if(no_intercept) {
    x_standardized_aug = x_standardized;
  } else {
    x_standardized_aug = cbind(1,x_standardized);
  }
  
  if(diff(range(y)) < .Machine$double.eps^0.5) {
    warning("No variation in 'y' means that complete separation is trivially satisfied"); 
    return_vals = list(statistic = NA,
                       included_subset = NA,
                       linear_separator = NA,
                       solitary_separators = NA);
  } else if(p == 1) {#If single x caused separation, then, by definition, it was solitary
    model_uni = glm.fit(x_standardized_aug, y, family = binomial());
    fitted_uni = model_uni$linear.predictors;
    #THIS LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_uni, FUN=sum);
    foo = data.frame(y = y,fitted_uni = fitted_uni,const = 1);
    aggregate_y = foo$y;
    aggregate_fitted_uni = foo$fitted_uni;
    aggregate_count = foo[,3];
    order_increasing_y = order(aggregate_fitted_uni,aggregate_y);
    order_decreasing_y = order(aggregate_fitted_uni,-aggregate_y);
    aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_uni[order_increasing_y]));
    ncomplete_uni1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
    ncomplete_uni2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
    min_ub = min(ncomplete_uni1$ub,ncomplete_uni2$ub);
    min_plane = model_uni$coefficients;
    if(ncomplete_uni1$ub <= ncomplete_uni2$ub) {
      min_subset = sort(order_increasing_y[ncomplete_uni1$sep_subseq]);
    } else {
      min_subset = sort(order_decreasing_y[ncomplete_uni2$sep_subseq]);
    }
    return_vals = list(statistic = min_ub,
                       included_subset = min_subset,
                       linear_separator = min_plane,
                       solitary_separators = 1 * (min_ub == 0));
  } else {#Otherwise assume that x_standardized is a matrix and go through each column at a time
    if(no_intercept) {
      min_plane = numeric(p);
      names(min_plane) = paste0("x_standardized",colnames(x_standardized));
    } else {
      min_plane = numeric(p+1);
      names(min_plane) = c("(Intercept)",paste0("x_standardized",colnames(x_standardized)));
    }
    for(i in 1:p) {
      model_uni = glm.fit(x_standardized_aug[,c(!no_intercept,i+!no_intercept),drop = F], y, family = binomial());
      fitted_uni = model_uni$linear.predictors;
      if(diff(range(fitted_uni)) < .Machine$double.eps^0.5) {next;}
      #THE FOLLOWING LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_uni, FUN=sum);
      foo = data.frame(y = y,fitted_uni = fitted_uni,const = 1);
      aggregate_y = foo$y;
      aggregate_fitted_uni = foo$fitted_uni;
      aggregate_count = foo[,3];
      order_increasing_y = order(aggregate_fitted_uni,aggregate_y);
      order_decreasing_y = order(aggregate_fitted_uni,-aggregate_y);
      aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_uni[order_increasing_y]));
      ncomplete_uni1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
      ncomplete_uni2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
      curr_min_ub = min(ncomplete_uni1$ub,ncomplete_uni2$ub);
      if(curr_min_ub < min_ub) {
        min_ub = curr_min_ub;
        min_plane = min_plane * 0;
        #cat("new min_ub =", curr_min_ub,"\n");
        if(no_intercept) {
          min_plane[i] = model_uni$coefficients;
        } else {
          min_plane[c(1, 1 + i)] = model_uni$coefficients;
        }
        if(ncomplete_uni1$ub <= ncomplete_uni2$ub) {
          min_subset = sort(order_increasing_y[ncomplete_uni1$sep_subseq]);
        } else {
          min_subset = sort(order_decreasing_y[ncomplete_uni2$sep_subseq]);
        }
      }
      if(curr_min_ub == 0) {
        solitary_separators = c(solitary_separators,i);
      }
    }
    rm(foo,aggregate_y,aggregate_fitted_uni,aggregate_count,ncomplete_uni1,ncomplete_uni2,curr_min_ub)
    
    loop_stop = max(maxk,2);
    #If the smallest upper bound is 0, then this is trivially the lowest upperbound.
    lowest_ub_found = (min_ub == 0);
    #Now if needed go to iterative sampling of subset
    if(!lowest_ub_found & (loop_stop > 1)) {
      combined_data = cbind(x_standardized,y);
      which_unique = which(!duplicated(combined_data));
      fitted_heldout = numeric(length(which_unique));
      affiliated_obs = vector("list",length(which_unique));
      names(affiliated_obs) = names(fitted_heldout) = which_unique;
      for(j in which_unique) {
        curr_remove = as.numeric(which(rowSums(combined_data == rep(combined_data[j,], each = nrow(combined_data))) == ncol(combined_data)));
        affiliated_obs[[as.character(j)]] = curr_remove;
        curr_keep = setdiff(1:n,curr_remove);
        model_subset = glm.fit(x_standardized_aug[curr_keep,,drop = F], y[curr_keep], family = binomial());
        fitted_heldout[as.character(j)] = sum(x_standardized_aug[curr_remove[1],]*ifelse(is.na(model_subset$coefficients),0,model_subset$coefficients));
      }
      model_full = glm.fit(x_standardized_aug, y, family = binomial());
      fitted_full = model_full$linear.predictors[which_unique];
      weights_y = 1 / abs(fitted_heldout-fitted_full)^0.5;
      weights_y[y[which_unique]==1] = pmin(weights_y[y[which_unique]==1], sum(1-y) + 1e-3);
      weights_y[y[which_unique]==0] = pmin(weights_y[y[which_unique]==0], sum(y) + 1e-3);
      rm(model_full,fitted_full,fitted_heldout,model_subset,curr_keep,curr_remove);
      
      j = -1;
      n_unique = length(which_unique);
      n_most_likely = round(ndraws/2);#how many of the approximately most likely subsets to look at
      n_random_weighted = pmax(1, ndraws - n_most_likely);#how many to randomly sample based upon less uniform weights
      while(T) {
        j = j + 1;
        choose_n_nminusj = choose(n_unique, n_unique-j);
        ten_times_n_most_likely = pmin(choose_n_nminusj, 10 * n_most_likely);
        #cat(j,min_ub,"\n");
        if(choose_n_nminusj < ndraws) {#all possible combinations is feasible
          all_subsets = combn(n_unique,n_unique-j);
        } else {
          #combination of most likely weighted subsets
          #and weighted random sampling from among remaining subsets 
          if(choose_n_nminusj > ten_times_n_most_likely) {
            k = max((0:(n_unique-j))[mapply(FUN = "choose",n=n_unique-(0:(n_unique-j)),k=(n_unique-j)-(0:(n_unique-j))) >= ten_times_n_most_likely]);
            if(k > 0) {
              all_subsets = rbind(matrix(rep(1:k,times=ten_times_n_most_likely),nrow=k,ncol=ten_times_n_most_likely),k + combn(n_unique-k,n_unique-j-k)[,1:ten_times_n_most_likely]);
            } else {
              all_subsets = combn(n_unique,n_unique-j)[,1:ten_times_n_most_likely];
            }
          } else {
            all_subsets = combn(n_unique,n_unique-j);
          }
          all_subsets = matrix(order(weights_y,decreasing = T)[all_subsets],nrow=n_unique-j,ncol=ten_times_n_most_likely)
          all_subsets = all_subsets[,order(colSums(matrix(weights_y[all_subsets],nrow=n_unique-j,ncol=ncol(all_subsets))),decreasing=T)[1:n_most_likely]];
          all_subsets = cbind(all_subsets,
                              replicate(5 * n_random_weighted,sample(n_unique,n_unique-j,prob=weights_y)));
          #replicate(n_random_lessweighted,sample(n,n-j,prob=weights_y_alt)));
          #all_subsets = replicate(round(ndraws),sample(n,n-j,prob=weights_y));
          all_subsets = apply(all_subsets,2,sort);
          all_subsets = all_subsets[,!duplicated(all_subsets,MARGIN = 2),drop = F];
          all_subsets = all_subsets[,1:min(ndraws,ncol(all_subsets)),drop = F];
        }
        all_subsets_weight = matrix(weights_y[all_subsets],nrow=n_unique-j,ncol=ncol(all_subsets));
        all_subsets = all_subsets[,order(colSums(all_subsets_weight),decreasing = T),drop=F];
        i = 1;
        while(!lowest_ub_found & (i < ncol(all_subsets))) {
          curr_subset = sort(as.numeric(unlist(affiliated_obs[all_subsets[,i]])));
          model_subset = glm.fit(x_standardized_aug[curr_subset,,drop = F], y[curr_subset], family = binomial());
          fitted_subset = drop(x_standardized_aug%*%ifelse(is.na(model_subset$coefficients),0,model_subset$coefficients));
          if(diff(range(fitted_subset)) < .Machine$double.eps^0.5) {i = i+1; next;}
          #THIS LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_subset, FUN=sum);
          foo = data.frame(y = y,fitted_subset = fitted_subset,const = 1);
          aggregate_y = foo$y;
          aggregate_fitted_subset = foo$fitted_subset;
          aggregate_count = foo[,3];
          order_increasing_y = order(aggregate_fitted_subset,aggregate_y);
          order_decreasing_y = order(aggregate_fitted_subset,-aggregate_y);
          aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_subset[order_increasing_y]));
          ncomplete_subset1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
          ncomplete_subset2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
          curr_min_ub = min(ncomplete_subset1$ub,ncomplete_subset2$ub);
          if(curr_min_ub < min_ub) {
            #cat("new min_ub =", curr_min_ub,"\n");
            min_ub = curr_min_ub;
            min_plane = coef(model_subset);
            if(ncomplete_subset1$ub <= ncomplete_subset2$ub) {
              min_subset = sort(order_increasing_y[ncomplete_subset1$sep_subseq]);
            } else {
              min_subset = sort(order_decreasing_y[ncomplete_subset2$sep_subseq]);
            }
          }
          #If we've achieved an upperbound of 1 and have already checked that separation doesn't exist in 
          #the complete dataset, then an upper bound of 1 is the best we can do
          if(min_ub <= min(j,1)) {
            lowest_ub_found = T;#Flag that separation was identified
            break;
          }
          i = i + 1;
        }
        if(lowest_ub_found | j >= loop_stop) {
          break;
        }
      }
    }
    return_vals = list(statistic = min_ub,
                       included_subset = min_subset,
                       linear_separator = min_plane,
                       solitary_separators = solitary_separators);
  }
  return(return_vals)
}


calculate_sep_stat = function(x_standardized, y, which_stat = "comp", maxk = 10, ndraws = 1e3, seed = .Machine$integer.max) {
  set.seed(seed);
  if(which_stat == "comp") {
    #calculates number of observations needed to remove to induce complete separation
    overlap_allowed = F; 
    no_intercept = F; 
  } else if(which_stat == "piv") {
    #calculates number of observations needed to remove to induce pivotal separation
    overlap_allowed = F; 
    no_intercept = T; 
  } else if(which_stat == "over") {
    #calculates number of observations needed to remove to induce quasi-complete separation
    overlap_allowed = T; 
    no_intercept = F; 
  } else {
    stop("'which_stat' must be one of 'comp', 'piv', or 'over'");
  }
  
  n = length(y);
  p = ifelse(class(x_standardized) == "numeric", 1, ncol(x_standardized));
  separating_plane = numeric(p+1);
  min_ub = n;
  lowest_ub_found = F;
  solitary_separators = NULL;
  
  if(diff(range(y)) < .Machine$double.eps^0.5) {
    warning("No variation in 'y' means that complete separation is trivially satisfied"); 
    return_vals = list(statistic = NA,
                       included_subset = NA,
                       linear_separator = NA,
                       solitary_separators = NA);
  } else if(p == 1) {#If single x caused separation, then, by definition, it was solitary
    if(no_intercept) {
      model_uni = glm(y ~ -1 + x_standardized, family = "binomial");
    } else {
      model_uni = glm(y ~ x_standardized, family = "binomial");  
    }
    fitted_uni = predict(model_uni,type='link');    
    #THIS LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_uni, FUN=sum);
    foo = data.frame(y = y,fitted_uni = fitted_uni,const = 1);
    aggregate_y = foo$y;
    aggregate_fitted_uni = foo$fitted_uni;
    aggregate_count = foo[,3];
    order_increasing_y = order(aggregate_fitted_uni,aggregate_y);
    order_decreasing_y = order(aggregate_fitted_uni,-aggregate_y);
    aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_uni[order_increasing_y]));
    ncomplete_uni1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
    ncomplete_uni2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
    min_ub = min(ncomplete_uni1$ub,ncomplete_uni2$ub);
    min_plane = coef(model_uni);
    if(ncomplete_uni1$ub <= ncomplete_uni2$ub) {
      min_subset = sort(order_increasing_y[ncomplete_uni1$sep_subseq]);
    } else {
      min_subset = sort(order_decreasing_y[ncomplete_uni2$sep_subseq]);
    }
    return_vals = list(statistic = min_ub,
                       included_subset = min_subset,
                       linear_separator = min_plane,
                       solitary_separators = 1 * (min_ub == 0));
  } else {#Otherwise assume that x_standardized is a matrix and go through each column at a time
    if(no_intercept) {
      min_plane = numeric(p);
      names(min_plane) = paste0("x_standardized",colnames(x_standardized));
    } else {
      min_plane = numeric(p+1);
      names(min_plane) = c("(Intercept)",paste0("x_standardized",colnames(x_standardized)));
    }
    for(i in 1:p) {
      if(no_intercept) {
        model_uni = glm(y ~ -1 + x_standardized[,i,drop = F], family = "binomial");
      } else {
        model_uni = glm(y ~ x_standardized[,i,drop = F], family = "binomial");
      }
      fitted_uni = predict(model_uni,type='link');
      if(diff(range(fitted_uni)) < .Machine$double.eps^0.5) {next;}
      #THE FOLLOWING LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_uni, FUN=sum);
      foo = data.frame(y = y,fitted_uni = fitted_uni,const = 1);
      aggregate_y = foo$y;
      aggregate_fitted_uni = foo$fitted_uni;
      aggregate_count = foo[,3];
      order_increasing_y = order(aggregate_fitted_uni,aggregate_y);
      order_decreasing_y = order(aggregate_fitted_uni,-aggregate_y);
      aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_uni[order_increasing_y]));
      ncomplete_uni1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
      ncomplete_uni2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
      curr_min_ub = min(ncomplete_uni1$ub,ncomplete_uni2$ub);
      if(curr_min_ub < min_ub) {
        min_ub = curr_min_ub;
        min_plane = min_plane * 0;
        #cat("new min_ub =", curr_min_ub,"\n");
        if(no_intercept) {
          min_plane[i] = coef(model_uni);
        } else {
          min_plane[c(1, 1 + i)] = coef(model_uni);
        }
        if(ncomplete_uni1$ub <= ncomplete_uni2$ub) {
          min_subset = sort(order_increasing_y[ncomplete_uni1$sep_subseq]);
        } else {
          min_subset = sort(order_decreasing_y[ncomplete_uni2$sep_subseq]);
        }
      }
      if(curr_min_ub == 0) {
        solitary_separators = c(solitary_separators,i);
      }
    }
    rm(foo,aggregate_y,aggregate_fitted_uni,aggregate_count,ncomplete_uni1,ncomplete_uni2,curr_min_ub)
    
    loop_stop = max(maxk,2);
    #If the smallest upper bound is 0, then this is trivially the lowest upperbound.
    lowest_ub_found = (min_ub == 0);
    #Now if needed go to iterative sampling of subset
    if(!lowest_ub_found & (loop_stop > 1)) {
      data_x_standardized = data.frame(x_standardized);
      
      combined_data = cbind(x_standardized,y);
      which_unique = which(!duplicated(combined_data));
      fitted_heldout = numeric(length(which_unique));
      affiliated_obs = vector("list",length(which_unique));
      names(affiliated_obs) = names(fitted_heldout) = which_unique;
      for(j in which_unique) {
        curr_remove = as.numeric(which(rowSums(combined_data == rep(combined_data[j,], each = nrow(combined_data))) == ncol(combined_data)));
        affiliated_obs[[as.character(j)]] = curr_remove;
        curr_keep = setdiff(1:n,curr_remove);
        if(no_intercept) {
          model_subset = glm(y ~ -1 + x_standardized, subset = curr_keep, family = "binomial");
        } else {
          model_subset = glm(y ~ x_standardized, subset = curr_keep, family = "binomial");
        }
        fitted_heldout[as.character(j)] = predict(model_subset,newdata=data_x_standardized,type='link')[curr_remove[1]];
      }
      if(no_intercept) {
        model_full = glm(y ~ -1 + x_standardized, family = "binomial");
      } else {
        model_full = glm(y ~ x_standardized, family = "binomial");
      }
      fitted_full = predict(model_full,type='link')[which_unique];
      weights_y = 1 / abs(fitted_heldout-fitted_full)^0.5;
      weights_y[y[which_unique]==1] = pmin(weights_y[y[which_unique]==1], sum(1-y) + 1e-3);
      weights_y[y[which_unique]==0] = pmin(weights_y[y[which_unique]==0], sum(y) + 1e-3);
      rm(model_full,fitted_full,fitted_heldout,model_subset,curr_keep,curr_remove);
      
      j = -1;
      n_unique = length(which_unique);
      n_most_likely = round(ndraws/2);#how many of the approximately most likely subsets to look at
      n_random_weighted = pmax(1, ndraws - n_most_likely);#how many to randomly sample based upon less uniform weights
      while(T) {
        j = j + 1;
        choose_n_nminusj = choose(n_unique, n_unique-j);
        ten_times_n_most_likely = pmin(choose_n_nminusj, 10 * n_most_likely);
        #cat(j,min_ub,"\n");
        if(choose_n_nminusj < ndraws) {#all possible combinations is feasible
          all_subsets = combn(n_unique,n_unique-j);
        } else {
          #combination of most likely weighted subsets
          #and weighted random sampling from among remaining subsets 
          if(choose_n_nminusj > ten_times_n_most_likely) {
            k = max((0:(n_unique-j))[mapply(FUN = "choose",n=n_unique-(0:(n_unique-j)),k=(n_unique-j)-(0:(n_unique-j))) >= ten_times_n_most_likely]);
            if(k > 0) {
              all_subsets = rbind(matrix(rep(1:k,times=ten_times_n_most_likely),nrow=k,ncol=ten_times_n_most_likely),k + combn(n_unique-k,n_unique-j-k)[,1:ten_times_n_most_likely]);
            } else {
              all_subsets = combn(n_unique,n_unique-j)[,1:ten_times_n_most_likely];
            }
          } else {
            all_subsets = combn(n_unique,n_unique-j);
          }
          all_subsets = matrix(order(weights_y,decreasing = T)[all_subsets],nrow=n_unique-j,ncol=ten_times_n_most_likely)
          all_subsets = all_subsets[,order(colSums(matrix(weights_y[all_subsets],nrow=n_unique-j,ncol=ncol(all_subsets))),decreasing=T)[1:n_most_likely]];
          all_subsets = cbind(all_subsets,
                              replicate(5 * n_random_weighted,sample(n_unique,n_unique-j,prob=weights_y)));
          #replicate(n_random_lessweighted,sample(n,n-j,prob=weights_y_alt)));
          #all_subsets = replicate(round(ndraws),sample(n,n-j,prob=weights_y));
          all_subsets = apply(all_subsets,2,sort);
          all_subsets = all_subsets[,!duplicated(all_subsets,MARGIN = 2),drop = F];
          all_subsets = all_subsets[,1:min(ndraws,ncol(all_subsets)),drop = F];
        }
        all_subsets_weight = matrix(weights_y[all_subsets],nrow=n_unique-j,ncol=ncol(all_subsets));
        all_subsets = all_subsets[,order(colSums(all_subsets_weight),decreasing = T),drop=F];
        i = 1;
        while(!lowest_ub_found & (i < ncol(all_subsets))) {
          curr_subset = sort(as.numeric(unlist(affiliated_obs[all_subsets[,i]])));
          if(no_intercept) {
            model_subset = glm(y ~ -1 + x_standardized, subset = curr_subset, family = "binomial");
          } else {
            model_subset = glm(y ~ x_standardized, subset = curr_subset, family = "binomial");
          }
          fitted_subset = predict(model_subset, newdata = data_x_standardized,type='link');
          if(diff(range(fitted_subset)) < .Machine$double.eps^0.5) {i = i+1; next;}
          #THIS LINE IS STILL IN DEVELOPMENT; #foo = aggregate(rep(1,n) ~ y + fitted_subset, FUN=sum);
          foo = data.frame(y = y,fitted_subset = fitted_subset,const = 1);
          aggregate_y = foo$y;
          aggregate_fitted_subset = foo$fitted_subset;
          aggregate_count = foo[,3];
          order_increasing_y = order(aggregate_fitted_subset,aggregate_y);
          order_decreasing_y = order(aggregate_fitted_subset,-aggregate_y);
          aggregate_fitted_label = cumsum(!duplicated(aggregate_fitted_subset[order_increasing_y]));
          ncomplete_subset1 = min_to_separate(aggregate_y[order_increasing_y],aggregate_fitted_label,aggregate_count[order_increasing_y], overlap_allowed);
          ncomplete_subset2 = min_to_separate(aggregate_y[order_decreasing_y],aggregate_fitted_label,aggregate_count[order_decreasing_y], overlap_allowed);
          curr_min_ub = min(ncomplete_subset1$ub,ncomplete_subset2$ub);
          if(curr_min_ub < min_ub) {
            #cat("new min_ub =", curr_min_ub,"\n");
            min_ub = curr_min_ub;
            min_plane = coef(model_subset);
            if(ncomplete_subset1$ub <= ncomplete_subset2$ub) {
              min_subset = sort(order_increasing_y[ncomplete_subset1$sep_subseq]);
            } else {
              min_subset = sort(order_decreasing_y[ncomplete_subset2$sep_subseq]);
            }
          }
          #If we've achieved an upperbound of 1 and have already checked that separation doesn't exist in 
          #the complete dataset, then an upper bound of 1 is the best we can do
          if(min_ub <= min(j,1)) {
            lowest_ub_found = T;#Flag that separation was identified
            break;
          }
          i = i + 1;
        }
        if(lowest_ub_found | j >= loop_stop) {
          break;
        }
      }
    }
    return_vals = list(statistic = min_ub,
                       included_subset = min_subset,
                       linear_separator = min_plane,
                       solitary_separators = solitary_separators);
  }
  return(return_vals)
}


#Identifies smallest non-trivial separation of an ordered binary sequence
#Specifically, it identifies the the smallest number of elements of the sequence
#that, if removed, would result in a monotonic subsequence of the remaining elements
#that still contains at least one element of each class
min_to_separate = function(numeric_seq, label_seq, aggregated = NULL, overlap_allowed = T) {
  stopifnot(class(numeric_seq) %in% c("integer","numeric"));
  stopifnot(min(diff(label_seq)) >= 0 & max(diff(label_seq)) > 0);
  if(is.null(aggregated)) {aggregated = 1 + numeric(length(numeric_seq));}
  if(diff(range(aggregated)) > .Machine$double.eps ^ 0.5) {
    stop("In the current version, each element of the sequence must be equally weighted");
  };
  
  #Identify shortest non-trivial, non-decreasing sequence
  first_min = min(which(numeric_seq == min(numeric_seq)));
  last_max = max(which(numeric_seq == max(numeric_seq)));
  upper_seq = first_min:last_max;
  bound1 = nondecreasing_subsequence(numeric_seq[upper_seq]);
  bound1_weight = sum(aggregated[upper_seq[bound1[[2]]]]);
  if(!overlap_allowed & any(duplicated(label_seq))) {
    label_subseq = label_seq[upper_seq[bound1[[2]]]];
    numeric_subseq = numeric_seq[upper_seq[bound1[[2]]]];
    foo = xtabs(~ numeric_subseq + label_subseq)
    overlap_by_label = colSums(foo > 0) > 1;
    if(sum(overlap_by_label) == 1) {
      if(which(overlap_by_label) == 1) {
        overlap_seq = which(label_seq[upper_seq]%in%names(which(overlap_by_label)) & !numeric_seq[upper_seq]%in%names(foo[,which(overlap_by_label)][1]));
      } else if(which(overlap_by_label) == ncol(foo)) {
        overlap_seq = which(label_seq[upper_seq]%in%names(which(overlap_by_label)) & !numeric_seq[upper_seq]%in%names(foo[,which(overlap_by_label)][nrow(foo)]));
      } else {
        overlap_seq = which(label_seq[upper_seq]%in%names(which(overlap_by_label)) & !numeric_seq[upper_seq]%in%names(which.max(foo[,which(overlap_by_label)])));
      }
      bound1[[2]] = setdiff(bound1[[2]],overlap_seq)
      bound1[[1]] = length(bound1[[2]]);
      bound1_weight = sum(aggregated[upper_seq[bound1[[2]]]]);
    } else if(sum(overlap_by_label) > 1) {
      stop("check code: minimum non-increasing sequence has more than one point of overlap remaining");
    }
  }
  
  #Identify shortest non-trivial, non-increasing sequence
  first_max = min(which(numeric_seq == max(numeric_seq)));
  last_min = max(which(numeric_seq == min(numeric_seq)));
  lower_seq = first_max:last_min;
  bound2 = nondecreasing_subsequence(-1*numeric_seq[lower_seq]);
  bound2_weight = sum(aggregated[lower_seq[bound2[[2]]]]);
  if(!overlap_allowed & any(duplicated(label_seq))) {
    label_subseq = label_seq[lower_seq[bound2[[2]]]];
    numeric_subseq = numeric_seq[lower_seq[bound2[[2]]]];
    foo = xtabs(~ numeric_subseq + label_subseq)
    overlap_by_label = colSums(foo > 0) > 1;
    if(sum(overlap_by_label) == 1) {
      if(which(overlap_by_label) == 1) {
        overlap_seq = which(label_seq[lower_seq]%in%names(which(overlap_by_label)) & !numeric_seq[lower_seq]%in%names(foo[,which(overlap_by_label)][nrow(foo)]));
      } else if(which(overlap_by_label) == ncol(foo)) {
        overlap_seq = which(label_seq[lower_seq]%in%names(which(overlap_by_label)) & !numeric_seq[lower_seq]%in%names(foo[,which(overlap_by_label)][1]));
      } else {
        overlap_seq = which(label_seq[lower_seq]%in%names(which(overlap_by_label)) & !numeric_seq[lower_seq]%in%names(which.max(foo[,which(overlap_by_label)])));
      }
      bound2[[2]] = setdiff(bound2[[2]], overlap_seq)
      bound2[[1]] = length(bound2[[2]]);
      bound2_weight = sum(aggregated[lower_seq[bound2[[2]]]]);
    } else if(sum(overlap_by_label) > 1) {
      stop("check code: minimum non-increasing sequence has more than one point of overlap remaining");
    }
  }
  
  if(bound1_weight > bound2_weight) {
    #The longest subsequence, i.e. fewest number of elements to remove, is non-decreasing
    sep_subseq = upper_seq[bound1[[2]]];
    return_vals = list(ub = sum(aggregated) - bound1_weight,
                       sep_subseq = sep_subseq,
                       removed = setdiff(1:length(numeric_seq),sep_subseq));
  } else {
    #The longest subsequence, i.e. fewest number of elements to remove, is non-increasing
    sep_subseq = lower_seq[bound2[[2]]];
    return_vals = list(ub = sum(aggregated) - bound2_weight,
                       sep_subseq = sep_subseq,
                       removed = setdiff(1:length(numeric_seq),sep_subseq));  
  }
  return_vals;
}


#Adapted from wikipedia page
nondecreasing_subsequence = function(numeric_seq) {
  
  n = length(numeric_seq);
  P = numeric(n);
  M = numeric(n+1);
  
  L = 0;
  for (i in 1:n) {
    #// Binary search for the largest positive j <= L
    #// such that numeric_seq[M[j]] <= numeric_seq[i]
    lo = 1;
    hi = L;
    while(lo <= hi){
      mid = ceiling((lo+hi)/2);
      if (numeric_seq[M[mid+1]+1] <= numeric_seq[i]) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }
    #// After searching, lo is 1 greater than the
    #// length of the longest prefix of numeric_seq[i]
    newL = lo;
    
    #// The predecessor of numeric_seq[i] is the last index of 
    #// the subsequence of length newL-1
    P[i] = M[newL];
    M[newL+1] = i-1;
    
    #// If we found a subsequence longer than any we've
    #// found yet, update L
    if (newL > L)
      L = newL;
  }
  #// Reconstruct the index of the longest non-decreasing subsequence
  S = numeric(L);
  k = M[L+1]+1;
  for(i in L:1) {
    S[i] = k;
    k = P[k]+1;
  };
  list(L, S);
}



solve_for_exppow_scale = function(abs_x,mass_above_abs_x,shape) {
  stopifnot(mass_above_abs_x < 1 & mass_above_abs_x > 0 & abs_x >0);
  foo = function(t) {(pgamma((abs_x/sqrt(2)/t)^shape,1/shape,lower=F) - mass_above_abs_x);}
  uniroot(foo, lower=0, upper = 100)$root;
}


solve_for_logistic_scale = function(abs_x,mass_above_abs_x) {
  stopifnot(mass_above_abs_x < 1 & mass_above_abs_x > 0 & abs_x >0);
  abs_x/log((1-mass_above_abs_x/2)/(mass_above_abs_x/2));
}

solve_for_hiershrink_scale = function(target_mean,
                                      local_dof = 1,
                                      global_dof = 1,
                                      npar = 0,
                                      n,
                                      sigma = 2,
                                      tol = .Machine$double.eps^0.5,
                                      max_iter = 100, 
                                      nsim = 2e5,
                                      slab_precision = 1/(15^2)
) {
  
  
  stopifnot(target_mean > 0 & target_mean < npar);
  lambda = matrix(rt(nsim*npar,df=local_dof),nrow=nsim) * rt(nsim,df=global_dof);

  log_tau0 = diff_target = numeric(max_iter);
  log_tau0[1] = log(target_mean/(npar - target_mean)*sigma/sqrt(n));
  
  random_scales = 1 / (slab_precision + 1/(exp(2*log_tau0[1]) * lambda^2));
  kappa = 1/(1+n*random_scales/sigma^2);
  
  diff_target[1] = mean(rowSums(1-kappa)) - target_mean;
  log_tau0[2] = 0.02 + log_tau0[1];
  random_scales = 1 / (slab_precision + 1/(exp(2*log_tau0[2]) * lambda^2));
  
  kappa = 1/(1+n*random_scales/sigma^2);
  diff_target[2] = mean(rowSums(1-kappa)) - target_mean;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_tau0[i] = log_tau0[i-1] - diff_target[i-1]*(log_tau0[i-1]-log_tau0[i-2])/(diff_target[i-1]-diff_target[i-2]);
    random_scales = 1 / (slab_precision + 1/(exp(2*log_tau0[i]) * lambda^2));
    kappa = 1/(1+n*random_scales/sigma^2);
    diff_target[i] = mean(rowSums(1-kappa)) - target_mean;
    if(abs(diff_target[i]-diff_target[i-1])<tol) {break;}
  }
  return(list(tau0 = exp(log_tau0[i]),
              diff_from_target = diff_target[i],
              iter = i));
}
