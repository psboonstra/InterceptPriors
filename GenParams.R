#################
#April 15, 2016
#Phil Boonstra
#Code to generate a list of lists of simulation settings. 
################

################
#Common settings
if(!"priors_to_fit"%in%ls()){priors_to_fit=NULL;}
if(!"fit_methods"%in%ls()){fit_methods=T;}
if(!"calc_n_over"%in%ls()){calc_n_over=F;}
if(!"calc_n_comp"%in%ls()){calc_n_comp=T;}
if(!"calc_n_piv"%in%ls()){calc_n_piv=T;}
if(!"algorithm2_maxk"%in%ls()){algorithm2_maxk = 8;}
if(!"algorithm2_ndraws"%in%ls()){algorithm2_ndraws = 2.5e3;}
if(!"n_new"%in%ls()){n_new=500;}
if(!"mc_chains"%in%ls()){mc_chains = 2;}
if(!"mc_iter_after_warmup"%in%ls()){mc_iter_after_warmup = 2e3;}
if(!"mc_warmup"%in%ls()){mc_warmup = 2e3;}
if(!"mc_thin"%in%ls()){mc_thin = 1;}
if(!"mc_max_treedepth"%in%ls()){mc_max_treedepth = 18;}
if(!"mc_adapt_delta"%in%ls()){mc_adapt_delta = 0.9999;}
if(!"mc_stepsize"%in%ls()){mc_stepsize = 0.1;}
if(!"ntries_per_iter"%in%ls()){ntries_per_iter = 3;}
if(!"local_dof"%in%ls()){local_dof = 1;}
if(!"global_dof"%in%ls()){global_dof = 1;}
if(!"slab_precision"%in%ls()){slab_precision = (1/15)^2;}
if(!"exppow_prior_args"%in%ls()){
  exppow_prior_args = list(EP2=list(alpha_scale = NULL,alpha_power = 2,adaptive_scale = T),
                           EP10A=list(alpha_scale = NULL,alpha_power = 10,adaptive_scale = T),
                           #Fixed scale with 1% of prior mass following outside of |logit(1/1e-4)|=~9.21
                           EP10F=list(alpha_scale = 6.218102,alpha_power = 10,adaptive_scale = F));
}

#Settings that vary
if(!"all_settings_list"%in%ls()){
  all_settings_list = list(
    #Simple scenario for testing code
    #1
    list(scenario = 1, 
         hiershrink_prior_num_relevant = 2,
         true_alpha = -2.5, true_betas = rep(0.5,4),
         x_binom = 1:4, case_control_ratio = NA,
         n = 50 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.5,nrow=4,ncol=4)+diag(0.5,4))),
    list(scenario = 1, 
         hiershrink_prior_num_relevant = 2,
         true_alpha = -2.5, true_betas = rep(0.5,4),
         x_binom = 1:4, case_control_ratio = NA,
         n = 100 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.5,nrow=4,ncol=4)+diag(0.5,4))),
    list(scenario = 1, 
         hiershrink_prior_num_relevant = 2,
         true_alpha = -2.5, true_betas = rep(0.5,4),
         x_binom = 1:4, case_control_ratio = NA,
         n = 200 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.5,nrow=4,ncol=4)+diag(0.5,4))),
    list(scenario = 1, 
         hiershrink_prior_num_relevant = 2,
         true_alpha = -2.5, true_betas = rep(0.5,4),
         x_binom = 1:4, case_control_ratio = NA,
         n = 400 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.5,nrow=4,ncol=4)+diag(0.5,4))),
    #2
    list(scenario = 2, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = c(1.5,rep(0,24)), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 50 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    list(scenario = 2, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = c(1.5,rep(0,24)), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 100 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    list(scenario = 2, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = c(1.5,rep(0,24)), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 200 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    #3
    list(scenario = 3, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = rep(0.06,25), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 50 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    list(scenario = 3, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = rep(0.06,25), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 100 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    list(scenario = 3, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-2,true_betas = rep(0.06,25), 
         x_binom = 1:25,case_control_ratio = NA,
         n = 200 , mu_latent_x = rep(qnorm(0.25),length=25), chol_latent_x = chol(matrix(0.15,nrow=25,ncol=25)+diag(0.85,25))),
    #4
    list(scenario = 4, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-6.5,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 100 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    list(scenario = 4, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-6.5,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 200 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    list(scenario = 4, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-6.5,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 400 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    #5
    list(scenario = 5, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 50 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    list(scenario = 5, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 100 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    list(scenario = 5, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(rep(3,10),rep(0,15)),
         x_binom = 1:25,case_control_ratio = NA,
         n = 200 , mu_latent_x = rep(qnorm(.05),length=25), chol_latent_x = chol(matrix(0.3,nrow=25,ncol=25)+diag(0.7,25))),
    ##Cohort sampling; rare disease; many small OR; n=75;p=50; normal covariates
    #6
    list(scenario = 6, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(2,rep(0,74)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 100 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    list(scenario = 6, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(2,rep(0,74)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 200 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    list(scenario = 6, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-4,true_betas = c(2,rep(0,74)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 400 , mu_latent_x = NULL, chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    ##Cohort sampling; rare disease; many small OR; n=75;p=50; normal covariates
    #7
    list(scenario = 7, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3.5,true_betas = c(2,rep(0,74)), 
         x_binom = 1:75,case_control_ratio = NA,
         n = 100 , mu_latent_x = rep(qnorm(0.25),75), chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    list(scenario = 7, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3.5,true_betas = c(2,rep(0,74)), 
         x_binom = 1:75,case_control_ratio = NA,
         n = 200 , mu_latent_x = rep(qnorm(0.25),75), chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    list(scenario = 7, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3.5,true_betas = c(2,rep(0,74)), 
         x_binom = 1:75,case_control_ratio = NA,
         n = 400 , mu_latent_x = rep(qnorm(0.25),75), chol_latent_x = chol(matrix(0.30,nrow=75,ncol=75)+diag(0.70,75))),
    ##Cohort sampling; common disease; two small-moderate OR; many normal covariates with small effects
    #8
    list(scenario = 8, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = c(-0.5,-0.5,rep(0,148)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 100 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 8, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = c(-0.5,-0.5,rep(0,148)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 200 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 8, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = c(-0.5,-0.5,rep(0,148)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 400 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 8, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = c(-0.5,-0.5,rep(0,148)), 
         x_binom = NULL,case_control_ratio = NA,
         n = 600 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    ##Cohort sampling; common disease; many small risk factors
    #9
    list(scenario = 9, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = rep(-1/150,150),
         x_binom = NULL,case_control_ratio = NA,
         n = 100 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 9, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = rep(-1/150,150),
         x_binom = NULL,case_control_ratio = NA,
         n = 200 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 9, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = rep(-1/150,150),
         x_binom = NULL,case_control_ratio = NA,
         n = 400 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150))),
    list(scenario = 9, 
         hiershrink_prior_num_relevant = 10,
         true_alpha=-3,true_betas = rep(-1/150,150),
         x_binom = NULL,case_control_ratio = NA,
         n = 600 , mu_latent_x = NULL, chol_latent_x = chol(matrix(.1,nrow=150,ncol=150)+diag(.9,150)))
  );
}

if(!"arglist"%in%ls()){arglist = list();}

if(!"rep_settings"%in%ls()){rep_settings = c(rep(25,22),rep(50,8));}
stopifnot(length(all_settings_list) == length(rep_settings));
for(i in length(rep_settings):1) {
  all_settings_list = c(all_settings_list[-i],rep(all_settings_list[i],rep_settings[i]))
}
all_settings_list = rev(all_settings_list);

#sim_number = length(arglist)+1;
random_seed = sample(.Machine$integer.max - 1000,length(all_settings_list));

for(i in 1:length(all_settings_list)) {
  true_alpha = all_settings_list[[i]]$true_alpha;
  true_betas = all_settings_list[[i]]$true_betas;
  x_binom = all_settings_list[[i]]$x_binom;
  case_control_ratio = all_settings_list[[i]]$case_control_ratio;
  n = all_settings_list[[i]]$n;
  mu_latent_x = all_settings_list[[i]]$mu_latent_x;
  chol_latent_x = all_settings_list[[i]]$chol_latent_x;
  hiershrink_prior_num_relevant = all_settings_list[[i]]$hiershrink_prior_num_relevant;
  scenario = all_settings_list[[i]]$scenario;
  
  curr_params = list(niter = ifelse("niter"%in%ls(),niter,ifelse(length(true_betas) <= 100, 8, 4)),
                     n = n,
                     n_new = n_new,
                     true_alpha = true_alpha,
                     true_betas = true_betas,
                     x_binom = x_binom,
                     case_control_ratio = case_control_ratio,
                     mu_latent_x = mu_latent_x,
                     chol_latent_x = chol_latent_x,
                     priors_to_fit = priors_to_fit,
                     hiershrink_prior_num_relevant = hiershrink_prior_num_relevant,
                     local_dof = local_dof,
                     global_dof = global_dof,
                     slab_precision = slab_precision,
                     exppow_prior_args = exppow_prior_args,
                     stan_file_path = stan_file_path,
                     hsbeta_file_name = hsbeta_file_name,
                     logisbeta_file_name = logisbeta_file_name,
                     calc_n_over = calc_n_over,
                     calc_n_comp = calc_n_comp,
                     calc_n_piv = calc_n_piv,
                     algorithm2_maxk = algorithm2_maxk,
                     algorithm2_ndraws = algorithm2_ndraws,
                     mc_chains = mc_chains, 
                     mc_iter_after_warmup = mc_iter_after_warmup,
                     mc_warmup = mc_warmup,
                     mc_thin = mc_thin, 
                     mc_max_treedepth = mc_max_treedepth,
                     mc_adapt_delta = mc_adapt_delta,
                     mc_stepsize = mc_stepsize,
                     ntries_per_iter = ntries_per_iter,
                     random_seed = random_seed[i],
                     array_id = i,
                     scenario = scenario,
                     fit_methods = fit_methods,
                     data_analysis = F,
                     verbose = F)
  arglist = c(arglist,list(curr_params));
  rm(true_alpha,true_betas,x_binom,case_control_ratio,n,mu_latent_x,chol_latent_x,hiershrink_prior_num_relevant);
  #sim_number = sim_number + 1;
}
################


