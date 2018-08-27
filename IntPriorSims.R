#Code to run the simulation study:
#Calling 'GenParams.R' creates a list called 'arglist' of length 1625 for embarrassingly parallel execution:
#The 25 p<= 75 scenarios are replicated across 25 jobs, with 8 iterations per job
#The 8 p == 150 scenarios are replicated across 100 jobs, with 2 iterations per job
#The simulation study in Boonstra and Barbaro was conducted in two batches of 625 jobs followed by another 800 (because the cluster only accepts 
#jobs in batches of size up to 1000). 'which_run' (line 16) indicates whether this is the first batch (= 1 ) or the second batch ( = 2)

rm(list=ls(all=TRUE));
library(rstan);
library(e1071)
library(pROC);
library(shinystan);
#Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = T;

which_run = 1;

if(my_computer) {
  array_id = 1;
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}

#Recommended options from rstan:
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

stan_file_path = "";#where are the stan files stored?
source_file_path = "";#where are the other R scripts stored?
hsbeta_file_name = "HSBeta.stan";
logisbeta_file_name = "LogisBeta.stan";

file_name = "Sim";
priors_to_fit = c(
  "StudTFixed",
  "NormalFixed",
  "EPTwoAdapt",
  "EPTenAdapt",
  "EPTenFixed",
  "LogisAdapt"
);


if(which_run == 1) {#first batch
  array_id_offset = 0;
} else {#second batch
  array_id_offset = 625;  
}

source(paste0(source_file_path,"GenParams.R"));

rm(list=setdiff(ls(),c("arglist","my_computer","file_name","array_id","array_id_offset","source_file_path")))
source(paste0(source_file_path,"Functions.R"));

curr_args = arglist[[array_id + array_id_offset]];

timings  = system.time(assign(paste("sim",array_id_offset + array_id,sep=""),do.call("InterceptSim",args=curr_args)));
cat(array_id,":",as.numeric(timings[3])/3600,"hours\n",file="results.txt",append=T);

save(list=paste("sim",array_id_offset + array_id,sep=""),file=paste(file_name,array_id_offset + array_id,".RData",sep=""));
