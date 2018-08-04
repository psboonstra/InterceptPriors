#Code to run the simulation study:
#Calling 'GenParams.R' creates a list called 'arglist' of length 950:
#The 22 p<= 75 scenarios are replicated 25 times (for embarrassingly parallel execution)
#The 8 p== 150 scenarios are replicated 50 times (for embarrassingly parallel execution)
#Choose any value of 1<=array_id<=950 to run that particular scenario

rm(list=ls(all=TRUE));
library(rstan);
library(pROC);
library(shinystan);
my_computer = T;
if(my_computer) {
  array_id = 1;
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}

rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

stan_file_path = "";#where are the stan files stored?
source_file_path = "";#where are the other R scripts stored?
hsbeta_file_name = "HSBeta.stan";
logisbeta_file_name = "LogisBeta.stan";

file.name = "EstInterceptNewMetrics";
priors_to_fit = c(
  "StudT10",
  "Normal10",
  "EP2",
  "EP10A",
  "EP10F"
);

#arglist is a list of lists, each of which contains
#different simulation settings. 

source(paste0(source_file_path,"GenParams.R"));

rm(list=setdiff(ls(),c("arglist","my_computer","file.name","array_id","source_file_path")))
source(paste0(source_file_path,"Functions.R"));

curr_args = arglist[[array_id]];
curr_args$n = 600;

timings  = system.time(assign(paste("sim",array_id,sep=""),do.call("InterceptSim",args=curr_args)));
cat(array_id,":",as.numeric(timings[3])/3600,"hours\n",file="results.txt",append=T);

save(list=paste("sim",array_id,sep=""),file=paste(file.name,array_id,".RData",sep=""));
