# Intercept Priors

### Current Suggested Citation

Boonstra, Philip S., Barbaro, Ryan P., and Sen, Ananda. "Default Priors for the Intercept Parameter in Logistic Regressions" (March 2018). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 123.
https://biostats.bepress.com/umichbiostat/paper123

## Executive Summary

The functions <samp>glm_HSBeta</samp> and <samp>glm_LogisBeta</samp> implement the primary methodological development in this manuscript. They are wrappers for code written in STAN that implements Bayesian logistic regressions using horseshoe-type or logistic-type priors on the regression coefficients $\beta$, respectively. Most importantly for this paper, however, they allow the user to choose from a number of priors on the intercept parameter $\alpha$, including a student-$t$, logistic, or exponential-power prior. 

## Further details

In more detail, there are six files included in this repository (in addition to this README): one text file (ending in <samp>.txt</samp>), three <samp>R</samp> scripts (ending in <samp>.R</samp>), and two STAN functions (ending in <samp>.stan</samp>). The simulation studies reported in Boonstra and Barbaro were run using commit 3.

### Text file

<samp>runIntPriorSims.txt</samp> is the script for submitting parallel runs of <samp>IntPriorSims.R</samp> (described below) to a cluster that is running SLURM. The following command does this:

<code> sbatch runIntPriorSims.txt </code>

### <samp>R</samp> files

<samp>Functions.R</samp> provides all of the necessary functions to fit the methods described in the paper, including but not limited to code for fitting the logistic regressions under a variety of priors on $\beta$ and $\alpha$ and code for executing both Algorithms proposed in the manuscript. 

<samp>GenParams.R</samp> constructs inputs for running the simulation study. As described in the script's documentation and the language below, these inputs can be overwritten by the user.

<samp>IntPriorSims.R</samp> is the script to conduct the large-scale simulation study described in the manuscript. On a local machine, the user may choose a specific <samp>array_id</samp> (as described in this script's documentation) and run the code locally on his/her machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously (as described above). 


### STAN files
The STAN files are described below. These currently all implement a logistic link; changing to a non-logistic link (i.e. log, probit, etc.) would be relatively easy. Upon using these for the first time, <samp>R</samp> will need to compile these programs, creating an <samp>R</samp> data object file (ending in <samp>.rds</samp>) in the current working directory. Re-compilation of the STAN files are not necessary as long as they are unchanged.

<samp>HSBeta.stan</samp> implements a Bayesian logistic regression with a regularized horseshoe prior applied to the regression coefficients and one of several possible priors on the intercept, depending on the user's selections. An <samp>R</samp> user calls this with <samp>glm_HSBeta</samp> in <samp>Functions.R</samp>. 

<samp>LogisBeta.stan</samp> implements a Bayesian logistic regression with independent, standard logistic priors on each of the regression coefficients and one of several possible priors on the intercept, depending on the user's selections. An <samp>R</samp> user calls this with <samp>glm_LogisBeta</samp> in <samp>Functions.R</samp>. 
