# Intercept Priors

### Current Suggested Citation

Boonstra, Philip S., Barbaro, Ryan P., and Sen, Ananda. "Default Priors for the Intercept Parameter in Logistic Regressions" (March 2018). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 123.
https://biostats.bepress.com/umichbiostat/paper123

## Executive Summary

The functions <samp>glm_HSBeta</samp> and <samp>glm_LogisBeta</samp> implement the primary methodological development in this manuscript. They are wrappers for code written in STAN that implements Bayesian logistic regressions using horseshoe-type or logistic-type priors on the regression coefficients $\beta$, respectively. Most importantly for this paper, however, they allow the user to choose from a number of priors on the intercept parameter $\alpha$, including a student-$t$, logistic, or exponential-power prior. 

## Further details

In total there are 

### Text file

### <samp>R</samp> files

<samp>Functions.R</samp> provides all of the necessary functions to fit the methods described in the paper. 

<samp>GenParams.R</samp> constructs inputs for running the simulation study. As described in the script's documentation and the language below, these inputs can be overwritten by the user.

### STAN files
