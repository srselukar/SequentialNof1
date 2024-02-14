##############################################
# Functions to Run Simulation Study
# Effect of Boundary Shape on Bias/MSE of
# Combined Meta-Analysis Estimator
# 
# Subodh Selukar
# 2021-09-03
##############################################

# -------------------------------- Description ------------------------------------------------ #

##### This simulation study is intended to describe 
##### early stopping with sequential monitoring boundaries 
##### for a series of N-of-1 trials with 2 treatments

#### Goal is to describe combined point estimator for trials with OBF boundary shape
#### across various effect sizes and numbers of trials

##### Analysis is done with linear mixed effects model with inverse-variance weights
##### with random intercepts (i.e., each trial has its own true effect)
##### Point estimate is the fixed intercept

##### True model is actually fixed-effects model

##### NOTE: Simulation is required rather than numerical results 
##### because very complex data generation mechanisms

##### VREY SPECIFIC DATA STRUCTURES EXPECTED
##### Every setting must have an identical data frame of estimates/standard error from OBF-monitored trials



#### Requires packages
library(metafor)
library(gtools)

# -------------------------------- Functions ------------------------------------------------ #

### Function to read in OBF datasets and output estimates/variances for each series of N-of-1 trials
# Inputs:
# 1. obfDat: OBF dataset (data set with 12 columns, columns 1 and 7 should be MUE.SM and SE(MUE.SM), respectively)
# 2. sizeSeries: number of trials in series (assumed to be compatible such that nrow(obfDat)/sizeSeries is an integer)

# Outputs:
# vector of combined point estimates with length=nrow(obfDat)/sizeSeries
simSeries <- function(obfDat,
                      sizeSeries,numOBF){
    numSeries <- nrow(obfDat)/sizeSeries # groups the independent trials from input according to desired series size 
    
    out <- rep(NA,numSeries)
    for (i in 1:numSeries){ # definitely a slow way to do it: there's probably a way to use closed-form estimator from meta-analysis model + sample from all possible permutations...
        
        tmpDat <- obfDat[(sizeSeries*(i-1)+1):(i*sizeSeries),c(6,12)] # use the naive estimator
        
        tmpOut <- rma(yi=tmpDat[,1],sei=tmpDat[,2], control=list(maxiter=500,stepadj=0.5)) # defaults to using inverse variance weighting with REML
        
        out[i] <- tmpOut$beta[1]
    }
    
    return(out)
}
