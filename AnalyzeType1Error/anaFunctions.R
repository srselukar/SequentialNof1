##############################################
# Functions to Analyze Simulation Study
# Number of Blocks for Sequential Monitoring
# 
# Subodh Selukar
# 2021-03-03
##############################################

# -------------------------------- Description ------------------------------------------------ #

##### This simulation study is intended to describe 
##### early stopping with sequential monitoring boundaries 
##### for N-of-1 trials with 2 treatments

#### Goal is to assess number of blocks/periods 
#### needed for consequential early stopping 
#### across various effect sizes and monitoring rules (stop times and boundary types)

##### Analysis is done with linear mixed effects model
##### with random intercepts (i.e., exchangeable correlation)
##### Uses Wald statistic from REML (less biased estimate of SE than ML) 

##### That model is true: 
##### data are generated as independent blocks
##### each with input L periods

##### The periods may have correlation, so blcoks are generated 
##### via mvtnorm package

#### NOTE: Simulation is required rather than numerical results 
#### because (1) treatment assignment is random 
#### (not computationally efficient to look across
#### all possible combinations of within-block and between-block
#### treatment sequences)
#### otherwise, could use numerical integration of density
#### and (2) using observed information to calculate test statistics
#### other results show that, with small samples, observed
#### vs. expected information may greatly affect properties

#### Requires package
# gsDesign (monitoring boundaries)

# -------------------------------- Functions ------------------------------------------------ #

### Function to create monitoring boundaries (requires gsDesign)
## Right now only supporting symmetric bounds

# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. totBlocks: total number of blocks (max "sample size")
# 3. boundType: type of boundary (defaults to Eff; supporting OBF, Pocock and effect size-based boundaries)
# 4. alpha: one-sided alpha level (defaults to 0.025, not planning to change)
# 5. effParams: vector of length 3 of the target effect, error variance and 
#   number of periods per block when using the Eff or stdEff boundaries 

# Outputs: 
# matrix with nrow=number of stopping times and 2 columns with the monitoring boundaries at each stopping time
makeBounds <- function(stopTimes,
                       totBlocks,
                       boundType="Eff",
                       alpha=0.025,
                       effParams=c(1,1,4)){
    
    if (length(stopTimes)==1){ # need to have special case when there's only one look
        out <- matrix(qnorm(c(alpha,1-alpha)),
                      nrow=1,ncol=2)
        return(out)
    }
    
    if (boundType=="Poc"){
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="Pocock")$upper$bound
        
        out <- cbind(-upBound,upBound)
    } else if (boundType=="OBF"){
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="OF")$upper$bound
        
        out <- cbind(-upBound,upBound)
    } else if (boundType=="Asy"){
        loBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="Pocock")$lower$bound
        
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="OF")$upper$bound
        
        out <- cbind(loBound,upBound)
    } else if (boundType=="stdEff") { # stdEff is a Pocock boundary, but specifying the effect size rather than alpha
        stdBound <- effParams[1]/sqrt(effParams[2]*2/effParams[3]) # this standardizes the effect size to the scale of the treatment effect estimator with expected info
        
        out <- cbind(rep(-stdBound,length(stopTimes)), # lower bound
                     rep(stdBound,length(stopTimes))) # upper bound
    } else if (boundType=="Eff") { # uses raw effect size rather than standardizing by expected information
        stdBound <- effParams[1]
        
        out <- cbind(rep(-stdBound,length(stopTimes)), # lower bound
                     rep(stdBound,length(stopTimes))) # upper bound
    } else stop("Not supported")
    
    return(out)
}

### Function to convert from z-boundaries to (approximate) t-boundaries (per Nikolakopoulos 2016)
# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. zIn: input Z critical values, as made by makeBounds()

# Outputs: 
# matrix with nrow=number of stopping times and 2 columns with the monitoring boundaries at each stopping time
makeZtoTBounds <- function(stopTimes,
                       zIn){
    
    out <- cbind(
        qt(pnorm(zIn[,1]),stopTimes-2),
        qt(pnorm(zIn[,2]),stopTimes-2)
    )
    out[,1] <- ifelse(is.nan(out[,1]),-Inf,out[,1]) # need to deal with NaN
    out[,2] <- ifelse(is.nan(out[,2]),Inf,out[,2])
    
    return(out)
}



### Function to use input monitoring boundary and simulation results to derive stopping and error probabilities

## Helper function for logical testing that trial stopped at this block (and not before)
# Inputs: 
# 1. x: sequence of estimates (scale depends on boundary scale)
#   length(x) > 1 (do not use if only 1 estimate); assumes sequential blocks! (can't be out of order)
# 2. lower: sequence of lower bounds with length=length(x)
# 3. upper: sequence of upper bounds with length=length(x)
logReached <- function(x,lower,upper){
    nowBlock <- length(x)
    all(x[-nowBlock] > lower[-nowBlock] & x[-nowBlock] < upper[-nowBlock]) & 
        (x[nowBlock] < lower[nowBlock] | x[nowBlock] > upper[nowBlock])
}

## Inputs:
# 1. dataIn: simulation results - a matrix of the estimates and standard errors with 
#   nrow() = number of simulations
#   ncol() = 2*(number of trial blocks-1) with 
#   odd columns the estimates after each block (no analysis after block 1) and 
#   even columns the (observed) standard error
# 2. stopTimes: vector of stopping times (i.e., blocks after which data would have been analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 3. boundIn: a matrix of the monitoring boundaries with 
#   nrow() = length(stopTimes)
#   ncol() = 2 for upper and lower boundaries (currently only supporting symmetric rules)
# 4. effBound: indicator of using effect size boundary; i.e., NOT on Z-scale (don't compute Wald statistic)

## Outputs: list of
# 1. error probability under truth (i.e., power with truEff != 0 or type-1 error with truEff==0) and under alternatives
# 2. stopping probabilities at each stopping time
#  these are conditional on continuing to the given stopping time

## Important: block 1 estimate/SE not included in simulation results - impacts the indexing

anaSims <- function(dataIn,stopTimes,boundIn,effBound=FALSE){
    
    ### Calculate error and stopping probabilities
    # Pick out the statistics at stopping times
    
    if (length(stopTimes)==1) { # need to deal specially with 1 look because later code expects a matrix
        statIn <- matrix(
            dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)], # transform the simulation data to wald statistic for monitoring boundary
            nrow=nrow(dataIn),ncol=length(stopTimes)
        )
    } else {
        statIn <- dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)] # transform the simulation data to wald statistic for monitoring boundary
    }
    
    
    stopProbs <- rep(NA,length(stopTimes))
    fullStop <- rep(FALSE,nrow(statIn)) # will hold if the trial crossed boundary; starts with no 
    for (j in 1:length(stopTimes)){ # loop over each stopping time
        tmp <- statIn[,1:j] # statistics until current block/possible stopping time j
        
        if (j==1) {
            tmpStop <- (tmp < boundIn[j,1] | tmp > boundIn[j,2])
        } else {
            tmpStop <- apply(X=tmp,MARGIN=1,FUN=logReached,
                             lower=boundIn[1:j,1],upper=boundIn[1:j,2])
        }
        stopProbs[j] <- mean(tmpStop) # mean of how many trials crossed boundary at this stopping time
        
        fullStop <- fullStop | tmpStop # did trial stop now or had it before?
    }
    errProb <- mean(fullStop) # mean of how many trials crossed boundary at any point
    
    return(list(
        errProb,
        stopProbs
    ))
}


