##############################################
# Functions to Analyze Simulation Study
# Average Sample Number, 
# Probability of Stopping Early,
# Probability of Correct Decision
# 
# Subodh Selukar
# 2021-10-12
##############################################

### 2023-04-19: update Z to T function to b-2 df, consistent with Nikolopaulos 2016

# -------------------------------- Description ------------------------------------------------ #

##### This simulation study is intended to describe 
##### early stopping with sequential monitoring boundaries 
##### for N-of-1 trials with 2 treatments

#### Goal is to assess 
#### Average Sample Number, 
#### Probability of Stopping Early,
#### Probability of Correct Decision
#### across various trial sizes and monitoring rules (stop times and boundary types)

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
## Right now, only selected boundary shapes

# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. totBlocks: total number of blocks (max "sample size")
# 3. boundType: type of boundary (supporting symmetric OBF, symmetric Pocock and asymmetric with OBF upper and Pocock lowerr)
# 4. alpha: one-sided alpha level (defaults to 0.025, not planning to change)

# Outputs: 
# matrix with nrow=number of stopping times and 2 columns with the monitoring boundaries at each stopping time
makeBounds <- function(stopTimes,
                       totBlocks,
                       boundType="symOBF",
                       alpha=0.025){
    
    if (length(stopTimes)==1){ # need to have special case when there's only one look
        out <- matrix(qnorm(c(alpha,1-alpha)),
                      nrow=1,ncol=2)
        return(out)
    }
    
    if (boundType=="symPoc"){
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="Pocock")$upper$bound
        
        out <- cbind(-upBound,upBound)
    } else if (boundType=="symOBF"){
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="OF")$upper$bound
        
        out <- cbind(-upBound,upBound)
    } else if (boundType=="asym") { # asymmetric with OBF upper and Poc lower; NOTE: allows for conclusion of null at last analysis (if no crossing at last stage)
        upBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="OF")$upper$bound
        loBound <- gsDesign(k=length(stopTimes),
                            timing=stopTimes/totBlocks,
                            test.type=2,
                            alpha=alpha,
                            sfu="Pocock")$lower$bound
        
        
        out <- cbind(loBound, # lower bound
                     upBound) # upper bound
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


### Function to use input monitoring boundary and simulation results to calculate
# 1. probability of correct decision/direction
# 2. probability of any early stopping (0% if no seq monitoring)
# 3. average sample number

## Helper function for logical testing that trial stopped at this block (and not before)
# Inputs: 
# 1. x: sequence of estimates (scale depends on boundary scale)
#   length(x) > 1 (do not use if only 1 estimate)
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
# 4. truEffect: true effect that generated the data

## Outputs: list of
# 1. overall stopping probability
# 2. stopping probabilities at each stopping time
#   conditional on continuing to the given stopping time

## Important: block 1 estimate/SE not included in simulation results - impacts the indexing

anaSims <- function(dataIn,stopTimes,boundIn,truEffect){
    
    # Pick out the statistics at stopping times
    if (length(stopTimes)==1) { # need to deal specially with 1 look because later code expects a matrix
        statIn <- matrix(
            dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)], # transform the simulation data to wald statistic for monitoring boundary
            nrow=nrow(dataIn),ncol=length(stopTimes)
        )
    } else {
        statIn <- dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)] # transform the simulation data to wald statistic for monitoring boundary
    }
    
    # pre-allocate vectors of quantities of interest; rows will be NA until the loop iteration for the trial stopping
    stopSign <- # will hold direction at trial stop (-1, 0, 1)
        stopEarly <- # will hold if the trial stopped early
        stopSize <- rep(NA,nrow(statIn))  # will hold number of blocks at trial stop
    
    # loop over each stopping time
    for (j in 1:length(stopTimes)){ 
        tmp <- statIn[,1:j] # statistics until current block/possible stopping time j
        
        # check if stopped at jth analysis
        if (j==1) {
            tmpStop <- (tmp < boundIn[j,1] | tmp > boundIn[j,2])
        } else {
            tmpStop <- apply(X=tmp,MARGIN=1,FUN=logReached,
                             lower=boundIn[1:j,1],upper=boundIn[1:j,2])
        }
        
        # modify stopEarly if j < J and crossed boundary, modify stopSign when stopEarly with direction or with 0 or direction when stopEarly==0
        if (j < length(stopTimes)) {
            stopEarly[is.na(stopEarly) & tmpStop==1] <- TRUE
            
            stopSign[is.na(stopSign) & tmpStop==1] <- sign(statIn[is.na(stopSign) & tmpStop==1,j])
        } else if (j == length(stopTimes)){
            stopEarly[is.na(stopEarly)] <- FALSE
            
            stopSign[is.na(stopSign)] <- ifelse(tmpStop[is.na(stopSign)]==0, # if not stopped yet, check if not crossed the boundary
                                                0, # conclude no treatment if no boundary crossing at the last look
                                                sign(statIn[is.na(stopSign),j])) # if boundary crossed, use sign
        }
        
        # modify stopSize if trial stopped early due to boundary or last analysis (only modify for those newly stopped)
        stopSize[is.na(stopSize) & 
                     ((tmpStop==1 & j < length(stopTimes)) | 
                          j==length(stopTimes))] <- stopTimes[j]
    }
    
    ### Calculate probabilities/average
    correctProb <- mean(stopSign == sign(truEffect))
    earlyProb <- mean(stopEarly)
    asn <- mean(stopSize)
    
    asnCI <- quantile(stopSize,probs=c(0.025,0.975))
    asnSD <- sd(stopSize)
    
    return(c(correctProb,earlyProb,asn,asnCI,asnSD))
}








