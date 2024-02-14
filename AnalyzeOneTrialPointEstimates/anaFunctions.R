##############################################
# Functions to Analyze Simulation Study
# Number of Blocks for Sequential Monitoring
# 
# Subodh Selukar
# 2021-05-18
##############################################

# -------------------------------- Description ------------------------------------------------ #

##### This simulation study is intended to describe 
##### early stopping with sequential monitoring boundaries 
##### for N-of-1 trials with 2 treatments

#### Goal is to evaluate different estimator choices 
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
# RCTdesign

# -------------------------------- Functions ------------------------------------------------ #

### Function to create monitoring boundaries (requires RCTdesign)
# only 3 specific boundary shapes supported currently

# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. totBlocks: total number of blocks (max "sample size")
# 3. boundType: type of boundary (supporting symmetric OBF, symmetric Pocock and asymmetric with OBF upper and Pocock lowerr)
# 4. alpha: one-sided alpha level (defaults to 0.025, not planning to change)

# Outputs: 
# seqDesign object
makeBounds <- function(stopTimes,
                       boundType="symOBF",
                       alpha=0.025){
    
    if (boundType=="symPoc"){
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         design.family="Pocock",
                         test.type="two.sided")
    } else if (boundType=="symOBF"){
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         design.family="OBF",
                         test.type="two.sided")
    } else if (boundType=="asym") { # asymmetric with OBF upper and Poc lower; NOTE: allows for conclusion of null at last analysis (if no crossing at last stage)
        upBound <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                             arms=1,
                             design.family="OBF",
                             test.type="two.sided")$boundary[,4]
        loBound <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                             arms=1,
                             design.family="Pocock",
                             test.type="two.sided")$boundary[,1]
        
        
        out <- seqDesign(nbr.analysis=length(stopTimes),sample.size=stopTimes,
                         arms=1,
                         exact.constraint=cbind(loBound,0,0,upBound),
                         test.type="two.sided")
    } else stop("Not supported")
    
    return(out)
}


### Function to convert from z-boundaries to (approximate) t-boundaries (per Nikolakopoulos 2016)
## Edited for seqDesign
# Inputs: 
# 1. stopTimes: vector of stopping times (i.e., blocks after which data are analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 2. zIn: input Z critical values, as made by makeBounds()

# Outputs: 
# matrix with nrow=number of stopping times and 2 columns with the monitoring boundaries at each stopping time
makeZtoTBounds <- function(stopTimes,
                           boundObj){
    
    zIn <- seqBoundary(boundObj,scale="Z")
    
    out <- cbind(
        qt(pnorm(zIn[,1]),stopTimes-2),
        NA,
        NA,
        qt(pnorm(zIn[,4]),stopTimes-2) # zIn is matrix with 4 columns
    )
    out[,1] <- ifelse(is.nan(out[,1]),-Inf,out[,1]) # need to deal with NaN
    out[,4] <- ifelse(is.nan(out[,4]),Inf,out[,4])
    
    out[nrow(out),2:3] <- out[nrow(out),c(1,4)] # need to specify the equivalence boundaries for the last analysis are = to the rejection boundaries
  
    desOut <- seqDesign(
        #boundObj,
                     exact.constraint=seqBoundary(out,scale="Z"),design.family="Z",#, # update the Z boundaries
                     sample.size=boundObj$specification$sample.size,
                     nbr.analyses = length(boundObj$specification$sample.size),
                     arms=1
                     # alt.hypothesis=hi$hypothesis$theta.alt,
                     # G=hi$parameters$G,
                     # power="calculate"
    )
    
    return(desOut)
}



### Function to return bias-adjusted estimates from a naive estimate at stopping (requires RCTdesign)
## Edited to remove CBAM
# Inputs: 
# 1. boundObj: the seqDesign object used for monitoring
# 2. stopIdx: index of stopping
# 3. est.Naive: naive estimate at stopping
# 4. se.Naive: naive estimate of standard error at stopping
# 5. h: increment size for numerical derivative

# Outputs: list of two vectors of length 6
# 1. vector of naive + bias-adjusted estimates
# 2. vector of corresponding estimated standard errors

calcEst <- function(boundObj,stopIdx,est.Naive,se.Naive,h){
    newBound <- update.seqDesign(boundObj,
                                 exact.constraint=seqBoundary(boundObj,scale="Z"),design.family="Z", # otherwise, update may change the Z boundaries! (esp for asymmetric bounds)
                                 sd=se.Naive*sqrt(boundObj$specification$sample.size[stopIdx]) # need to update the input standard deviation (note that seqDesign expects SD not SE)
                                 )
    
    ### RCTdesign natively returns all bias-adjusted estimates except conditional BAM
    infOut <- seqInference(newBound,observed=est.Naive,analysis.index = stopIdx,inScale="X") 
    estMUE.SM <- infOut[,4]
    estMUE.AT <- infOut[,5]
    estMUE.LR <- infOut[,6]
    est.BAM <- infOut[,7]
    
    ## 2023-05-01: CBAM was removed from calculations; some remnants left alone to allow backward compatibility
    
    ## collect all estimates
    
    estOut <- cbind(
        estMUE.SM,estMUE.AT,estMUE.LR,
        est.BAM,
        NA, #est.CBAM,
        est.Naive
    )
    
    ### Use Z-estimator large-sample theory to produce estimated variance/SE
    
    pSeqOut <- cbind(
        (pSeq(newBound,observed=estOut[,1],theta=estOut[,1]-h,analysis.index=stopIdx)[,4]-
             pSeq(newBound,observed=estOut[,1],theta=estOut[,1]+h,analysis.index=stopIdx)[,4])/(2*h), # MUE.SM
        (pSeq(newBound,observed=estOut[,2],theta=estOut[,2]-h,analysis.index=stopIdx)[5]-
             pSeq(newBound,observed=estOut[,2],theta=estOut[,2]+h,analysis.index=stopIdx)[,5])/(2*h), # MUE.AT
        (pSeq(newBound,observed=estOut[,3],theta=estOut[,3]-h,analysis.index=stopIdx)[,6]-
             pSeq(newBound,observed=estOut[,3],theta=estOut[,3]+h,analysis.index=stopIdx)[,6])/(2*h) # MUE.LR
    )
    
    meanSeqOut <- (meanSeq(newBound,estOut[,4]+h)[2]-meanSeq(newBound,estOut[,4]-h)[2])/(2*h)
    
    B <- cbind(
        pSeqOut,
        meanSeqOut[[1]],
        NA # leftover from CBAM
    )
    
    A <- (cbind( # delta method 
        (pSeq(newBound,observed=estOut[,1]-h,theta=estOut[,1],analysis.index=stopIdx)[,4]-
             pSeq(newBound,observed=estOut[,1]+h,theta=estOut[,1],analysis.index=stopIdx)[,4])/(2*h), # MUE.SM
        (pSeq(newBound,observed=estOut[,2]-h,theta=estOut[,2],analysis.index=stopIdx)[5]-
             pSeq(newBound,observed=estOut[,2]+h,theta=estOut[,2],analysis.index=stopIdx)[,5])/(2*h), # MUE.AT
        (pSeq(newBound,observed=estOut[,3]-h,theta=estOut[,3],analysis.index=stopIdx)[,6]-
             pSeq(newBound,observed=estOut[,3]+h,theta=estOut[,3],analysis.index=stopIdx)[,6])/(2*h), # MUE.LR
        
        1, # no need for Delta method for BAM or CBAM because derivative is 1
        NA #1 # leftover from CBAM
    )**2)*
        ((se.Naive**2)) # NOTE: naive SE takes into account sample size
    
    seOut <- cbind(
        sqrt( # SE instead of variance
            (A/(B**2)) # large-sample approximation for variance of Z-estimators (NOTE: naive SE already took into account sample size)
        ),
        se.Naive
    )
    
    return(list(
        estOut,seOut
    ))
    
}

### Function to use input monitoring boundary and simulation results to output bias-adjusted estimates with standard errors

## Inputs:
# 1. dataIn: simulation results - a matrix of the estimates and standard errors with 
#   nrow() = number of simulations
#   ncol() = 2*(number of trial blocks-1) with 
#   odd columns the estimates after each block (no analysis after block 1) and 
#   even columns the (observed) standard error
# 2. stopTimes: vector of stopping times (i.e., blocks after which data would have been analyzed)
#   stopTimes[1] >= 2 and stopTimes[length(stopTimes)] == last block of trial
# 3. boundIn: seqDesign object
# 4. h: increment size for numerical derivative

## Outputs: list of 2 matrices (numSims rows and 6 columns, one for each estimator)
# 1. estimates 
# 2. estimated standard errors

## Important: block 1 estimate/SE not included in simulation results - impacts the indexing

anaSims <- function(dataIn,stopTimes,boundIn,
                    h=1e-6){
    
    numBlocks <- stopTimes[length(stopTimes)]
    
    ### First find the block the trial stops at
    bounds <- matrix(
        seqBoundary(boundIn,scale="Z")[,c(1,4)],
        ncol=2)
    stopFun <- function(x) min(which(x < bounds[,1] | x > bounds[,2]),nrow(bounds)) # return logical index of stop
    
    if (length(stopTimes) > 1){
        statIn <- dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)] # transform the simulation data to wald statistic for monitoring boundary
          
    } else {
        statIn <- matrix(
            dataIn[,2*(stopTimes-1)-1]/dataIn[,2*(stopTimes-1)], # transform the simulation data to wald statistic for monitoring boundary
            ncol=length(stopTimes)
        )
    }
    
    
        
    stopIdx <- apply(statIn,1, # look at each row's Wald statistics
                       FUN = stopFun # find the index of stopping for that trial
                       )
    
    stopBlock <- stopTimes[stopIdx] # return the stopping block
    
    ### Find the corresponding naive estimate
    
    est.Naive <-  dataIn[,2*(stopTimes-1)-1][cbind(1:nrow(dataIn),stopIdx)] # this picks out the specific stopIdx column of each row to get the resultant naive estimates
    se.Naive <- dataIn[,2*(stopTimes-1)][cbind(1:nrow(dataIn),stopIdx)] # this repeats to return the naive SE
    
    ### Loop over each trial to return the desired estimates/standard errors 
    # Need to run this loop because each bias-adjusted estimate requires an update to the trial's SD
    
    estOut <- seOut <- matrix(ncol=6,
                              nrow=nrow(dataIn))
    
    for (i in 1:nrow(dataIn)){
        tmp <- tryCatch(calcEst(boundIn,stopIdx[i],est.Naive[i],se.Naive[i],h),error=function(e) return(NA))
        
        if (length(tmp)==2){ # some minor error handling if a given design did not work
            estOut[i,] <- tmp[[1]]
            seOut[i,] <- tmp[[2]]
        } 
        
    }
    
    ### Report out estimates with estimated standard errors
    
    return(list(
        estOut,
        seOut
    ))
}
