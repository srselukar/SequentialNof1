This folder contains the simulation and analysis functions used to conduct the simulation study described in the article "A Framework for Sequential Monitoring of Individual N-of-1 Trials and Combining Results across a Series of Sequentially-Monitored N-of-1 Trials"

Each subfolder is dedicated to the functions for a component of the study

GenerateNof1Trials: simFunctions.R
- Main function: runSimSet(), returns the point estimates and variances at every analysis block from independent N-of-1 trials under a user-specified setting
- Helper functions generate the trial data and compute the estimates

AnalyzeType1Error: anaFunctions.R
- Main function: anaSims(), summarizes the stopping probability at each potential stopping time and the overall stopping probability; under the null, these are stagewise and overall type-1 error probabilities; user specifies the boundary shape and stopping times and provides the output from one run of runSimSet() above
- Helper functions construct sequential boundaries, convert Z-quantiles to T-quantiles   and test if a monitoring boundary is crossed

AnalyzeSeqProperties: anaFunctions.R
- Main function: anaSims(), summarizes the power, proportion of trials stopped before planned end, average block number at stopping (along with the SD and 2.5%-97.5% quantiles of block number at stopping); user specifies the boundary shape, stopping times and true effect size and provides the output from one run of runSimSet() above
- Helper functions construct sequential boundaries (note that some boundary options were dropped after the type-1 error analyses), convert Z-quantiles to T-quantiles   and test if a monitoring boundary is crossed

AnalyzeOneTrialPointEstimates: anaFunctions.R
- Main function: anaSims(), calculates the naive and bias-adjusted point estimates for  an input from runSimSet() 
- Helper functions construct sequential boundaries (note that some boundary options were dropped after the type-1 error analyses), convert Z-quantiles to T-quantiles   and test if a monitoring boundary is crossed
- Note that point estimation is computationally intensive; the RCTdesign package is required and not available on CRAN, but is available by request at http://www.rctdesign.org/Software.html

AnalyzeSeries: anaFunctions.R
- Function: simSeries(), computes the combined point estimates from series of trials; uses an input dataset of input point estimates from independent trials and groups according to series size 
- This function expects the output from anaSims() in the script AnalyzeOneTrialPointEstimates/anaFunctions.R; that data structure is required (though would be easy to modify)

