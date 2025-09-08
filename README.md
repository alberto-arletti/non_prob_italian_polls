# non_prob_italian_polls
This repository hosts part of the code and the data relative to the working paper "Adjusting Selection Bias in Non-Probability Samples: the Case of Italian Electoral Polls". 
dataset_5.csv file contains the observations for the corresponding dataset for the study variables used for adjustment. 
The census folder contains the different Italian census tables that can be used for adjustment in .csv format. GT.csv contains the true election results for the 2022 national elections at a country-wide level. 
The helper.R file contains the main functions used for adjustment: Raking, neural network post-stratification, doubly robust post-stratification, multilevel regression and post-stratification, propensity score based inverse probability weighting
