The R package is forthcoming. However, you can still use COSCI for your analysis. This folder has an R file that holds two functions.

### 1. cosci_is

This is the main function that ranks the p features in an n by p design matrix where n represents the sample size and p is the number of features. This function takes three inputs.

i)   dat - n by p data matrix

ii)  min.alpha - the smallest threshold (typically set to 0)

iii) small.perturbation - a small positive number to remove ties (typicall set to 10^(-6))

The output will be a p vector of scores. 

### 2. cosci_is_select

Once you have the feature scores from (1), you can

(i.) select the features based on a pre-defined threshold, 

(ii.) use table 1 in our paper to determine an appropriate threshold or 

(iii.) use the data driven approach described in our paper to select the features and obtain an implicit threshold value. 

The function 'cosci_is_select' implements option (iii.). This function takes two inputs (a) the p vector of scores and (b) gamma.
Please refer to our paper for what gamma stands for in this context. If your sample size n is smaller than 100, setting gamma = 0.85 is recommended. 
