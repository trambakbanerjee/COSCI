The R package is forthcoming. However, you can still use COSCI for your analysis. This folder has an R file that holds two functions.

### 1. cosci_is

This is the main function that ranks the p features in an n by p design matrix where n represents the sample size and p is the number of features. This function takes three inputs.

i)  dat - n by p data matrix

ii) min.alpha - the smallest threshold (typically set to 0)

iii)small.perturbation - a small positive number to remove ties (typicall set to 10^(-6))

2. cosci_is_select:
