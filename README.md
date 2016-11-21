## What is COSCI?

COSCI (COnvex Screening for Cluster Information) is a non-parametric method for ranking and screening non-informative features in large scale cluster analysis problems. Unlike the non-parametric density estimation based screening techniques, COSCI is very
scalable, and can successfully handle datasets with more than one million observations. COSCI discards non-informative features by first computing a clustering score for the clustering tree constructed for each feature, and then thresholds the resulting values.

### Perfect Screening Property

COSCI produces a score for each feature, which reflects its relative importance for clustering, and then screens out the features with lower scores. It enjoys a perfect screening property in the sense that under mild regularity conditions on the densities of the features, COSCI screens out all the non-informative features with high probability.

### Reference
[1.] Banerjee, T., Mukherjee, G. and Radchenko, P. “Feature Screening in Large Scale Cluster Analysis.” 

## How to use this repository?

This reporsitory has two main folders: (1) R Package - this is forthcoming but you can still use the COSCI functions placed inside this folder for your analysis. (2) Tables and Figures in the paper - this folder holds the scripts that reproduce the analysis in the paper [1].
