What is COSCI?
======

COSCI (COnvex Screening for Cluster Information) is a non-parametric method for ranking and screening non-informative features in large scale cluster analysis problems. Unlike the non-parametric density estimation based screening techniques, COSCI is very
scalable, and can successfully handle datasets with more than one million observations. COSCI discards non-informative features by first computing a clustering score for the clustering tree constructed for each feature, and then thresholds the resulting values.

Perfect Screening Property
---------

COSCI produces a score for each feature, which reflects its relative importance for clustering, and then screens out the features with lower scores. It enjoys a perfect screening property in the sense that under mild regularity conditions on the densities of the features, COSCI screens out all the non-informative features with high probability.

How to use this repository?
----------

This repository has two main folders: (1) R Package - this is forthcoming but you can still use the COSCI functions placed inside this folder for your analysis. (2) Tables and Figures in the paper - this folder holds the scripts that reproduce the analysis in the paper [1]. Send me an email if anything does not work as expected.

References
=======
[1.] [Feature Screening in Large Scale Cluster Analysis](http://www.sciencedirect.com/science/article/pii/S0047259X17300271)    
Banerjee, T., Mukherjee, G. and Radchenko P.  *Journal of Multivariate Analysis (2017)*

[2.] [Convex clustering via ℓ1 fusion penalization](10.1111/rssb.12226)   
Radchenko P., Mukherjee G.   *J. R. Stat. Soc. Ser. B Stat. Methodol. (2017)*

