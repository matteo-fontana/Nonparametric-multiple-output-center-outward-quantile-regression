# Nonparametric multiple-ouput center-outward quantile regression with R

This folder contains the file 'CenterOutwardRegression.R', with the main functions for reproducing the figures in the manuscript "Nonparametric multiple-output center-outward quantile regression".
The function 'QuantileRegressionOT' implements the quantile regression method with $k$-nearest neighbors weights (see Appendix C.2 for details). For the sake of visualization, in the simulated examples the regressor is univariate and the response is bivariate. The function 'QuantileRegressionOT' is written for ths setup, but the code can be easily modified to deal with higher dimensional cases.

 The files with a name starting with 'Clover_shaped_model' correspond to model (1.8), 'Mildly_nonconvex_model.R' correspond to model (1.7) and 
