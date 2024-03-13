# Nonparametric multiple-ouput center-outward quantile regression with R

This folder contains the file `CenterOutwardRegression.R`, with the main functions for reproducing the figures in the manuscript "Nonparametric multiple-output center-outward quantile regression".
The function `QuantileRegressionOT` implements the quantile regression method with $k$-nearest neighbors weights (see Appendix C.2 for details). For the sake of visualization, in the simulated examples the regressor is univariate and the response is bivariate. The function `QuantileRegressionOT` is written for this setup, but the code can be easily modified to deal with higher dimensional cases. The function `plotQuantileRegressionOT3D` gives 3D plots of the conditional center-outward quantile contours and rays estimated with `QuantileRegressionOT`.

The files with a name starting with `Clover_shaped_model` correspond to model (1.8) (Fig. 2, right, and Fig. 6). `Mildly_nonconvex_model.R` corresponds to model (1.7) (Fig. 2, left, and Fig. 5) and `Elliptic_model.R` to model (4.1) in the main manuscript (Fig. 3 and Fig. 4).
