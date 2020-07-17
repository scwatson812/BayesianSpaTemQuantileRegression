# BayesianSpaTemQuantileRegression
Contains the code used for the analysis presented in 'Identifying Meteorological Drivers of PM2.5 Levels via a Bayesian Spatial Quantile Regression
The file ‘PManalysisUnstandardizedVars.RData’ contains the PM2.5 data and the meteorological factors.
The ‘states.xxx’ files are the GIS file for creating the covariate effect maps.
The file ‘Model Fit.R’ performs model fitting for the summer and winter data.
The main functions are ‘mcmcFitSummer.R’ and ‘mcmcFitWinter.R’, which perform model fitting for the summer and winter data respectively.  
Each function takes 5 inputs: the number of MCMC iterations (G), the number of MCMC iterations to discard as burn in (burn), the quantile of interest (r), a string specifying the directory to which the map images will be written (path) and the shapefile of US states used for making the map images (m1).  

The other functions are simply helper functions used to fit the model.
