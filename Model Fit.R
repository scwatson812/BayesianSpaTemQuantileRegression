###Load Packages
library(fields)
library(GIGrvg)
library(truncnorm)
library(Matrix)
library(MASS)
library(invgamma)
library(spam)
library(splines)
library(rgdal)
library(raster)
library(RColorBrewer)

###Load the Datafile from Github
path.data = 'path/to/data/'
load(paste0(path.data,'PManalysisUnstandardizedVars.RData'))
###Load the shapefile for the maps
us.states = readOGR(paste0(path.source,'states'))

###Specify the path to the functions
path.source = 'path/to/source/'
#Specify a file path to write the output will be written
path.out = "path/to/output/"

source(paste0(path.source,'mvnsamp.R'))
source(paste0(path.source,'colfunc.R'))
source(paste0(path.source,'rowfunc.R'))
source(paste0(path.source,'sparseW.R'))
source(paste0(path.source,'sparseWtime.R'))
source(paste0(path.source,'zsamp.R'))
source(paste0(path.source,'mcmcFitSummer.R'))
source(paste0(path.source,'mcmcFitWinter.R'))

fit.summer = mcmcFitSummer(20000,15000,0.5,path.out,us.states)
fit.winter = mcmcFitWinter(20000,15000,0.5,path.out,us.states)

