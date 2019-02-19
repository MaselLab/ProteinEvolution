library(MASS)
library(nlme)
library(lme4)
library(MuMIn)

# Author : Sara Willis
# Date   : February 18, 2019


# This is a short R script that is used to determine the optimal lambda value for 
# a box-cox transform 

df <- read.csv(file = "",header = T)
bc <- boxcox(df$ISD~1, lambda = seq(.1,.7,0.01))
lambda <- bc$x[which.max(bc$y)]
lambda


