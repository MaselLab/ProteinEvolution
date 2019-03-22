library(MASS)


# This is a short script that is used to determine the optimal lambda value for 
# a box-cox transform 

df <- read.csv(file = "",header = T)
bc <- boxcox(df$ISD~1, lambda = seq(.1,.7,0.01))
lambda <- bc$x[which.max(bc$y)]
lambda


