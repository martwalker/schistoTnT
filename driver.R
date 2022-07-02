rm(list = ls())
##########################
## load packages
###########################
library(deSolve)

####################################
## load model & associated functions
####################################
source("TnTmodel.R")

####################################
## load default parameter values
####################################
source("parms.R")

#######################################
## run model to baseline 
#######################################
base <- runmod(parameters)

#######################################
## run an intervention from baseline
#######################################
parameters["dotx"] <- 1
parameters["start.tx"] <- 1
parameters["stop.t"] <- 9
parameters["n.tx"] <- 5

interv <- runmod(parameters, inits = base[nrow(base),2])

plot(interv$prev~interv$time, type="l", ylim=c(0,1), ylab="prevalence", xlab="time")

