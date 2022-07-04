## Loading the model

This repository contains the R code (1) and associated data for the
schistosomiasis transmission dynamics model described in Díaz et
al. 2022 (2). The model requires installation of the deSolve (3) package
and can be loaded, with default parameter values, as follows:

    library(deSolve)
    source("TnTmodel.R")
    source("parms.R")

## Running the model

The model is best used by running simulations to equilibrium and running
interventions using the equilibrium (baseline) value for the mean number
of worms per host as the initial value. The below code illustrates this
approach, simulating an intervention of 5 annual treatments:

    base <- runmod(parameters)
    parameters["dotx"] <- 1
    parameters["start.tx"] <- 1
    parameters["stop.t"] <- 9
    parameters["n.tx"] <- 5
    parameters["freq.tx"] <- 1
    interv <- runmod(parameters, inits = base[nrow(base),2])

![](README_files/figure-markdown_strict/plotmod-1.png)

# Simulation data

All simulation results presented in Díaz et al. (2) are stored in the
.rds files and can be loaded as follows:

    baseres <- readRDS("baselinesims.rds")
    intervres <- readRDS("interventionsims.rds")
    cyres <- readRDS("caseyears.rds")

# References

<span class="csl-left-margin">1. </span><span class="csl-right-inline">R
Core Team. *R: A language and environment for statistical computing*.
Vienna, Austria: R Foundation for Statistical Computing (2021).
<https://www.R-project.org/></span>

<span class="csl-left-margin">2. </span><span
class="csl-right-inline">Díaz AV, Lambert S, M. Inês Neve an AB, Léger
E, Diouf ND, Sène M, Webster JP, Walker M. Modelling livestock
test-and-treat: A novel one health strategy to control schistosomiasis
and mitigate against drug resistance. *Front Trop Dis* (2022) in
press:</span>

<span class="csl-left-margin">3. </span><span
class="csl-right-inline">Soetaert K, Petzoldt T, Setzer RW. Solving
differential equations in R: Package deSolve. *J Stat Softw* (2010)
33:1–25. doi:
[10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09)</span>
