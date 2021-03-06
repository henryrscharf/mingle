###############################################################################
###############################################################################
## Title: Simulation for "Dynamic social networks based on movement"
##
## Author: Henry Scharf
##
## Description: This script fits the latent social network model for movement 
## to simulated data.
###############################################################################
###############################################################################
## install mingle ----
## Our package contains functions used to simulate the data, and sample from 
## the posterior. Uncomment the following to install.

# devtools::install_github(repo = "henryrscharf/mingle/pkg")

## preliminaries ----
library(mingle)              
N <- 2E5                     ## number of iterations in the MCMC chain
n.indiv <- 6                 ## number of simulated individuals
TIME <- 100                  ## number of time steps
seed <- 1985                 ## seed used for reproducibility at various stages

## model parameters + simulated data ----
params <- list("alpha" = 0.9, "beta" = 0.5,
               "sigsq" = 1, "c" = 1/3,
               "phi" = 0.95, "p1" = 0.2)
params0 <- list("sigsq0" = 30^2)       ## variation in initial position
data <- sim_mingle(n.indiv = n.indiv,  ## simulate the data
                   T = TIME, params, 
                   params0, 
                   seed = seed)
input <- list("params" = params,       ## create list of inputs
              "params0" = params0,
              "data" = data)

## hyperparameters + tuning ----
## p1
alpha.p1 <- 1
beta.p1 <- 1

## phi
alpha.phi <- 100             ## mean of 0.9
beta.phi <- 100/9            ## approx 95% of mass between 0.85 and 0.94

## sigma
a <- 1E-1
b <- 1E-3

## c
b.c <- 1.5                   ## approx 89% of mass between 0 and 1
a.c <- 3*b.c - 1             ## mode of 1/3

## alpha
alpha.alpha <- 1
beta.alpha <- 1

## beta
sigsq.beta <- 1000
mu.beta <- 0

## tuning
sigsq.alpha.tune <- 0.1^2    ## smaller -> higher AR
beta.phi.tune <- 7           ## larger -> higher AR
beta.p1.tune <- 24           ## larger -> higher AR

priors <- list("alpha.p1" = alpha.p1,
               "beta.p1" = beta.p1,
               "beta.p1.tune" = beta.p1.tune,
               "alpha.phi" = alpha.phi,
               "beta.phi" = beta.phi,
               "beta.phi.tune" = beta.phi.tune,
               "alpha.alpha" = alpha.alpha,
               "beta.alpha" = beta.alpha,
               "sigsq.alpha.tune" = sigsq.alpha.tune,
               "mu.beta" = mu.beta,
               "sigsq.beta" = sigsq.beta,
               "a.sigsq" = a,
               "b.sigsq" = b,
               "a.c" = a.c,
               "b.c" = b.c)

## random starting values + storage ----
set.seed(seed)
## network
W <- array(0, dim = c(n.indiv, n.indiv, TIME, N))
for (tt in 1:TIME){
  W[, , tt, 1][upper.tri(W[, , tt, 1])] <-
    sample(1:0, size = n.indiv * (n.indiv - 1) / 2,
           replace = TRUE,
           prob = c(params$p1, 1 - params$p1)) #(n.indiv, n.indiv, T, N)
  W[, , tt, 1][lower.tri(W[, , tt, 1])] <-
    t(W[, , tt, 1])[lower.tri(W[, , tt, 1])]
}
W[, , , 2] <- W[, , , 1]
## params
alpha <- rep(runif(1, min = -1, max = 1), N)
beta <- rep(rnorm(1, sd = 10^1), N)
sigsq <- rep(1/rgamma(1, 1E-1, 1E-3), N)
p1 <- rep(runif(1, min = 0.25, max = 0.75), N)
c <- rep(runif(1), N)
phi <- rep(runif(1, min = 0.85, max = 0.94), N)          

## Fixed parameters can be specified by defining the first value as NA, and the 
## second value as the used value. For example:
# phi[1] <- NA

## fit_mingle() is built for multiple imputation with the fourth dimension of 
## the data array indexing K, so set K = 1.
if (length(dim(data$mu)) == 3)
  mu.aug <- array(data$mu, c(dim(data$mu), 2)); data$mu <- mu.aug

## mcmc ----
## WARNING--SLOW. This takes about ~ 0.3 seconds/iteration, or ~17 hours for 2e5 iterations and can be memory-intensive ~14GB (this storage demand could be reduced in the future with more careful coding.)]
system.time({
  output <- fit_mingle(N_iterations = N, mu = data$mu, 
                       phi, alpha, beta, p1, sigsq, c, W,
                       alpha.p1, beta.p1, alpha.phi, beta.phi,
                       alpha.alpha, beta.alpha,
                       mu.beta, sigsq.beta, a, b, a.c, b.c,
                       beta.p1.tune, beta.phi.tune, sigsq.alpha.tune)
})

## save results ----
## You'll likely want to run this as a BATCH job, so this saves the output.
date <- format(Sys.time(), "%Y_%m_%d_%X")
results <- list("input" = input, "priors" = priors, "output" = output,
                "seed" = seed, "date" = date)
save(results, file = paste("../data/simulation_mcmc_output_",
                           date, ".RData", sep = ""))

## plots ----
source("simulation_plots.R")
