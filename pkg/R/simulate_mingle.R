###############################################################################
###############################################################################
## Henry Scharf
##
## These are functions used to generate mu from the process model.
###############################################################################
###############################################################################

#' making transition matrix for network edge probabilities
#'
#' @param p1 the expected density of the network
#' @param phi a number between 0 and 1 indicated the stability of the network. 1 results in a fully static network, 0 results in no Markovian dependence at all.
#'
#' @return 2x2 matrix filled columnwise with \eqn{p_{0|0}, p_{1|0}, p_{0|1}, p_{1|1}}
#' @export
make.cond.M.mat <- function(p1 = 0.2, phi = 0.95) {
  p0g0 <- 1 - (1 - phi)*p1
  p1g0 <- (1 - phi)*p1
  p0g1 <- (1 - phi)*(1 - p1)
  p1g1 <- 1 - (1 - phi)*(1 - p1)
  return(matrix(c(p0g0, p1g0, p0g1, p1g1), 2, 2))
}
#' compute \eqn{\overline{\boldsymbol{\mu}}(t)} for a given network and set of positions \eqn{\boldsymbol{\mu}(t)}. Intended mostly for internal use in simulation.
#'
#' @param W matrix of dimension \code{c(n.indiv, n.indiv)}
#' @param mu matrix of dimension \code{c(n.indiv, 2)}
#'
#' @return matrix of dimension \code{c(n.indiv, 2)}
#' @export
make.mubar.sim <- function(W, mu){
  w.it.sums <- rowSums(diag(dim(W)[1]) + W)
  ans <- (diag(dim(W)[1]) + W)%*%mu / w.it.sums
  return(ans)
}
#' compute \eqn{\tilde{\boldsymbol{\mu}}(t)} for a given network and set of positions \eqn{\boldsymbol{\mu}(t)}. Intended mostly for internal use in simulation.
#'
#' @inheritParams make.mubar.sim
#'
#' @return matrix of dimension \code{c(n.indiv, 2)}
#' @export
#'
make.mutil.sim <- function(W, mu){
  ##
  ##
  mubar <- make.mubar.sim(W, mu)
  diff <- matrix(mubar - mu, ncol = 2)
  norm <- sqrt(diff[, 2]^2 + diff[, 1]^2)
  norm[norm == 0] <- 1
  mutil <- array(diff/norm, dim = dim(mu))
  return(mutil)
}
#' build the precision matrix for the distribution of \eqn{\boldsymbol{\mu}(t)|\boldsymbol{\mu}(t)\dots}. Intended mostly for internal use in simulation.
#'
#' @param alpha scalar. alignment parameter from \eqn{[0, 1)}
#' @param sigsq scalar. marginal residual uncertainty after accounting for social effects.
#' @param c scalar. effective ego-network size for unconnected individuals.
#' @param wn.t matrix of dimension \code{c(n.indiv, n.indiv)} indicating network at time \eqn{t}
#'
#' @return matrix of dimension \code{c(n.indiv, n.indiv)}
#' @export
#'
buildQ.t.sub <- function(alpha, sigsq, c, wn.t) {
  w.plus <- rowSums(wn.t)
  w.plus[w.plus == 0] <- c
  Q.t <- sigsq^{-1} * (diag(w.plus) - alpha * wn.t)
  return(Q.t)
}

#' simulate paths \eqn{\mu} for \code{n.indiv} individuals at \code{T} time points.
#'
#' @param n.indiv scalar. number of individuals in network
#' @param T number of time points
#' @param params list of model parameters used to generate data. List must contain elements \code{alpha, beta, sigsq, c, phi, p1}.
#' @param params0 list of parameters for initial position \eqn{\mu_0}. Currently only \eqn{\sigma^2_0}.
#' @param seed seed to use in \code{set.seed()} before simulating data.
#'
#' @return list of two elements. \code{mu} is an array of dimension \code{c(n.indiv, T, 2)}, and \code{W.sim} is an array of dimension \code{c(n.indiv, n.indiv, T)}.
#' @export
#' @examples  
#' simulation.default <- sim_mingle()
sim_mingle <- function(n.indiv = 6,
                       T = 100,
                       params = list("alpha" = 0.8,
                                     "beta" = 0.5,
                                     "sigsq" = 1,
                                     "c" = 1/3,
                                     "phi" = 0.95,
                                     "p1" = 0.2),
                       params0 = list("sigsq0" = 30^2),
                       seed = 1985) {
  set.seed(seed)
  ##
  ## params
  ##
  alpha <- params$alpha
  beta <- params$beta
  sigsq <- params$sigsq
  c <- params$c
  phi <- params$phi
  p1 <- params$p1
  sigsq0 <- params0$sigsq0
  nC2 <- n.indiv*(n.indiv - 1)/2
  ##
  ## simulate data
  ##
  ##
  ## starting location
  ##
  mu <- array(0, dim = c(n.indiv, T, 2)) #(n.indiv, time, lat/long)
  mutil <- mu
  mu[, 1, ] <- matrix(rnorm(n = 2 * n.indiv, sd = sqrt(sigsq0)), n.indiv, 2)
  wn <- array(0, dim = c(n.indiv, n.indiv, T))
  for (i in 1:(n.indiv - 1)) {
    for (j in (i + 1):n.indiv) {
      wn[i, j, 1] <- wn[j, i, 1] <- sample(0:1, 1, prob = c(1 - p1, p1))
    }
  }
  ##
  ## position at t
  ##
  cond.M.mat <- make.cond.M.mat(p1, phi)
  for (t in 1:(T - 1)) {
    for (i in 1:(n.indiv - 1)) {
      for (j in (i + 1):n.indiv) {
        if (wn[i, j, t] == 0) {
          p.ijt <- cond.M.mat[, 1]
        } else {
          p.ijt <- cond.M.mat[, 2]
        }
        wn[i, j, t + 1] <- wn[j, i, t + 1] <- sample(0:1, 1, prob = p.ijt)
      }
    }
    Q.tp1 <- kronecker(buildQ.t.sub(alpha, sigsq, c, wn[, , t + 1]), diag(2))
    mutil[, t, ] <- make.mutil.sim(wn[, , t], mu[, t, ])
    mu[, t + 1, ] <- matrix(mvrnorm(n = 1,
                                    mu = c(t(mu[, t, ] + beta * mutil[, t, ])),
                                    Sigma = solve(Q.tp1)), n.indiv, 2, byrow = T)
  }
  return(list("mu" = mu, "W.sim" = wn))
}
