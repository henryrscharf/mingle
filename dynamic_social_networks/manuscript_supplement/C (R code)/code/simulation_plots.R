###############################################################################
###############################################################################
## Title: Diagnostic plots for "Dynamic social networks based on movement"
##
## Author: Henry Scharf
##
## Description: This script makes diagnostic plots for MCMC output.
###############################################################################
###############################################################################
## plot colors, linetypes, and characters ----
library(scales)              ## for plotting with transparency
truth.col <- "black"         ## adjust colors, line types, etc. as desired
truth.lty <- 2
truth.jitter <- 0.07         ## this makes the lines decipherable near 0 and 1
fit.col <- "black"
fit.lty <- 1
fit.lwd <- 1.7
fit.pch <- 1
naive.col <- "black"
naive.lty <- 3
naive.pch <- 2
naive.jitter <- 0.035

## load MCMC Rdata file by uncommenting below ----
## The working directory is assumed to be "code/"

# load("../data/simulation_mcmc_output_2016_07_21_08:01:19 PM.RData")
# list2env(results, envir = environment())

N <- length(output$phi)      ## map values from loaded data to more convenient names
TIME <- dim(output$W[[1]])[3]
n.indiv <- dim(results$output$W[[1]])[1]
nC2 <- n.indiv*(n.indiv - 1)/2
mu <- input$data$mu
W.sim <- input$data$W.sim
p1 <- input$params$p1

## proximity-based analysis ----
## get distances between all individuals at all time points
distances <- apply(mu, 2, function(x) dist(x))
## explore a range of radii
R <- seq(min(distances) + 1, max(distances), l = 1E3)
## build estimated network
p.R <- NULL
for (r in R) {
    w.index <- apply(distances, 2, function(x) {
        dist.mat <- matrix(NA, n.indiv, n.indiv)
        dist.mat[lower.tri(dist.mat)] <- x
        which(dist.mat < r, arr.ind = TRUE)})
    w.proxim.t <- array(0, dim = c(n.indiv, n.indiv, TIME))
    for(t in 1:TIME) {
        w.proxim.t[, , t][w.index[[t]]] <- 1
        w.proxim.t[, , t] <- w.proxim.t[, , t] + t(w.proxim.t[, , t])
    }
    p.R <- c(p.R, sum(w.proxim.t)/prod(dim(w.proxim.t)))
}
## pick the value of R that generates network with known p1
R.ideal <- R[max(which(p.R <= p1))]
## uncomment the lines below to make a plot of R vs. p1
# plot(R, p.R, type = "l", lwd = 2)
# points(R.ideal, input$params$p1, lwd = 2, cex = 2)

## generate the network corresponding to R.ideal
w.index <- apply(distances, 2, function(x) {
    dist.mat <- matrix(NA, n.indiv, n.indiv)
    dist.mat[lower.tri(dist.mat)] <- x
    which(dist.mat < R.ideal, arr.ind = TRUE)})
w.proxim.t <- array(0, dim = c(n.indiv, n.indiv, TIME))
for(t in 1:TIME) {
    w.proxim.t[, , t][w.index[[t]]] <- 1
    w.proxim.t[, , t] <- w.proxim.t[, , t] + t(w.proxim.t[, , t])
}

## compute posterior statistics for model-based nework [SLOW] ----
burnt <- ceiling(N/2)        ## index of iterations to be discarded as used.iterations
thinned <- 5                 ## thinning gap
##
used.iterations <- ceiling((burnt/thinned):(N/thinned)*thinned)
## compute acceptance rates based on unthinned iterations
ar.alpha <- sum(output$alpha[burnt:N] != output$alpha[burnt:N - 1]) /
    (N - burnt)
ar.phi <- sum(output$phi[burnt:N] != output$phi[burnt:N - 1]) /
    (N - burnt)
ar.p1 <- sum(output$p1[burnt:N] != output$p1[burnt:N - 1]) /
    (N - burnt)
## compute network summary 
w <- array(0, dim = c(n.indiv, n.indiv, TIME, N))
for (n in 1:N) {
     w[, , , n] <- output$W[[n]]
}
w.mean.t <- apply(w[, , , used.iterations], c(1, 2, 3), mean)
w.sd.t <- sqrt(w.mean.t * (1 - w.mean.t))
w.mode.t <- w.mean.t > 0.5

################## ----
## plot network ##
################## ----
## filename ----
n.columns <- 2
n.pairs.shown <- 8
pdf(file = paste("../fig/simulation_W_", date, ".pdf", sep = ""),
    height = 8*(n.pairs.shown/(4*n.columns)), width = 4*n.columns)
## plot ----
index <- expand.grid(1:n.indiv, 1:n.indiv)
index <- index[index[, 1] < index[, 2], ]
if (n.pairs.shown < nC2){
  index <- index[sample(1:dim(index)[1], n.pairs.shown), ]
}
## choose subset of pairs from manuscript
index <- data.frame("Var1" = c(1, 2, 3, 3, 4, 1, 2, 1),
                    "Var2" = c(5, 3, 4, 6, 5, 4, 5, 6))
par(mfcol = c(ceiling(n.pairs.shown/n.columns), n.columns),
    mar = c(2, 2.2, 2.6, 0.6), bg = "white")
## first plot
plot(1:length(w.mode.t[index[1, 1], index[1, 2], ]),
     w.mode.t[index[1, 1], index[1, 2], ],
     type = "l", lwd = 2, lty = fit.lty,
     main = paste(index[1, 1], index[1, 2], sep = "-"), cex.main = 1.5,
     ylim = c(-0.2, 1.2),
     col = alpha(fit.col, alpha = 0),
     yaxt = "n")
axis(2, at = c(0, 1), labels = c("no (0)", "yes (1)"), cex.axis = 1.5)
polygon(c(1:length(w.mode.t[index[1, 1], index[1, 2], ]),
          length(w.mode.t[index[1, 1], index[1, 2], ]):1),
        sapply(c(w.mean.t[index[1, 1], index[1, 2], ] -
                   1.96*w.sd.t[index[1, 1], index[1, 2], ],
                 rev(w.mean.t[index[1, 1], index[1, 2], ] +
                       1.96*w.sd.t[index[1, 1], index[1, 2], ])),
               function(x){min(max(0, x), 1)}),
        col = "gray87", border = NA)
lines(1:length(w.mode.t[index[1, 1], index[1, 2], ]),
      w.mean.t[index[1, 1], index[1, 2], ], type = "l",
      col = alpha(fit.col, alpha = 0.9), lwd = fit.lwd, lty = fit.lty, pch = fit.pch)
lines(1:length(w.mode.t[index[1, 1], index[1, 2], ]),
      (1 + 2 * truth.jitter) * W.sim[index[1, 1], index[1, 2], ] - truth.jitter,
      col = alpha(truth.col, alpha = 0.9), lwd = 2, lty = truth.lty)
lines(1:length(w.proxim.t[index[1, 1], index[1, 2], ]),
      (1 + 2 * naive.jitter) * w.proxim.t[index[1, 1], index[1, 2], ] - naive.jitter, 
      type = "l", col = alpha(naive.col, alpha = 0.9), lwd = 2, 
      pch = naive.pch, lty = naive.lty)
abline(h = 0.5, col = alpha("black", 0.25))
## legend
legend("topleft", lwd = 2, lty = c(fit.lty, truth.lty, naive.lty),
       col = c(fit.col, truth.col, naive.col),
       legend = c("model", "truth", "proximity-based"), bty = "n")
## rest of plots
for (ind in 2:dim(index)[1]) {
    i <- index[ind, 1]
    j <- index[ind, 2]
    plot(1:length(w.mode.t[i, j, ]), w.mode.t[i, j, ],
         type = "l", lwd = 2,
         main = paste(i, j, sep = "-"), cex.main = 1.5,
         ylim = c(-0.2, 1.2),
         col = alpha(fit.col, alpha = 0),
         yaxt = "n", xaxt = "n")
    polygon(c(1:length(w.mode.t[i, j, ]), length(w.mode.t[i, j, ]):1),
            c(w.mean.t[i, j, ] - 1.96*w.sd.t[i, j, ],
              rev(w.mean.t[i, j, ] + 1.96*w.sd.t[i, j, ])),
            col = "gray87", border = NA)
    lines(1:length(w.mode.t[i, j, ]), w.mean.t[i, j, ],
          col = alpha(fit.col, alpha = 0.9), lwd = fit.lwd, lty = fit.lty)
    lines(1:length(w.mode.t[i, j, ]), 
          (1 + 2 * truth.jitter) * W.sim[i, j, ] - truth.jitter,
          col = alpha(truth.col, alpha = 0.9), lwd = 2, lty = truth.lty)
    lines(1:length(w.proxim.t[i, j, ]), 
          (1 + 2 * naive.jitter) * w.proxim.t[i, j, ] - naive.jitter,
          col = alpha(naive.col, alpha = 0.9), lwd = 2, lty = naive.lty)
    abline(h = 0.5, col = alpha("black", 0.25))
}
## dev.off ----
dev.off()

############################## ----
## trace plots + histograms ##
############################## ----
## filename ----
pdf(file = paste("../fig/simulation_trace_", date, ".pdf", sep = ""),
    height = 6.5, width = 6.5)
## plot ----
par(mfcol = c(6, 2), oma = c(0, 0, 0, 0), mar = c(2, 3, 4, 1))
## beta
plot(output$beta[used.iterations], type = "l")
## alpha
plot(output$alpha[used.iterations], type = "l",
     main = paste("acceptance rate = ", ar.alpha))
## p1
plot(output$p1[used.iterations], type = "l",
     main = paste("acceptance rate = ", ar.p1))
## phi
plot(output$phi[used.iterations], type = "l",
     main = paste("acceptance rate = ", ar.phi))
## sigsq
plot(output$sigsq[used.iterations], type = "l")
## c
plot(output$c[used.iterations], type = "l")
## histograms
## beta
hist(output$beta[used.iterations], breaks = 50, probability = TRUE, 
     main = expression(beta))
abline(v = input$params$beta, col = alpha(truth.col, alpha = 0.8), lwd = 2)
## alpha
hist(output$alpha[used.iterations], breaks = 50, probability = TRUE, 
     main = expression(alpha))
dalpha <- function(xx, alpha = 1, beta = 1) {
    dx <- xx[-1] - xx[-length(xx)]
    density <- (1 + xx)^(alpha - 1)*(1 - xx)^(beta - 1)
    density <- density / (sum(density[-1] * dx))
    return(density)
}
xx <- seq(-1, 1, 1E-3)
lines(xx, dalpha(xx, results$priors$alpha.alpha, results$priors$beta.alpha), lty = 2)
abline(v = input$params$alpha, col = alpha(truth.col, alpha = 0.8), lwd = 2)
## p1
hist(output$p1[used.iterations], breaks = 50, probability = TRUE, 
     main = expression(p[1]))
abline(v = input$params$p1, col = alpha(truth.col, alpha = 0.8), lwd = 2)
xx <- seq(0, 1, length.out = 30)
lines(xx, dbeta(xx, results$priors$alpha.p1, results$priors$beta.p1), lty = 2)
## phi
hist(output$phi[used.iterations], breaks = 50, probability = TRUE, 
     xlim = range(c(quantile(output$phi, probs = c(0.001, 0.991)), 1)),
     main = expression(phi))
abline(v = input$params$phi, col = alpha(truth.col, alpha = 0.8), lwd = 2)
xx <- seq(0, 1, length.out = 30)
lines(xx, dbeta(xx, results$priors$alpha.phi, results$priors$beta.phi), lty = 2)
legend("topleft", lty = 2, legend = c("prior"), bty = "n")
## sigsq
hist(output$sigsq[used.iterations],
     breaks = 50, probability = TRUE, main = expression(sigma^2))
abline(v = input$params$sigsq, col = alpha(truth.col, alpha = 0.8), lwd = 2)
xx <- seq(par("usr")[1], par("usr")[2], length.out = 1000)
lines(xx, 1/dgamma(xx, results$priors$a.sigsq, results$priors$b.sigsq), lty = 2)
## c
hist(output$c[used.iterations], breaks = 50, probability = TRUE, main = expression(c))
abline(v = input$params$c, col = alpha(truth.col, alpha = 0.8), lwd = 2)
lines(xx, 1/dgamma(xx, results$priors$a.c, results$priors$b.c), lty = 2)
## dev.off ----
dev.off()

######################### ----
## parameters CI table ##
######################### ----
CI <- as.data.frame(lapply(output[1:6],
                           FUN = function(x) quantile(x,
                               probs = c(0.025, 0.5, 0.975))))
CI <- rbind(CI, "truth" = t(as.matrix(results$input$params[c(1:6)])))
t(CI)

################################### ----
## correlation between variables ##
################################### ----
## filename ----
pdf(file = paste("../fig/simulation_parameter_correlation_", date,
                 ".pdf", sep = ""), width = 6.6, height = 6.6)
## plot ----
require(corrplot)
param.corr <- cor(data.frame("p1" = output$p1,
                             "phi" = output$phi,
                             "sigsq" = output$sigsq, "c" = output$c,
                             "alpha" = output$alpha, "beta" = output$beta))
corrplot.mixed(param.corr, upper = "ellipse")
## dev.off ----
dev.off()