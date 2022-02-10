# jonashaslbeck@gmail.com; Jan 14, 2022

# ----------------------------------------------------------------------
# ----- What is Happening Here? ----------------------------------------
# ----------------------------------------------------------------------

# Show in two additional analyses that pairwise KLD does not fully
# determine the recovery performance/probability of K (this is for
# the appendix)

# ----------------------------------------------------------------------
# ----- Load packages & source aux functions ---------------------------
# ----------------------------------------------------------------------

library(RColorBrewer)
library(mclust)
library(MASS)

# Functions for data generation and estimation
source("aux_functions.R")


# ----------------------------------------------------------------------
# ----- Case 1: Differences due to Means vs. Covariances ---------------
# ----------------------------------------------------------------------

N <- 1000

# --- GMM 1: mean differences ---
p <- 2
Sigma <- diag(p)
mean1 <- c(0,2)
mean2 <- c(2,0)
data1 <- mvrnorm(n = N, mu = mean1, Sigma = Sigma)
data2 <- mvrnorm(n = N, mu = mean2, Sigma = Sigma)
data_GMM1 <- rbind(data1, data2)
GMM1 <- list("means" = list(mean1, mean2),
             "sigma" = list(Sigma, Sigma))

# Calc KLD
KLD(mu1 = GMM1[[1]][[1]],
    mu2 = GMM1[[1]][[2]],
    sig1 = GMM1[[2]][[1]],
    sig2 = GMM1[[2]][[2]])

# --- GMM 2: Cov differences ---
p <- 2
mean1 <- mean2 <- c(0,0)
Sigma1 <- Sigma2 <- diag(p)
co <- 0.8165
Sigma1[1,2] <- Sigma1[2,1] <- co
Sigma2[1,2] <- Sigma2[2,1] <- -co
data1 <- mvrnorm(n = N, mu = mean1, Sigma = Sigma1)
data2 <- mvrnorm(n = N, mu = mean2, Sigma = Sigma2)
data_GMM2 <- rbind(data1, data2)
GMM2 <- list("means" = list(mean1, mean2),
             "sigma" = list(Sigma1, Sigma2))

# Calc KLD
KLD(mu1 = GMM2[[1]][[1]],
    mu2 = GMM2[[1]][[2]],
    sig1 = GMM2[[2]][[1]],
    sig2 = GMM2[[2]][[2]])


# --- Make Picture ---
sc <- .7
pdf("figures/Figure_APP_KLDnotall_1.pdf", width = 10*sc, height = 5*sc)
par(mfrow=c(1,2), mar=c(3,3,2,1))
PlotCont(GMM1[[1]], GMM1[[2]])
title(main="Difference in Means", font.main = 1, cex.main=.9)
PlotCont(GMM2[[1]], GMM2[[2]])
title(main="Difference in Covariance", font.main = 1, cex.main=.9)
dev.off()

# --- Small simulation determining performance difference ---

nIter <- 100
N <- 1000
m_out <- matrix(NA, nIter, 2)

set.seed(1)

for(i in 1:nIter) {

  # GMM1
  data <- DataGen_GMM2(mu = GMM1[[1]],
                       Sigma = GMM1[[2]],
                       cats = "cont",
                       pi = c(.5, .5),
                       N = N)

  out1 <- Est_GGM(data = data$data, Kseq = 1:7)
  m_out[i, 1] <- out1$K

  # GMM1
  data <- DataGen_GMM2(mu = GMM2[[1]],
                       Sigma = GMM2[[2]],
                       cats = "cont",
                       pi = c(.5, .5),
                       N = N)

  out2 <- Est_GGM(data = data$data, Kseq = 1:7)
  m_out[i, 2] <- out2$K

  print(i)

} # end for: iter

apply(m_out==2, 2, mean)


# ----------------------------------------------------------------------
# ----- Case 2: Vary Distribution of Differences in Means --------------
# ----------------------------------------------------------------------

p <- 10

# --- GMM 1: All in 1 dimension ---
mu <- list(c(2.8281, rep(0,p-1)),
           c(0, rep(0,p-1)))
sig <- list(diag(p),
            diag(p))
GMM1 <- list(mu, sig)

KLD(mu1 = GMM1[[1]][[1]],
    mu2 = GMM1[[1]][[2]],
    sig1 = GMM1[[2]][[1]],
    sig2 = GMM1[[2]][[2]])


# --- GMM 2: Equal across p dimensions ---
theta <- .8945 # p=10
mu <- list(rep(theta, p),
           rep(0,p))
sig <- list(diag(p),
            diag(p))
GMM2 <- list(mu, sig)

KLD(mu1 = GMM2[[1]][[1]],
    mu2 = GMM2[[1]][[2]],
    sig1 = GMM2[[2]][[1]],
    sig2 = GMM2[[2]][[2]])


# --- Small simulation determining performance difference ---

nIter <- 100
N <- 10000
m_out <- matrix(NA, nIter, 2)

set.seed(1)

for(i in 1:nIter) {

  # GMM1: one dim
  data <- DataGen_GMM2(mu = GMM1[[1]],
                       Sigma = GMM1[[2]],
                       cats = "cont",
                       pi = c(.5, .5),
                       N = N)

  out1 <- Est_GGM(data = data$data, Kseq = 1:7)
  m_out[i, 1] <- out1$K

  # GMM2: distributed
  data <- DataGen_GMM2(mu = GMM2[[1]],
                       Sigma = GMM2[[2]],
                       cats = "cont",
                       pi = c(.5, .5),
                       N = N)

  out2 <- Est_GGM(data = data$data, Kseq = 1:7)
  m_out[i, 2] <- out2$K

  print(i)

} # end for: iter

apply(m_out==2, 2, mean)










