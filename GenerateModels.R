# jonashaslbeck@gmail.com; Nov 25, 2021

# ----------------------------------------------------------------------
# ----- What are we doing here? ----------------------------------------
# ----------------------------------------------------------------------

# Here we generate the GMMs varied acros ...
# - K
# - pairwise KLD
# - number of variables
# ... which I use in the simulation study.


# ----------------------------------------------------------------------
# ----- Load Packages --------------------------------------------------
# ----------------------------------------------------------------------

library(MASS)
library(mclust)
library(scales)
library(RColorBrewer)
library(mvtnorm) # for densities for contours

source("aux_functions.R")


# ----------------------------------------------------------------------
# ----- Define Gaussian Mixtures ---------------------------------------
# ----------------------------------------------------------------------

# ----- Make List Structure -----

# 3 times nested: K, p, dist (and then later added: mu and Sigma)
l_p <- vector("list", length = 10)
l_Ks <- vector("list", length = 3)
l_GMMs <- vector("list", length = 3)
l_Ks <- lapply(l_Ks, function(x) l_p)
l_GMMs <- lapply(l_GMMs, function(x) l_Ks)


# ----- For p>2 we use numerical optimization -----

p_seq <- 2:10
k_seq <- 2:4
d_seq <- c(2, 3.5, 5)

# Use different intial values for different KL-divergences
m_init_ranges <- rbind(c(0, 1),
                       c(0, 1),
                       c(0, 1))

set.seed(1)

for(pi in 2:9) { # p=2 is specified below
  for(k in 1:3) {
    for(d in 1:3) {
      l_GMMs[[k]][[d]][[pi]] <- makeGGM_KLD(K = k_seq[k],
                                            p = p_seq[pi],
                                            target_KLD = d_seq[d],
                                            maxIter = 50,
                                            init_range = m_init_ranges[d, ],
                                            method = "Nelder-Mead",
                                            tol = 0.01,
                                            verbose = FALSE) # This does the job
    }
  }
  print(pi)
}

# ----- For p=2 we choose geometrically nice looking arrangements -----

sig <- 0.5

# K = 2
a_seq <- c(0.5, .662, 0.791)
Sigma <- list(diag(2)*sig, diag(2)*sig)

for(d in 1:3) {
  a <- a_seq[d]
  mu <- list(c(-a, a),
             c(a, -a))
  l_GMMs[[1]][[d]][[1]] <- list(mu, Sigma)
}

# K = 3
a_seq <- c(0.5, .662, 0.791)
Sigma <- list(diag(2)*sig, diag(2)*sig,  diag(2)*sig)

for(d in 1:3) {
  a <- a_seq[d]
  mu <- list(c(-a, a),
             c(a, -a),
             c(a*sqrt(3), a*sqrt(3)))
  l_GMMs[[2]][[d]][[1]] <- list(mu, Sigma)
}


# K = 4
v_sig_offd <- c(.195, .2133, .223)
a_add <- c(.09, .1289, .16)
a_seq_K4 <- a_seq + a_add

sig_diag <- 0.5

for(d in 1:3) {
  a <- a_seq_K4[d]
  mu <- list(c(-a, a),
             c(a, -a),
             c(-a, -a),
             c(a, a))

  sig_offd <- v_sig_offd[d]
  # Define Sigma for components on diagonal
  Sigma1 <- matrix(sig_diag, 2, 2)
  Sigma1[1, 2] <- Sigma1[2, 1] <- -sig_offd
  Sigma2 <- Sigma1
  # Define Sigma for components on off-diagonal
  Sigma3 <- matrix(sig_diag, 2, 2)
  Sigma3[1, 2] <- Sigma3[2, 1] <- sig_offd
  Sigma4 <- Sigma3
  # Combine in list
  Sigma <- list(Sigma1, Sigma2, Sigma3, Sigma3)


  l_GMMs[[3]][[d]][[1]] <- list(mu, Sigma)
}



# ----- Save GMMs for Simulation -----


saveRDS(l_GMMs, file="files/l_GMMs.RDS")











