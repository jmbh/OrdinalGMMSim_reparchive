# jonashaslbeck@gmail.com; Feb 09, 2022

# --------------------------------------------------------------
# ---------- Get Iteration Number ------------------------------
# --------------------------------------------------------------

# !/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
print(iter)
iter <- as.numeric(iter)

# ----------------------------------------------------------------------
# ----- Load Packages --------------------------------------------------
# ----------------------------------------------------------------------

library(MASS) # Data generation
library(mclust) # Estimation

# Functions for data generation and estimation
source("aux_functions.R")

# Parallel
library(foreach)
library(parallel)
library(doParallel)


# ----------------------------------------------------------------------
# ----- Load Models ----------------------------------------------------
# ----------------------------------------------------------------------

l_GMMs <- readRDS(file="l_GMMs.RDS")


# ----------------------------------------------------------------------
# ----- Iterate --------------------------------------------------------
# ----------------------------------------------------------------------

# ----- Some Specification for Simulation -----

n_seq  <- c(1000, 2500, 10000) # fixed Sample size
n_n_seq <- length(n_seq)
cats_seq <- c(2:12, "cont")


# ----- Iterating -----

# Setup parallelization
cluster <- 12
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)

timer_total <- proc.time()[3]

set.seed(iter)

out_iter <- foreach(cat = 1:12,
                    .packages = c("mclust", "mvtnorm", "MASS"),
                    .export = c("l_GMMs", "DataGen_GMM2", "Est_GGM"),
                    .verbose = TRUE) %dopar% {


                      # ----- Create Storage ------

                      # Structure for full models
                      l_mclust_out_i <- l_GMMs # Just to copy list structure
                      for(k in 1:3) for(d in 1:3) for(p in 1:9) l_mclust_out_i[[k]][[d]][[p]] <- vector("list", length=n_n_seq)

                      # Structure for K-estimates
                      a_Kest <- array(NA, dim = c(3,3,9,3)) # dim: K, dist, p, N

                      # ----- Loop Through Remaining factors ------

                      for(k in 1:3) {
                        for(d in 1:3) {
                          for(p in 1:9) {

                            # Load GMM
                            GMM_fix <- l_GMMs[[k]][[d]][[p]]

                            for(n in 1:3) {

                              # Generate Data
                              data <- DataGen_GMM2(mu = GMM_fix[[1]],
                                                   Sigma = GMM_fix[[2]],
                                                   cats = cats_seq[cat],
                                                   pi = rep(1/(k+1), k+1),
                                                   N = n_seq[n])

                              # Recover GMM
                              # !!! Adjust here to model="EEE" for constrained estimation (VVV=unconstrained estimation)
                              out <- Est_GGM(data = data$data, Kseq = 1:7, model="VVV")

                              # out$K

                              # Save Results
                              # K
                              a_Kest[k, d, p, n] <- out$K

                              # Parameters / full model
                              l_mclust_out_i[[k]][[d]][[p]][[n]] <- out$object$parameters # only save parameters to save space

                            } # end for: n
                          }

                        }
                      } # end for: k


                      # ----- Define Outlist & Return ------

                      outlist <- list("a_Kest" = a_Kest,
                                      "l_mclust_out_i" = l_mclust_out_i)

                      return(outlist)

                    } # end foreach


# print total time of nodes
print(paste0("Full Timing Iteration ", iter, ":"))
proc.time()[3] - timer_total

stopCluster(cl)


# ----------------------------------------------------------------------
# ----- Export ---------------------------------------------------------
# ----------------------------------------------------------------------


# ----- Combine cats-variations -----

# Array with K-estimates
a_Kest_w_cats <- array(NA, dim=c(3,3,9,12,3))
for(cat in 1:12) a_Kest_w_cats[, , , cat, ] <- out_iter[[cat]]$a_Kest

# List with everything
l_mclust_out <- list()
for(cat in 1:12) l_mclust_out[[cat]] <- out_iter[[cat]]$l_mclust_out_i

# Combine into one object
outlist_i <- list("K" =  a_Kest_w_cats,
                  "models" = l_mclust_out)

# Save
saveRDS(outlist_i, file = paste0("Simres_Iter", iter, ".RDS"))




