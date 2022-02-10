# jonashaslbeck@gmail.com; Feb 09, 2022

# ----------------------------------------------------------------------
# ----- Load Packages --------------------------------------------------
# ----------------------------------------------------------------------

library(RColorBrewer)
library(scales)
library(combinat)

source("aux_functions.R")

# ----------------------------------------------------------------------
# ----- Select Simulation: EEE (constrained) or VVV (unconstrained) ----
# ----------------------------------------------------------------------

# Unconstrained Estimation
# simResDir <- "output_VVV/"
# version <- "VVV"

# # Constrained Estimation
simResDir <- "output_EEE/"
version <- "EEE"


# ----------------------------------------------------------------------
# ----- Load Results ---------------------------------------------------
# ----------------------------------------------------------------------


v_files <- list.files(simResDir)
n_files <- length(v_files)

l_files <- list()
for(i in 1:n_files) l_files[[i]] <- readRDS(paste0(simResDir, v_files[i]))

# ----- Split up into K-array and model-list -----

# Model list
l_results <- list()
for(i in 1:n_files) l_results[[i]] <- l_files[[i]]$models

# K-estimate array
nIter <- n_files
a_Kest <- array(NA, dim=c(nIter, 3, 3, 9, 12, 3)) # now 3 n-variations
for(i in 1:nIter) a_Kest[i, , , , , ] <- l_files[[i]]$K

# saveRDS(a_Kest, file="files/a_Kest_VVV.RDS")
# saveRDS(l_results, file="files/l_results.RDS")

# ----------------------------------------------------------------------
# ----- Pre-process K-recovery Results ---------------------------------
# ----------------------------------------------------------------------

# Get Indicator correct classification
a_Kest_ind <- a_Kest
for(k in 1:3) a_Kest_ind[, k, , , , ] <- a_Kest[, k, , , , ] == k+1

# Average over iterations
a_Kest_indP <- apply(a_Kest_ind, 2:6, mean)
dim(a_Kest_indP)

mean(a_Kest_indP) # mean accuracy over everything



# ----------------------------------------------------------------------
# ----- Pre-process Estimation Error Results ---------------------------
# ----------------------------------------------------------------------

# --------- Load true graphs ----------
l_GMMs <- readRDS(file="files/l_GMMs.RDS")


# --------- Compute Errors for Means, Variances, Covariances ----------

# Storage: error array
a_errors <- array(NA, dim=c(nIter, 12, 3, 3, 9, 3, 3)) # iter categories, K, dist, p, n, type error (mean, var, cov)

for(i in 1:nIter) {
  for(cat in 1:12) {
    for(k in 1:3) {
      for(d in 1:3) {
        for(n in 1:3) {
          for(p in 1:9) {

            # True pars
            GMM_fix <- l_GMMs[[k]][[d]][[p]]
            mean_true <- GMM_fix[[1]]
            cov_true <- GMM_fix[[2]]

            # Est pars
            mean_est <- l_files[[i]]$models[[cat]][[k]][[d]][[p]][[n]]$mean
            mean_est
            cov_est <- l_files[[i]]$models[[cat]][[k]][[d]][[p]][[n]]$variance$sigma

            # Compute Error
            # (only defined if K estimated correctly)
            if(a_Kest[i, k, d, p, cat, n] == k+1) {

              ## Need to see which components correspond
              # I do this by computing all possible correspondences, and taking the one with the smallest error in means
              # There are K! possibilities

              k_perm <- permn(k+1)
              n_perm <- length(k_perm)
              m_error_perm <- matrix(NA, nrow=k+1, ncol=n_perm)

              for(k_i in 1:(k+1)){
                for(m in 1:n_perm) {
                  m_error_perm[k_i, m] <- mean(abs(mean_true[[k_i]] - mean_est[ ,k_perm[[m]][k_i]])) # absolute error
                }
              }
              mean_perm_err <- colMeans(m_error_perm)
              perm_match <- which.min(mean_perm_err) # this is the best fitting permutation

              ## OK, now we can compute the errors
              v_k_means <- v_k_vars <- v_k_covs <- rep(NA, k+1)
              for(k_i in 1:(k+1)){

                # Means
                v_k_means[k_i] <- mean(abs(mean_true[[k_i]] - mean_est[ ,k_perm[[perm_match]][k_i]])) # mean absolute error
                # Variances
                v_k_vars[k_i] <- mean(abs(diag(cov_true[[k_i]]) - diag(cov_est[ , , k_perm[[perm_match]][k_i]]))) # mean absolute error
                # Covariances
                v_k_covs[k_i] <- mean(abs(cov_true[[k_i]][upper.tri(cov_true[[k_i]])] - cov_est[ , ,k_perm[[perm_match]][k_i]][upper.tri(cov_est[ , ,k_perm[[perm_match]][k_i]])])) # mean absolute error

              } # end for: loop Ks

              if(is.na(mean(v_k_means))) stop("This shouldn't happen ....")

              ## Average over Ks and save
              a_errors[i, cat, k, d, p, n, 1] <- mean(v_k_means)
              a_errors[i, cat, k, d, p, n, 2] <- mean(v_k_vars)
              a_errors[i, cat, k, d, p, n, 3] <- mean(v_k_covs)

            } # end if: K estimated correctly
          }
        }
      }
    }
  }
  print(i)
} # end for: iterations


# --------- Aggregate across Iterations ----------

a_errors_agg <- apply(a_errors, 2:7, function(x) mean(x, na.rm = TRUE))

# --------- Inspect ----------

## Sanity Check:
tb <- table(is.na(a_errors))
round(tb / sum(tb), 4)
# This should match with:
mean(a_Kest_indP) # mean accuracy over everything
# OK!

# And this should also be at least as big was the above numbers
tb2 <- table(is.na(a_errors_agg))
round(tb2 / sum(tb2), 4)




# ----------------------------------------------------------------------
# ----- Main Results Figure --------------------------------------------
# ----------------------------------------------------------------------

# Select gradient colors for p-variation
cols <- brewer.pal(9, "Oranges")
# display.brewer.pal(9, "Blues")
cats_seq <- c(2:12, "ct")
n_cats <- length(cats_seq)

for(n in 1:3) {

  sc <- .6
  pdf(paste0("figures/Figure_X_MainResults_",version,"_n=", n,".pdf"), width = 10*sc, height = 9*sc)

  # ----- Define Layout -----

  lmat <- rbind(c(0, 1:3),
                c(4, 7:9),
                c(5, 10:12),
                c(6, 13:15))

  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1),
               heights = c(.15, 1, 1, 1))

  # layout.show(lo)

  # ----- Plot Labels -----

  # Top row
  plotLabels(expression("D"["KL"]*"(A || B) = 2"))
  plotLabels(expression("D"["KL"]*"(A || B) = 3.5"))
  plotLabels(expression("D"["KL"]*"(A || B) = 5"))

  # First column
  plotLabels("K = 2", srt=90)
  plotLabels("K = 3", srt=90)
  plotLabels("K = 4", srt=90)


  par(mar=c(3.2,3,1,1))

  # ----- Plot Results -----

  for(k in 1:3) {

    for(d in 1:3) {

      # Make Plot area
      plot.new()
      plot.window(xlim=c(1, 12), ylim=c(0, 1))
      axis(1, labels=cats_seq, at=1:n_cats, las=2)
      axis(2, las=2)

      for(p in 1:9) {
        lines(1:12, a_Kest_indP[k, d, p, , n], col = cols[p])
        # lines(1:12, a_Kest_indP[k, d, p, , 1], col = p)
      } # end for: p

      # add x-axis
      if(k==3) title(xlab="Number of categories", line=2)


      # --- Legend ---

      if(k==1) if(d==3)  legend("bottomright", legend=paste0("p = ", 2:10),
                                col=cols, lty=rep(1, 10), bty="n", cex=.8)


    } # end for: dist
  } # end for: k

  dev.off()


} # end for: n-variations


# ----------------------------------------------------------------------
# ----- Boxplots -------------------------------------------------------
# ----------------------------------------------------------------------

v_var_labels <- c("1k", "2.5k", "10k")

# Loop over three n-variations
for(n in 1:3) {

  sc <- 0.75
  pdf(paste0("figures/Fig_paper_Boxplots_",version,"_", v_var_labels[n] ,".pdf"), width = 10*sc, height = 9*sc)

  # ----- Define Layout -----
  par(mar=c(3.5,3.5,2,1), mfrow=c(3,3))


  # ----- Plot Results -----

  k <- 1
  d <- 1

  for(p in 1:9) {
    boxplot(a_Kest[, k, d, p, , n], ylim=c(0,7), axes=FALSE)
    title(main=paste0(p+1, " Variables"), font.main=1, cex.main=1, line=0.75)
    if(p %in% c(1, 4, 7)) title(ylab=expression(hat(K)), line=1.75)
    if(p %in% 7:9) title(xlab="Number of categories", line=2.2)


    axis(1, labels = cats_seq, at=1:12, las=1, cex.axis=1, las=2)
    axis(2, las=2)
  }

  dev.off()

} # end for: n-variations


# ----------------------------------------------------------------------
# ----- Show exact estimates for p=2 (indead Boxplots) -----------------
# ----------------------------------------------------------------------

# ## Note: Abandoned because it looks terrible
#
# v_var_labels <- c("1k", "2.5k", "10k")
#
# # Loop over three n-variations
# for(n in 1:3) {
#
#   sc <- 0.75
#   pdf(paste0("figures/Fig_paper_ExEst_p=2_",version,"_", v_var_labels[n] ,".pdf"), width = 10*sc, height = 9*sc)
#
#   # ----- Define Layout -----
#   par(mar=c(3.5,3.5,2,1), mfrow=c(3,3))
#
#   # ----- Plot Results -----
#
#   k <- 1
#   d <- 1
#
#   for(p in 1:9) {
#
#     # Layout
#     plot.new()
#     plot.window(xlim=c(1,12), ylim=c(0, 8))
#     axis(1, labels = cats_seq, at=1:12, las=1, cex.axis=1, las=1)
#     axis(2, las=2)
#     title(main=paste0(p+1, " Variables"), font.main=1, cex.main=1, line=0.75)
#     if(p %in% c(1, 4, 7)) title(ylab=expression(hat(K)), line=1.75)
#     if(p %in% 7:9) title(xlab="Number of categories", line=2.2)
#
#     # Data
#     for(cat in 1:12) points(rep(cat, nIter), a_Kest[, k, d, p, cat, n],
#                             pch=20, cex=2, col=alpha(colour = "black", alpha=0.05))
#
#   }
#
#   dev.off()
#
# } # end for: n-variations


# ----------------------------------------------------------------------
# ----- Heatplot Overview Figure ---------------------------------------
# ----------------------------------------------------------------------

# Change Jan 14, 2022: also include continuous condition

v_var_labels <- c("1k", "2.5k", "10k")

# ------ Create Color Gradient ------

color.gradient <- function(x, colors=c("red","blue"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
x <- 1:101
grad <- color.gradient(x)
grad <- alpha(grad, .75)

# Loop over three n-variations
for(n in 1:3) {

  # ------ Average across KDL and K, for N=10000 ------
  m_acc_agg_fixN <- apply(a_Kest_indP[, , , , n], 3:4, mean)
  m_acc_agg_fixN <- round(m_acc_agg_fixN, 2)


  # ------ Plotting ------

  library(RColorBrewer)
  library(scales)

  sc <- 0.9
  pdf(paste0("figures/Fig_Summary_",version,"_",v_var_labels[n],".pdf"), width=7.5*sc, height=7*sc)

  # Make canvas
  plot.new()
  plot.window(xlim=c(0,12), ylim=c(0, 9))

  # Auxiliary plotting variables
  seq_mp_x <- seq(.5, 11.5, length=12)
  seq_mp_y <- seq(.5, 8.5, length=9)
  sfm <- 1/2 # distance: midpoint to side

  # Plot Axis
  axis(1, labels = c(2:12, "ct"), at=seq_mp_x)
  title(xlab="Number of Ordinal Categories")
  axis(2, labels = 2:10, at=seq_mp_y, las=2)
  title(ylab="Number of Variables")

  # Title
  title(main = expression("Accuracy averaged across "*"D"["KL"]*" and K"),
        line=1.5, cex.main=1.5)


  # Plot Data
  for(i in 1:12) {
    for(j in 1:9) {

      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_y[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_y[j]+sfm,
           # col=alpha("darkgreen", m_acc_agg_fixN[j, i]),
           col = grad[m_acc_agg_fixN[j, i] * 100 + 1]
      )

      text(seq_mp_x[i], seq_mp_y[j], m_acc_agg_fixN[j,i], cex=.8)

    }
  }

  dev.off()

} # end for: n-variations



# ----------------------------------------------------------------------
# ----- Heatplot Overview Figure [1k and 2.5k next to each other] ------
# ----------------------------------------------------------------------

sc <- 0.9
pdf(paste0("figures/Fig_Summary_CMB_1and2.5_", version,".pdf"), width=7.5*sc*2, height=7.5*sc)

# Layout
par(mfrow=c(1,2))

# Data
for(n in 1:2) {

  # Aggregate for fixed n
  m_acc_agg_fixN <- apply(a_Kest_indP[, , , , n], 3:4, mean)
  m_acc_agg_fixN <- round(m_acc_agg_fixN, 2)

  # Make canvas
  plot.new()
  plot.window(xlim=c(0,12), ylim=c(0, 9))

  # Auxiliary plotting variables
  seq_mp_x <- seq(.5, 11.5, length=12)
  seq_mp_y <- seq(.5, 8.5, length=9)
  sfm <- 1/2 # distance: midpoint to side

  # Plot Axis
  axis(1, labels = c(2:12, "ct"), at=seq_mp_x, cex.axis=1.2)
  axis(2, labels = 2:10, at=seq_mp_y, las=2, cex.axis=1.2)

  title(xlab="Number of Ordinal Categories", cex.lab=1.3)
  title(ylab="Number of Variables", cex.lab=1.3)

  # Title
  if(n==1) title(main = expression("Accuracy across "*"D"["KL"]*" and K (N = 1000)"),
                 line=1.5, cex.main=1.7)
  if(n==2) title(main = expression("Accuracy across "*"D"["KL"]*" and K (N = 2500)"),
                 line=1.5, cex.main=1.7)


  # Plot Data
  for(i in 1:12) {
    for(j in 1:9) {

      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_y[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_y[j]+sfm,
           # col=alpha("darkgreen", m_acc_agg_fixN[j, i]),
           col = grad[m_acc_agg_fixN[j, i] * 100 + 1]
      )

      text(seq_mp_x[i], seq_mp_y[j], m_acc_agg_fixN[j,i], cex=.8)

    }
  }

} # end for N

dev.off()




# ----------------------------------------------------------------------
# ----- Evaluate Estimation Error --------------------------------------
# ----------------------------------------------------------------------

# ------ Create Color Gradient ------

color.gradient <- function(x, colors=c("blue","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
x <- 1:101
grad_EE <- color.gradient(x)
grad_EE <- alpha(grad_EE, .75)



# --------- Plot Main EE Results Figure ----------

sc <- 1.2
pdf(paste0("figures/Fig_EE_abserr_", version, ".pdf"), width = 9*sc, height = 8*sc)

# ----- Define Layout -----

lmat <- rbind(c(0, 1:3),
              c(4, 7:9),
              c(5, 10:12),
              c(6, 13:15))

lo <- layout(mat = lmat,
             widths = c(.15, 1, 1, 1),
             heights = c(.15, 1, 1, 1))

# layout.show(lo)

# ----- Plot Labels -----

# Top row
ypos <- 0.45
cex <- 2
plotLabels("Means", y=ypos, cex=cex)
plotLabels("Variances", y=ypos, cex=cex)
plotLabels("Covariances", y=ypos, cex=cex)

# First column
plotLabels("N = 1000", srt=90, cex=cex)
plotLabels("N = 2500", srt=90, cex=cex)
plotLabels("N = 10000", srt=90, cex=cex)


# ----- Plot Data/Heatplots -----

# tb <- table(is.na(a_errors_agg))
# round(tb/sum(tb), 2)

range(a_errors_agg, na.rm=TRUE)


for(n in 1:3) {

  # ------ Average across KDL and K, for n-vars ------

  a_errors_agg_agg2_fixN <- apply(a_errors_agg[, , , , n, ], c(1,4,5), function(x) mean(x, na.rm = TRUE))
  a_errors_agg_agg2_fixN <- round(a_errors_agg_agg2_fixN, 2)


  for(etype in 1:3) {

    # Make canvas
    par(mar=c(4,4,0,0))
    plot.new()
    plot.window(xlim=c(0,12), ylim=c(0, 9))

    # Auxiliary plotting variables
    seq_mp_x <- seq(.5, 11.5, length=12)
    seq_mp_y <- seq(.5, 8.5, length=9)
    sfm <- 1/2 # distance: midpoint to side

    # Plot Axis
    axis(1, labels = c(2:12, "ct"), at=seq_mp_x, cex.axis=1)
    if(n==3) title(xlab="Number of Ordinal Categories", cex.lab=1.5, line=2.7)
    axis(2, labels = 2:10, at=seq_mp_y, las=2, cex=1)
    if(etype==1) title(ylab="Number of Variables", line=2.2, cex.lab=1.5)

    # Plot Data
    for(i in 1:12) {
      for(j in 1:9) {

        if(!is.na(a_errors_agg_agg2_fixN[i, j, etype])) {
          rect(xleft = seq_mp_x[i]-sfm,
               ybottom = seq_mp_y[j]-sfm,
               xright = seq_mp_x[i]+sfm,
               ytop = seq_mp_y[j]+sfm,
               # col=alpha("darkgreen", m_acc_agg_fixN[j, i]),
               col = grad_EE[a_errors_agg_agg2_fixN[i, j, etype] * 2.5 * 100 + 1] # times four because max error = approx 0.25
          )
        }

        if(!is.na(a_errors_agg_agg2_fixN[i, j, etype])) {
          text(seq_mp_x[i], seq_mp_y[j], a_errors_agg_agg2_fixN[i, j, etype], cex=.6, col="black")
        }

      }
    }

  } # end for: error-type

} # end for: n-var

dev.off()





