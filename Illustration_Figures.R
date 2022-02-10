# jonashaslbeck@gmail.com; Feb 6, 2022

# ----------------------------------------------------------------------
# ----- Load Packages & Aux functions ----------------------------------
# ----------------------------------------------------------------------

library(scales)
library(RColorBrewer)
library(MASS)
library(mclust)

source("aux_functions_GaussianMix.R")


# ----------------------------------------------------------------------
# ----- Illustration A: GMM Arrangements in p = 2 ----------------------
# ----------------------------------------------------------------------

# ----- Load Models -----

l_GMMs <- readRDS(file="files/l_GMMs.RDS")

sc <- .6
pdf("figures/Figure2_Kvar_SimSetupV3.pdf", width = 10*sc, height = 9*sc)

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
plotLabels(expression("D"["KL"]*" = 2"))
plotLabels(expression("D"["KL"]*" = 3.5"))
plotLabels(expression("D"["KL"]*" = 5"))

# First column
plotLabels("K = 2", srt=90)
plotLabels("K = 3", srt=90)
plotLabels("K = 4", srt=90)


par(mar=c(3,3,1,1))

# ----- K = 2 -----

for(i in 1:3) {

  # Get Models
  GMM_K2_i <- l_GMMs[[1]][[i]][[1]]

  # Plotting
  PlotCont(GMM_K2_i[[1]], GMM_K2_i[[2]])

} # end for: dist


# ----- K = 3 -----

for(i in 1:3) {

  # Get Models
  GMM_K3_i <- l_GMMs[[2]][[i]][[1]]

  # Plotting
  PlotCont(GMM_K3_i[[1]], GMM_K3_i[[2]])

} # end for: dist



# ----- K = 4 -----

for(i in 1:3) {

  # Get Models
  GMM_K4_i <- l_GMMs[[3]][[i]][[1]]

  # Plotting
  PlotCont(GMM_K4_i[[1]], GMM_K4_i[[2]])

} # end for: dist

dev.off()


# ----------------------------------------------------------------------
# ----- Illustration B: Mapping to Ordinal -----------------------------
# ----------------------------------------------------------------------

# Load Model
# Get Models
GMM_K2_i <- l_GMMs[[1]][[3]][[1]] # use d=3 (=5)

# Specifcy model
mu <- GMM_K2_i[[1]]
Sigma <- GMM_K2_i[[2]]

# Category variations
v_cats <- c("cont", 12, 10, 8, 6, 5, 4, 3, 2)
v_alpha <- c(.05, rep(.0075, 8)) # stronger alpha needed for continuous

# Plotting
sc <- 0.9
pdf("figures/Sim1_Illu_1NEW3.pdf", 8*sc, 8*sc)

# Layout
par(mfrow=c(3,3), mar=c(4,4,2,1))

for(cat in 1:9) {

  # Generate Data
  set.seed(1)
  dlist <- DataGen_GMM2(mu = mu,
                        Sigma = Sigma,
                        cats = v_cats[cat],
                        pi = rep(1/2,2),
                        N = 10000)

  # Plotting
  if(v_cats[cat]=="cont") {
    main_cat <- "Continuous"
  } else {
    main_cat <- paste0(v_cats[cat], " Categories")
  }

  # Subsample for continuous case so that picture becomes smaller
  if(cat==1) dlist$data <- dlist$data[sample(1:10000, size=1000), ]


  plot.new()
  plot.window(xlim=c(-4, 4), ylim=c(-4, 4))
  axis(1, seq(-4, 4, length=5))
  axis(2, seq(-4, 4, length=5), las=2)

  points(dlist$data[, 1], dlist$data[, 2],
         pch=20, cex=2, col=alpha("black", v_alpha[cat]),
         xlab="", ylab="")

  title(main=main_cat, font.main = 1)

}

dev.off()


# ----------------------------------------------------------------------
# ----- Illustration C: Mapping to Ordinal + Example Run ---------------
# ----------------------------------------------------------------------

# Load Model
GMM_K2_i <- l_GMMs[[1]][[3]][[1]]

# Specifcy model
mu <- GMM_K2_i[[1]]
Sigma <- GMM_K2_i[[2]]

# Plotting
sc <- 0.9
pdf("figures/Sim1_Illu_1NEW3_wEst.pdf", 8*sc, 8*sc)

# Layout
par(mfrow=c(3,3), mar=c(4,4,2,1))

for(cat in 1:9) {

  # Generate Data
  set.seed(1)
  dlist <- DataGen_GMM2(mu = mu,
                        Sigma = Sigma,
                        cats = v_cats[cat],
                        pi = rep(1/2,2),
                        N = 10000)

  # Plotting
  if(v_cats[cat]=="cont") {
    main_cat <- "Continuous"
  } else {
    main_cat <- paste0(v_cats[cat], " Categories")
  }

  # Subsample for continuous case so that picture becomes smaller
  if(cat==1)  dlist$data <- dlist$data[sample(1:10000, size=1000), ]


  plot.new()
  plot.window(xlim=c(-4, 4), ylim=c(-4, 4))
  axis(1, seq(-4, 4, length=5))
  axis(2, seq(-4, 4, length=5), las=2)

  points(dlist$data[, 1], dlist$data[, 2],
         pch=20, cex=2, col=alpha("black", v_alpha[cat]),
         xlab="", ylab="")

  title(main=main_cat, font.main = 1)

  # Plot estimated component means
  out <- Est_GGM(data = dlist$data, Kseq = 1:7)

  means <- out$object$parameters$mean

  points(means[1, ], means[2, ], col="red", pch=4, cex=2.5)

  if(cat==1) text(0, -3.8, "X = Estimated Component Means", col="red")


}

dev.off()



