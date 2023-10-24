# jonashaslbeck@gmail.com; Feb 09, 2022

# ----------------------------------------------------------------------
# ----- Generate GMMs with given K, p and pw KLD -----------------------
# ----------------------------------------------------------------------

makeGGM_KLD <- function(K,
                        p,
                        target_KLD,
                        maxIter=20,
                        init_range = c(-1, 1),
                        sig = 0.25,
                        method = "Nelder-Mead",
                        tol = .001,
                        verbose = FALSE) {

  # --- Error Function ---

  fn <- function(par, target_KLD, K, p, sig) {

    # Create data structure for mixtures
    Sigma <- list()
    for(k in 1:K) Sigma[[k]] <- diag(p)*sig

    m_mu <- matrix(par, K, p, byrow=TRUE)
    mu <- list()
    for(k in 1:K) mu[[k]] <- m_mu[k, ]

    # Compute all pairwise KLD
    v_KLD <- rep(NA, K*(K-1)/2)
    count <- 1
    for(i in 1:K) {
      for(j in i:K) {
        if(j!=i) {
          v_KLD[count] <-  KLD(mu1 = mu[[i]],
                               mu2 = mu[[j]],
                               sig1 = Sigma[[i]],
                               sig2 = Sigma[[j]])
          count <- count + 1
        }
      }
    }

    # Compute Error
    v_error <- (v_KLD - target_KLD)^2

    return(sum(v_error))

  } # eof

  # --- Call Optim ---

  conv_err <- 20
  counter <- 1
  while(conv_err > tol) {

    n_pars <- K*p # number of mean-parameters
    par <- runif(n_pars, min = init_range[1], max = init_range[2])

    out <- optim(par = par,
                 fn = fn,
                 target_KLD = target_KLD,
                 K = K,
                 p = p,
                 sig = sig,
                 method = method)


    conv_err <- out$value
    if(verbose) print(out$value)
    counter <- counter + 1

    if(counter > maxIter) stop(paste0("No convergence after ", maxIter," tries."))

  }

  # --- Put Parameters in List Structure ---

  # Check results:
  m_mu_optim <- matrix(out$par, K, p, byrow=TRUE)
  mu <- list()
  for(k in 1:K) mu[[k]] <- m_mu_optim[k, ]
  Sigma <- list()
  for(k in 1:K) Sigma[[k]] <- diag(p)*sig

  outlist <- list("mu" = mu,
                  "Sigma" = Sigma,
                  "fin.error" = out$value,
                  "restarts" = counter)

  return(outlist)


} # eoF


# ----------------------------------------------------------------------
# ----- Plotting Labels in Designs -------------------------------------
# ----------------------------------------------------------------------

plotLabels <- function(tex, srt=0, x=0.5, y=0.5, cex=1.2) {

  par(mar=rep(0, 4))
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))
  text(x, y, adj=0.3, labels=tex, srt=srt, cex=cex)

}



# ----------------------------------------------------------------------
# ----- Make Contour plots of bivariate Gaussian Mixtures --------------
# ----------------------------------------------------------------------

PlotCont <- function(mu, Sigma, drawlabels=FALSE) {

  # Number of components
  K <- length(mu)

  # Select colors
  n <- ifelse(K<3, 3, K)
  cols <- brewer.pal(n = n, name = "Set1")[1:K]

  # Make Plot area
  plot.new()
  plot.window(xlim=c(-4, 4), ylim=c(-4, 4))
  axis(1, seq(-4, 4, length=5))
  axis(2, seq(-4, 4, length=5), las=2)

  # Plot contours
  x.points <- seq(-4,4,length.out=100)
  y.points <- x.points
  z <- matrix(0,nrow=100,ncol=100)
  for(k in 1:K) {
    mu_k <- mu[[k]]
    sigma_k <- Sigma[[k]]

    for (i in 1:100) {
      for (j in 1:100) {
        z[i,j] <- mvtnorm::dmvnorm(c(x.points[i],y.points[j]),
                                   mean=mu_k,
                                   sigma=sigma_k)
      }
    }

    contour(x.points,y.points,z, add=T, col=cols[k], drawlabels = drawlabels)

  } # end for: k

} # eoF




# ----------------------------------------------------------------------
# ----- Compute KL-Divergence for Multivariate Gaussian ----------------
# ----------------------------------------------------------------------
# ----- General Formula -----

KLD <- function(mu1, mu2, sig1, sig2) {

  d <- length(mu1)
  sig1 <- matrix(sig1, d, d)
  sig2 <- matrix(sig2, d, d)
  mu1 <- matrix(mu1, nrow=d)
  mu2 <- matrix(mu2, nrow=d)

  kld <- 0.5 * (    log( det(sig2)/det(sig1) )
                    - d
                    + sum(diag(solve(sig2) %*% sig1))
                    + t(mu1 - mu2) %*% solve(sig2) %*% (mu1 - mu2)
  )

  as.numeric(kld)
}


# ----------------------------------------------------------------------
# ----- Function to Generate Data from GMM  [v2 with quantil def] ------
# ----------------------------------------------------------------------

DataGen_GMM2 <- function(mu, # list of mean vectors
                         Sigma, # list of covariance matrices
                         cats = "cont", # number of Likert-style categories; if NULL-> continuous
                         pi, # mixing vector
                         N,
                         qntl_outside=0.01) {


  K <- length(mu) # Number of components

  # Sample from individual components

  l_data <- list()
  l_trueclasses <- list()

  for(k in 1:K) {

    N_k <- round(N * pi[k])

    l_data[[k]] <- mvrnorm(n = N_k,
                           mu = mu[[k]],
                           Sigma = Sigma[[k]])

    l_trueclasses[[k]] <- rep(k, N_k)

  } # loop over K

  # Collapse
  data_out <- data <- do.call(rbind, l_data)
  trueclasses <- do.call(base::c, l_trueclasses)

  # Map to Likert-scale
  if(cats != "cont") {
    cats <- as.numeric(cats)
    p <- length(mu[[1]])
    data_L <- data_Lrec <- matrix(NA, nrow(data), p)

    for(i in 1:p) {

      X_i <- data[, i]
      qtls <- quantile(X_i, probs = c(qntl_outside, 1-qntl_outside))
      threshs <- seq(qtls[1], qtls[2], length = cats+1)

      # Define category levels within range, so that axes don't change; specifically: mean between thresholds
      cat_labels <- (threshs[-(cats+1)] + threshs[-1]) /2

      # Apply thresholding
      data_L[, i] <- cut(data[, i], breaks=threshs, labels=FALSE, include.lowest = TRUE)

      # Collapse extreme values in lowest/highest category; this is necessary for creating reasonable thresholds for Gaussian's which have infinite support
      data_L[, i][data[, i] <= qtls[1]] <- 1
      data_L[, i][data[, i] >= qtls[length(qtls)]] <- cats

    }


    # Code to category levels
    for(j in 1:cats) data_Lrec[data_L==j] <- cat_labels[j]

    data_out <- data_Lrec
  } # End: ordinal scaling


  # Return
  outlist <- list("data" = data_out,
                  "trueclasses" = trueclasses)

  return(outlist)

} # eoF


# ----------------------------------------------------------------------
# ----- Function to Estimate GMM with mclust ---------------------------
# ----------------------------------------------------------------------

Est_GGM <- function(data, Kseq = 1:5, model="EEE") {

  # TODO later:
  # - Model selection with different methods

  # Call mclust
  mod <- Mclust(data, G=Kseq, modelNames = model, verbose=F)

  # Later: return parameter estimates
  # pars <- summary(mod2, parameters=TRUE)
  # pars

  outlist <- list("object" = mod,
                  "K" = mod$G)

  return(outlist)

} # eoF


