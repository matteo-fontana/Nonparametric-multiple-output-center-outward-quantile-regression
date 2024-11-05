# Quantile Regression with Optimal Transport
# This implementation provides tools for computing and visualizing center-outward regression
# using optimal transport methods for multivariate quantile regression.

#' Compute center-outward regression using optimal transport
#' 
#' @param x Numeric vector - the regressor/predictor variable
#' @param y,z Numeric vectors - the two output/response variables
#' @param nR Integer - first parameter controlling the number of nearest neighbors
#' @param nS Integer - second parameter controlling the number of nearest neighbors
#'        (total neighbors k = nR × nS)
#' @param x_tick Numeric vector - points at which conditional quantiles are computed
#' @param ContLength Integer - number of points for evaluating contours (0.2, 0.4, 0.8)
#' @return List containing center points, quantiles, and input data

QuantileRegressionOT <- function(x, y, z, n.R, n.S, x_tick, ContLength) {
  library(transport)
  library(doFuture)

  # Setup parallel processing
  grid <- length(x_tick)
  registerDoFuture()
  nc <- makeCluster(min(grid, detectCores()))
  plan(cluster,workers = nc)
  
  n <- n.R * n.S #n.R number of radii, n.S number of spheres
  d <- 2  # dimension of the response
  
  # Generate uniform grid on unit circle
  # Creates nS equally spaced points on the unit circle
  disc.sphere <- cbind(
    cos((2*pi)*(0:(n.S-1))/n.S),
    sin((2*pi)*(0:(n.S-1))/n.S)
  )
  
  # Create nested uniform grid
  # Scales the unit circle points by factors of i/(n.R+1) for i in 1:n.R
  U.grid <- NULL
  for (i in 1:n.R) {
    U.grid <- rbind(U.grid, i/(n.R+1)*disc.sphere)
  }
  
  
  
  # Compute quantiles for each x_tick point in parallel
  p = progressor(along = grid)
  quantiles <- foreach(u=1:grid,.options.future = list(seed = TRUE)) %dofuture% {
    # Select n nearest neighbors to current x_tick point
    ord <- sort(abs(x - x_tick[u]), index.return=TRUE)$ix[1:n]
    X <- x[ord]
    Y <- y[ord]
    Z <- z[ord]
    random.sample <- cbind(Y, Z)
    
    # Normalize the response data
    normas.y <- rep(0,n)
    for (i in 1:n) normas.y[i] <- sqrt(sum(random.sample[i,]^2))
    ysup <- max(normas.y)
    
    # Step 1: Optimal Transport Assignment
    # Uses the auction algorithm to find optimal mapping between
    # uniform grid and normalized response data
    assignment <- transport(pp(U.grid), pp(random.sample), p=2, method='auctionbf')
    assignment <- assignment$to #assuming unitary mass... maybe we can relax this?
    
    xxx <- U.grid/ysup
    yyy <- random.sample[assignment,]/ysup
    
    # Step 2: Solve Linear Program using Karp's Algorithm
    # Computes optimal cyclically monotone extension
    cij <- list(apply(xxx*yyy, 1, sum))
    cij <- do.call(cbind, rep(cij,n))
    cij <- cij - xxx%*%t(yyy)
    ind.diag <- cbind(1:n, 1:n)
    cij[ind.diag] <- rep(Inf,n)
    
    # Dynamic programming matrix
    dkv <- matrix(Inf, nrow=n, ncol=n+1)
    dkv[1,1] <- 0
    
    # Fill dynamic programming table
    for (k in 2:(n+1)) {
      aux.mat <- list(dkv[,(k-1)])
      aux.mat <- do.call(cbind,rep(aux.mat,n)) + cij
      dkv[,k] <- apply(aux.mat, 2, min)
    }
    
    # Compute optimal value and potential function
    dndk <- list(dkv[,(n+1)])
    dndk <- do.call(cbind,rep(dndk,n)) - dkv[,1:n]
    denom <- list(n:1)
    denom <- do.call(rbind,rep(denom,n))
    dndk <- dndk/denom
    d.max <- apply(dndk, 1, max)
    mu.star <- min(d.max)
    
    # Compute potential function
    dkv.mu <- dkv - mu.star*matrix(1,ncol=1,nrow=n)%*%matrix(0:n,ncol=n+1,nrow=1)
    shortest.distances <- apply(dkv.mu, 1, min)
    psi <- (-shortest.distances + shortest.distances[1])*ysup^2
    e0 <- abs(mu.star)*ysup^2
    
    # Step 3: Compute cyclically monotone interpolation
    xxx <- U.grid
    yyy <- random.sample[assignment,]
    
    # Define transport map
    T.0 <- function(z) {
      auxfun <- function(i) { return(sum(z*yyy[i,]) - psi[i]) }
      scores <- sapply(1:n, auxfun)
      indice <- which.max(scores)
      return(c(x_tick[u], yyy[indice,]))
    }
    
    # Compute center point and quantile contours
    centr <- T.0(c(0,0))
    disc.sphere <- cbind(
      cos(2*pi*(0:ContLength)/ContLength),
      sin(2*pi*(0:ContLength)/ContLength)
    )
    
    # Compute quantile contours at different levels
    quant1 <- t(apply(0.2*disc.sphere, 1, T.0))  # 20% contour
    quant2 <- t(apply(0.4*disc.sphere, 1, T.0))  # 40% contour
    quant3 <- t(apply(0.8*disc.sphere, 1, T.0))  # 80% contour
    
    # Compute directional rays at different angles
    # These help visualize the spread of the conditional distribution
    rays <- list()
    for(theta in seq(0, 7*pi/4, by=pi/4)) {
      ru <- cbind(cos(theta)*(0:8)/10, sin(theta)*(0:8)/10)
      rays[[length(rays) + 1]] <- t(apply(ru, 1, T.0))
    }
    
    # Return results for this x_tick point
    c(
      list(center=centr, quantile1=quant1, quantile2=quant2, quantile3=quant3),
      setNames(rays, paste0("ray", 1:8))
    )
    
  }
  
  # Collect center points across all x_tick values
  center <- NULL
  for (a in 1:grid) center <- rbind(center, quantiles[[a]]$center)
  
  # Return complete results
  return(list(
    center=center,
    quantiles=quantiles,
    x=x, y=y, z=z,
    grid=grid,
    ContLength=ContLength,
    x_tick=x_tick
  ))
}

#' Create 3D visualization of quantile regression results
#' 
#' @param QRObject Output from QuantileRegressionOT function
#' @param contours Logical - whether to plot quantile contours
#' @param rays Logical - whether to plot directional rays
#' @return None (creates 3D plot)
plotQuantileRegressionOT3D <- function(QRObject, contours=TRUE, rays=FALSE) {
  library(rgl)
  
  # Plot center curve in red
  center <- QRObject$center
  lines3d(
    x=center[,1], y=center[,2], z=center[,3],
    col='red', lwd=2
  )
  
  # Set aspect ratio and high-level plotting parameters
  aspect3d(1,1,1)
  highlevel()
  
  # Plot quantile contours if requested
  if(contours == TRUE) {
    for(i in 1:QRObject$grid) {
      # 80% contour in light green
      lines3d(
        x=QRObject$quantiles[[i]]$quantile3[,1],
        y=QRObject$quantiles[[i]]$quantile3[,2],
        z=QRObject$quantiles[[i]]$quantile3[,3],
        col='darkolivegreen1', add=TRUE, lwd=2
      )
      
      # 40% contour in medium green
      lines3d(
        x=QRObject$quantiles[[i]]$quantile2[,1],
        y=QRObject$quantiles[[i]]$quantile2[,2],
        z=QRObject$quantiles[[i]]$quantile2[,3],
        col='darkolivegreen', add=TRUE, lwd=2
      )
      
      # 20% contour in dark green
      lines3d(
        x=QRObject$quantiles[[i]]$quantile1[,1],
        y=QRObject$quantiles[[i]]$quantile1[,2],
        z=QRObject$quantiles[[i]]$quantile1[,3],
        col='darkgreen', add=TRUE, lwd=2
      )
    }
  }
  
  # Plot directional rays if requested
  if(rays == TRUE) {
    for(i in 1:QRObject$grid) {
      # Plot 8 rays at 45-degree intervals
      for(j in 1:8) {
        ray <- QRObject$quantiles[[i]][[paste0("ray",j)]]
        lines3d(
          ray[,1], ray[,2], ray[,3],
          col='magenta4', add=TRUE, lwd=2
        )
      }
    }
  }
  
  # Add grid, axes, and labels
  grid3d(c('x','y+','z'), at=NULL, col="gray", lwd=1, lty=1, n=5)
  axes3d(
    edges=c('x','y+','z'),
    labels=TRUE, tick=TRUE, nticks=5,
    box=FALSE, expand=1.05
  )
  title3d(
    main=NULL, sub=NULL,
    xlab='X',
    ylab=expression(Y[1]),
    zlab=expression(Y[2]),
    line=NA
  )
}