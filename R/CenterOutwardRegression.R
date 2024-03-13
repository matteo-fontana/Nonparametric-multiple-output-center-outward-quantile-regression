## QuantileRegressionOT(x, y, z, nR, nS, x_tick,ContLength)
## main function for computing center-outward regression
## x,y,z column vectors; (y,z) output variables, x regressor
## nR, nS control the number of nearest neighbors (k = nR x nS)
## x_tick is the vector of x's at which conditional quantiles are computed
## ContLength is the number of points at which the 0.2, 0.4 and 0.8 contours are evaluated

QuantileRegressionOT<-function(x, y, z, nR, nS, x_tick,ContLength){
  library(transport)
  library(doParallel)
  grid<-length(x_tick)
  nc<-min(grid,detectCores())
  n.R<-nR
  n.S<-nS
  registerDoParallel(nc)
  n<-n.R*n.S
  d<-2
  ## generate uniform grid
  disc.sphere<-cbind(cos((2*pi)*(0:(n.S-1))/n.S ),sin((2*pi)*(0:(n.S-1))/n.S ))
  U.grid<-NULL
  for (i in 1:n.R) U.grid<-rbind(U.grid,i/(n.R+1)*disc.sphere)
  ## x mesh
  quantiles<-foreach (u=1:grid) %dopar%{
    ### select n-nearest neighbors
    ord = sort(abs(x - x_tick[u]),index.return=TRUE)$ix[1:n]
    X = x[ord]
    Y = y[ord]
    Z = z[ord]
    random.sample<-cbind(Y,Z)
    normas.y<-rep(0,n)
    for (i in 1:n) normas.y[i]<-sqrt(sum(random.sample[i,]^2))
    ysup<-max(normas.y)
    ##############################################
    ### Step 1: find optimal assignment
    ### We use transport (auctionbf algorithm)
    ### from 'transport' R package
    ##############################################
    asignacion<-transport(pp(U.grid),pp(random.sample),p=2,method='auctionbf')
    asignacion<-asignacion[,2]
    xxx<-U.grid/ysup
    yyy<-random.sample[asignacion,]/ysup
    ##############################################
    ### Step 2b: solve linear program (16)
    ###  via Karp's algorithm
    ##############################################
    cij<-list(apply(xxx*yyy,1,sum))
    cij<-do.call(cbind,rep(cij,n))
    cij<-cij-xxx%*%t(yyy)
    ind.diag<-cbind(1:n,1:n)
    cij[ind.diag]<-rep(Inf,n)
    dkv<-matrix(Inf,nrow=n,ncol=n+1)
    dkv[1,1]<-0
    start.time <- Sys.time()
    for (k in 2:(n+1)){
      aux.mat<-list(dkv[,(k-1)])
      aux.mat<-do.call(cbind,rep(aux.mat,n))+cij
      dkv[,k]<-apply(aux.mat,2,min)
    }
    dndk<-list(dkv[,(n+1)])
    dndk<-do.call(cbind,rep(dndk,n))-dkv[,1:n]
    denom<-list(n:1)
    denom<-do.call(rbind,rep(denom,n))
    dndk<-dndk/denom
    d.max<-apply(dndk,1,max)
    mu.star<-min(d.max)
    dkv.mu<-dkv-mu.star*matrix(1,ncol=1,nrow=n)%*%matrix(0:n,ncol=n+1,nrow=1)
    shortest.distances<-apply(dkv.mu,1,min)
    psi<-(-shortest.distances+shortest.distances[1])*ysup^2
    e0<-abs(mu.star)*ysup^2
    ##############################################
    ### Step 3: computation of cyclically monotone
    ### interpolation
    ##############################################
    xxx<-U.grid
    yyy<-random.sample[asignacion,]
    T.0<-function(z){
      auxfun<-function(i){return(sum(z*yyy[i,])-psi[i])}
      scores<-sapply(1:n,auxfun)
      indice<-which.max(scores)
      return(c(x_tick[u],yyy[indice,]))
    }
    centr<-T.0(c(0,0))
    disc.sphere<-cbind(cos(2*pi*(0:ContLength)/ContLength),sin(2*pi*(0:ContLength)/ContLength))
    quant1<-t(apply(0.2*disc.sphere,1,T.0))
    quant2<-t(apply(0.4*disc.sphere,1,T.0))
    quant3<-t(apply(0.8*disc.sphere,1,T.0))
    theta<-0
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray1<-t(apply(ru,1,T.0))
    theta<-pi/4
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray2<-t(apply(ru,1,T.0))
    theta<-pi/2
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray3<-t(apply(ru,1,T.0))
    theta<-3*pi/4
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray4<-t(apply(ru,1,T.0))
    theta<-pi
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray5<-t(apply(ru,1,T.0))
    theta<-5*pi/4
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray6<-t(apply(ru,1,T.0))
    theta<-3*pi/2
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray7<-t(apply(ru,1,T.0))
    theta<-7*pi/4
    ru<-cbind(cos(theta)*(0:8)/10,sin(theta)*(0:8)/10)
    ray8<-t(apply(ru,1,T.0))
    return(list(center=centr,quantile1=quant1,quantile2=quant2,quantile3=quant3,
                ray1=ray1,ray2=ray2,ray3=ray3,ray4=ray4,ray5=ray5,ray6=ray6,ray7=ray7,ray8=ray8))
  }
  center<-NULL
  for (a in 1:grid) center<-rbind(center,quantiles[[a]]$center)
  return(list(center=center,quantiles=quantiles,x=x,y=y,z=z,grid=grid,ContLength=ContLength,x_tick=x_tick))
}




plotQuantileRegressionOT3D<-function(QRObject,contours=TRUE,rays=FALSE){
  library(rgl)
  center<-QRObject$center
  lines3d(x=center[,1],y=center[,2],z=center[,3],col='red',lwd=2)
  aspect3d(1,1,1)
  highlevel()
  if(contours==TRUE){
  for(i in 1:QRObject$grid){
    lines3d(x=QRObject$quantiles[[i]]$quantile3[,1],
            y=QRObject$quantiles[[i]]$quantile3[,2],
            z=QRObject$quantiles[[i]]$quantile3[,3],
            col='darkolivegreen1',add=TRUE,lwd=2)
    lines3d(x=QRObject$quantiles[[i]]$quantile2[,1],
            y=QRObject$quantiles[[i]]$quantile2[,2],
            z=QRObject$quantiles[[i]]$quantile2[,3],
            col='darkolivegreen',add=TRUE,lwd=2)
    lines3d(x=QRObject$quantiles[[i]]$quantile1[,1],
            y=QRObject$quantiles[[i]]$quantile1[,2],
            z=QRObject$quantiles[[i]]$quantile1[,3],
            col='darkgreen',add=TRUE,lwd=2)
  }
  }
  if(rays==TRUE){
    for(i in 1:QRObject$grid){  
     lines3d(QRObject$quantiles[[i]]$ray1[,1],
             QRObject$quantiles[[i]]$ray1[,2],
             QRObject$quantiles[[i]]$ray1[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray2[,1],
             QRObject$quantiles[[i]]$ray2[,2],
             QRObject$quantiles[[i]]$ray2[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray3[,1],
             QRObject$quantiles[[i]]$ray3[,2],
             QRObject$quantiles[[i]]$ray3[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray4[,1],
             QRObject$quantiles[[i]]$ray4[,2],
             QRObject$quantiles[[i]]$ray4[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray5[,1],
             QRObject$quantiles[[i]]$ray5[,2],
             QRObject$quantiles[[i]]$ray5[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray6[,1],
             QRObject$quantiles[[i]]$ray6[,2],
             QRObject$quantiles[[i]]$ray6[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray7[,1],
             QRObject$quantiles[[i]]$ray7[,2],
             QRObject$quantiles[[i]]$ray7[,3],
             col='magenta4',add=TRUE,lwd=2)
     lines3d(QRObject$quantiles[[i]]$ray8[,1],
             QRObject$quantiles[[i]]$ray8[,2],
             QRObject$quantiles[[i]]$ray8[,3],
             col='magenta4',add=TRUE,lwd=2)
    }
  }
  grid3d(c('x','y+','z'), at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)
  axes3d(edges = c('x','y+','z'), labels = TRUE, tick = TRUE, nticks = 5, 
         box = FALSE,expand=1.05)
  title3d(main = NULL, sub = NULL, xlab = 'X', ylab = expression(Y[1]), 
          zlab = expression(Y[2]), line = NA) 
}
