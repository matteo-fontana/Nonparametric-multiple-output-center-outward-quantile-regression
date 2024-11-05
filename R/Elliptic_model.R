source('R/CenterOutwardRegression.R')

set.seed(1991)

n.R<-35
n.S<-35
nn<-n.R*n.S
n_sample <- 10000
x <- 2 * runif(n_sample) - 1
ss<-sqrt(1+(3/2)*sin(pi*x/2)^2)
ee<-matrix(rnorm(2*n_sample),ncol=2)
y<-sin((2/3)*pi*x)+0.575*ss*ee[,1]
z<-cos((2/3)*pi*x)+x^2+2.645*x^4+0.25*ee[,1]^3+ss*ee[,2]/2.3

xmin <- min(x)
xmax <- max(x)
grid <- 8
x_tick <- seq(xmin+0.5*(xmax-xmin)/grid,xmax-0.5*(xmax-xmin)/grid,by=(xmax-xmin)/grid)


fit<-QuantileRegressionOT(x,y,z,n.R,n.S, x_tick, ContLength = 50)

rotMat<-matrix(rbind(c(0.74,0.67,0.02,0),
                     c(-0.08,0.108,0.99,0),
                     c(0.67,-0.731,0.131,0),
                     c(0,0,0,1)))

plotQuantileRegressionOT3D(fit,contours=TRUE,rays=FALSE)
par3d(windowRect=c(100,100,700,700))
view3d(userMatrix = rotMat)
aspect3d(1,1,1)



movie3d( spin3d(rpm=3), duration=20,dir='/videos/', clean=FALSE )
