source('CenterOutwardRegression.R')


set.seed(532342)

n.R<-40
n.S<-40
nn<-n.R*n.S
n_sample<-50000

rmixture<-function(a){
  aux<-runif(1)
  o1<-c(0,0)
  o2<-10*c(cos(-pi/6),sin(-pi/6))
  o3<-10*c(cos(7*pi/6),sin(7*pi/6))
  o4<-c(0,10)
  result<-(rnorm(2)+(aux<=0.25)*o1+(aux>0.25)*(aux<=0.5)*o2+(aux>0.5)*(aux<=0.75)*o3+(aux>0.75)*o4)/10
  return(result)
}

rotation<-function(x){
  angle<-pi*x/2
  rotmatrix<-matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),ncol=2)
  return(rotmatrix)
}

x <- 4 * runif(n_sample) - 2
ss<-1+(3/2)*sin(pi*x)^2
ee<-t(sapply(1:n_sample,rmixture))
for (i in 1:n_sample) ee[i,]<-ee[i,]%*%rotation(x[i])
y<-x+ss*ee[,1]
z<-x^2+ss*ee[,2]


xmin <- min(x)
xmax <- max(x)
grid <- 8
x_tick <- seq(xmin+0.5*(xmax-xmin)/grid,xmax-0.5*(xmax-xmin)/grid,by=(xmax-xmin)/grid)



fit<-QuantileRegressionOT(x,y,z,n.R,n.S,x_tick,ContLength=180)

rotMat<-matrix(rbind(c(0.74,0.67,0.02,0),
                     c(-0.08,0.108,0.99,0),
                     c(0.67,-0.731,0.131,0),
                     c(0,0,0,1)))

plotQuantileRegressionOT3D(fit,contours=TRUE,rays=FALSE)
par3d(windowRect=c(100,100,700,700))
view3d(userMatrix = rotMat)
aspect3d(1,1,1)



movie3d(spin3d(rpm=3), duration=20, clean=FALSE)


