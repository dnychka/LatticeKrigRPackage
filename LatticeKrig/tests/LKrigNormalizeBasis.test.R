
suppressMessages(library( LatticeKrig))
#options( echo=FALSE)

##########################################
test.for.zero.flag<- 1

n<-500
set.seed(121)
x<-  cbind(runif(n, -1,1), runif(n,-1,1))
r<- ((x[,1]^2 + x[,2]^2)/2)
y<-  exp(-r*2) + .01*rnorm( n)


NC<- 5
a.wght<- 4.01
obj <- LatticeKrig( x, y, NC=NC, a.wght=a.wght, nlevel=1, 
                   alpha=1.0,
                   normalize=TRUE)

objFFT <- LatticeKrig( x, y, NC=NC,  a.wght=a.wght,
                       alpha=1.0,
                       normalize=TRUE,
                       normalizeMethod = "both",
                       nlevel=1)

gHat<- predictSurface( obj, nx=100, ny=100 )
gHatFFT <- predictSurface( objFFT, nx=100, ny=100 )

pred_diffs <- abs(c(gHat$z - gHatFFT$z))
errorStats<- stats( pred_diffs)


test.for.zero(errorStats["mean",] , 0, 
              relative = FALSE, tol = 8.5e-4)
test.for.zero(errorStats["median",] , 0, 
              relative = FALSE, tol = 6e-4)
test.for.zero(errorStats["max",] , 0, 
              relative = FALSE, tol = 7.8e-3)
test.for.zero(errorStats["min",] , 0, 
              relative = FALSE, tol = 3e-8)



