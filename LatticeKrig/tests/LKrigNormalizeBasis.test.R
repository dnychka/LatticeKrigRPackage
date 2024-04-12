
suppressMessages(library( LatticeKrig))
options( echo=FALSE)

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

pred_diffs <- na.omit(c(gHat$z - gHatFFT$z))
mean_error <- mean(pred_diffs)
max_error <- max(pred_diffs)
min_error <- min(pred_diffs)

test.for.zero(mean_error, 6.73e-5, tol = 1e-3)
test.for.zero(max_error, 0.00632, tol = 1e-3)
test.for.zero(min_error, -0.00583, tol = 1e-3)

