
suppressMessages(library( LatticeKrig))
options( echo=FALSE)

#--------------------------------------------------------------
test.for.zero.flag<- 1
#a bit of param setup
xr<- c( -1,1)
yr<- c( -1,1)
numLevel <- 1
numBasis <- 5
numBuffer <- 10
l <- numLevel
a.wght <- 4.05

#my LKinfo object for testing the base functions
LKinfo<- LKrigSetup(cbind(xr, yr), nlevel=numLevel, NC= numBasis, NC.buffer = numBuffer,
                    a.wght= a.wght, normalize= FALSE)

#grid 
gridList<- list( x= seq( -1,1,length.out=40),
                 y= seq( -1,1,length.out=40) )
x1 <- make.surface.grid(gridList)


#testing the low level functions and LKrig.cov ---------------
#LKrig.cov surface 
look <- LKrig.cov(x1, LKinfo = LKinfo, marginal=TRUE )

#fft surface 
lookFFT <- LKrigNormalizeBasisFFTInterpolate(LKinfo, Level=l, x1)

#kronecker surface 
lookKroneck <- LKrigNormalizeBasisFast( LKinfo, Level=l, x1)

#cov and kroneck should match up exactly
test.for.zero(look, lookKroneck)

#fft and kroneck should be close
test.for.zero(mean(lookFFT - lookKroneck), 0.002, tol = 1e-3)
test.for.zero(max(lookFFT - lookKroneck), 0.07759, tol = 1e-3)
test.for.zero(min(lookFFT - lookKroneck), -0.0805, tol = 1e-3)
#--------------------------------------------------------------

#testing the similarity of fits and predictions ---------------
#making some data
sideLength <-40
numSampled <- 400


sGridList<- list( x= seq( -1,1,length.out=sideLength),
                  y=seq( -1,1,length.out=sideLength))
set.seed(777)
sGrid<- make.surface.grid(sGridList)
ind<- sample( 1: nrow(sGrid),numSampled, replace = FALSE )
x<- sGrid[ind,]
r<- ((x[,1]^2 + x[,2]^2)/2)
y<-  exp(-r*2) + .01*rnorm( numSampled)

#exact kronecker normalization fit
obj<- LatticeKrig( x, y, 
                    NC=numBasis, a.wght=a.wght, nlevel=l, 
                    alpha=1.0,
                    normalize=TRUE,
                    normalizeMethod="exactKronecker", 
                    NC.buffer = numBuffer)


#fft normalization fit
objFFT <- LatticeKrig( x, y, 
                  NC=numBasis,  a.wght=a.wght, nlevel=l,
                  alpha=1.0,
                  normalize=TRUE,
                  normalizeMethod = "fftInterpolation",
                  NC.buffer = numBuffer)

#fft normalization fit (but using "both" and the selector function)
objBOTH <- LatticeKrig( x, y, 
                       NC=numBasis,  a.wght=a.wght, nlevel=l,
                       alpha=1.0,
                       normalize=TRUE,
                       normalizeMethod = "both",
                       NC.buffer = numBuffer)

#no normalization fit 
objNONE <- LatticeKrig( x, y, 
                      NC=numBasis,  a.wght=a.wght, nlevel=l,
                      alpha=1.0,
                      normalize=FALSE,
                      NC.buffer = numBuffer)



#exact kroneck normalization predictions
gHat<- predict( obj,sGrid)

# fft normalization predictions
gHatFFT<- predict( objFFT,sGrid)

# fft normalization predictions (using BOTH)
gHatBOTH<- predict( objBOTH,sGrid)

# no normalization predictions
gHatNONE<- predict( objNONE,sGrid)

#testing for differences between
test.for.zero(mean(gHat), mean(gHatFFT), tol = 1e-4)
test.for.zero(mean(gHat), mean(gHatNONE), tol = 1e-4)
test.for.zero(mean(gHatFFT), mean(gHatNONE), tol = 1e-5)
#the both and fft should be identical
test.for.zero(gHatBOTH, gHatFFT)

#testing for magnitudes of FFT error
errorFFT <- abs(gHatFFT - gHat)
test.for.zero(mean(errorFFT), 4.333e-4, 1e-4)
test.for.zero(min(errorFFT), 3.735e-7, 1e-3)
test.for.zero(max(errorFFT), 2.517e-3, 1e-4)

#--------------------------------------------------------------

#making sure that the different grid bug no longer exists
varFFT <- LKrigNormalizeBasisFFTInterpolate(objFFT$LKinfo, Level=1, sGrid)
test.for.zero(varFFT, lookFFT)
