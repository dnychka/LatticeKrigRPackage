
# This functions evaluates the variance of the basis functions on a coarser grid, 
# then uses 2D interpolation via FFT in order to smooth/interpolate the variance 
# up to the size of the original grid that we were working with. Should provide a 
# significant computational speedup. 


# NOTE: this function is reliant on the library(fftwtools). Needs to be installed 
#      with the LatticeKrig package itself. 




# Not sure if I named this correctly with the .LKRectangle
LKrigNormalizeBasisFFTInterpolate <- function(LKinfo, Level, x1){

  
  # Extracting important information from LKinfo 
  bounds <- LKinfo$x
  basisNum <- LKinfo$latticeInfo$mxDomain[Level,1]
  buffer <- LKinfo$NC.buffer
  alphaNum <- LKinfo$alpha[Level]
  awght <- LKinfo$a.wght[Level]
  
  #dimensions of original data, also helpful for shift parameter
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  
  # even number of rows or columns will result in asymmetric error due to fft assuming periodicity
  if (nr %% 2 == 0 || nc %% 2 == 0) {
    warning("Your input has an even number of rows or columns. 
          This will result in asymmetric error pattern, yet similar accuracy. 
          Considering using odd numbers for symmetric error.
            See LKrigNormalizeBasisFFTInterpolate help file.")
  }
  
  # Setting up a new LKrig object using the extracted info 
  LKinfoNew <- LKrigSetup(bounds, nlevel = 1, NC = basisNum, NC.buffer = buffer, 
                          alpha = alphaNum, a.wght = awght, normalize = FALSE) 
  
  # Setting a default coarse grid size based on the number of basis functions 
  # MINIMUM VALUE is 2 * basisNum + 2
  # NOTE can play with this for accuracy ----------------------------------------------------------- NOTE
  miniGridSize <- 4 * basisNum + 2
  
  # Creating the actual grid
  gridList<- list( x= seq( bounds[1,1],bounds[2,1],length.out = miniGridSize),
                   y= seq( bounds[1,2],bounds[2,2],length.out = miniGridSize) )
  
  sGrid<- make.surface.grid(gridList)
  
  
  # Calling LKrig.cov to evaluate the variance on our more coarse grid
  roughMat <- as.surface(sGrid, LKrig.cov(sGrid, LKinfo = LKinfoNew, marginal=TRUE ))[["z"]]
  
  # FFT step (reliant on fftwtools package)
  fftStep <- fftw2d((roughMat)/length(roughMat))
  
  # Helpful dimensions and shift parameter 
  snr <- nrow(fftStep)
  snc <- ncol(fftStep)
  
  yShift <- ceiling(((nr/snr) - 1)/2)
  xShift <- ceiling(((nc/snc) - 1)/2)
  
  # Instantiate empty matrix
  temp <- matrix(0, nrow = nr, ncol = nc)
  
  # Helpful for indexing later 
  bigindY <- 1:ceiling(snr/2)
  bigindX <- 1:ceiling(snc/2)
  indY <- 1:(snr/2)
  indX <- 1:(snc/2)
  bigOffsetY <- (nr - floor(snr/2))
  bigOffsetX <- (nc - floor(snc/2))
  smallOffsetY <- (snr - floor(snr/2))
  smallOffsetX <- (snc - floor(snc/2))
  
  # Stuffing the small FFT result into the large matrix of zeroes 
  temp[bigindY, bigindX] <- fftStep[bigindY, bigindX] #top left corner
  temp[indY, (indX + bigOffsetX)] <- fftStep[indY, (indX + smallOffsetX)] #top right corner
  temp[(indY + bigOffsetY), indX] <- fftStep[(indY + smallOffsetY), indX] #bottom left corner 
  temp[(indY + bigOffsetY), (indX + bigOffsetX)] <- fftStep[(indY + smallOffsetY), (indX + smallOffsetX)] #bottom right corner
  
  
  # Takes the IFFT of the modified big matrix to return a interpolated version of the input matrix
  # Again (reliant on fftwtools package)
  wght <- Re(fftw2d(temp, inverse = 1))
  
  # Shifting due to the periodicity assumed by the FFT
  wght <- LKrig.shift.matrix(wght, shift.row = yShift, shift.col = xShift, periodic = c(TRUE, TRUE))
  
  
  wght <- c(wght)
  
  # Returns matrix of weights (variances)
  return (wght)
  
}



