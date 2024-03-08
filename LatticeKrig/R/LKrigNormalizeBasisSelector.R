
# This functions choose which method to evaluates the variance of the 
# basis functions. It selects the fft method for lower levels where the
# number of basis functions is low enough, and then begins calling the
# Kronecker method for finer levels. This merges the two powerful 
# normalization methods in LatticeKrig

LKrigNormalizeBasisSelector <- function(LKinfo, Level, x1){
  
  # Extracting number of basis functions from LKinfo 
  basisNum <- LKinfo$latticeInfo$mxDomain[Level,2]
  
  #dimensions of original data
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  minDimension <- min(c(nr, nc))
  
  #coarse grid size for fft method
  miniGridSize <- 4 * basisNum
  
  #method selection
  # if coarse grid size is less than the size of the data, use FFT
  if (miniGridSize < minDimension){
    cat("Using FFT Interpolation method for level", Level, fill = TRUE)
    wght <- LKrigNormalizeBasisFFTInterpolate(LKinfo, Level, x1)
  }
  
  #when coarse grid size gets too big (too many basis functions)
  #switch over to the Kronecker method
  else {
    cat("Using Kronecker method for level", Level, fill = TRUE)
    wght <- LKrigNormalizeBasisFast(LKinfo,  Level,  x1)
  }
  
  return(wght)
  
}