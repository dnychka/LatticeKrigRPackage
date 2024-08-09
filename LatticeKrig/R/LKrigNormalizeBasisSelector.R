LKrigNormalizeBasisSelector <- function(LKinfo, Level, x1, verbose){
  
  # Extracting number of basis functions from LKinfo 
  basisNum <- max(LKinfo$latticeInfo$mxDomain[Level,1], LKinfo$latticeInfo$mxDomain[Level,2])
  
  # Dimensions of original data
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  maxDimension <- max(nr, nc)
  
  # Cutoff for where the fft becomes essentially less useful than Kronecker
  miniGridSize <- 4 * basisNum
  
  #method selection
  # if coarse grid size is less than the size of the data, use FFT
  if (Level < 2 && miniGridSize < maxDimension){
    if (verbose){
      cat("Using FFT Interpolation method for level", Level, fill = TRUE)
    }
    wght <- LKrigNormalizeBasisFFTInterpolate(LKinfo, Level, x1)
  }
  
  #when coarse grid size gets too big (too many basis functions)
  #switch over to the Kronecker method
  else {
    if (verbose){
      cat("Using Kronecker method for level", Level, fill = TRUE)
    }
    wght <- LKrigNormalizeBasisFast(LKinfo,  Level,  x1)
  }
  
  return(wght)
}
