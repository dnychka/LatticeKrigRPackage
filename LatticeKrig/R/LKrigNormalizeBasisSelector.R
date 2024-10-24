# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2024
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

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
