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
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2


LKrigNormalizeBasisFFTInterpolate <- function(LKinfo, Level, x1){
  # This functions evaluates the variance of the basis functions on a coarser grid, 
  # then uses 2D interpolation via FFT in order to smooth/interpolate the variance 
  # up to the size of the original grid that we were working with. Should provide a 
  # significant computational speedup. 

  # Extracting important information from LKinfo 
  bounds <- cbind(c(min(LKinfo$x[,1]), max(LKinfo$x[,1])), 
                  c(min(LKinfo$x[,2]), max(LKinfo$x[,2])))
  basisNum_big <- max(LKinfo$latticeInfo$mxDomain[Level,1], 
                      LKinfo$latticeInfo$mxDomain[Level,2])
  basisNum_small <- min(LKinfo$latticeInfo$mxDomain[Level,1], 
                        LKinfo$latticeInfo$mxDomain[Level,2])
  gridOrientation <- which.max(c(LKinfo$latticeInfo$mxDomain[Level,1],
                                 LKinfo$latticeInfo$mxDomain[Level,2]))

  buffer <- LKinfo$NC.buffer
  alphaNum <- LKinfo$alpha[Level]
  awght <- LKinfo$a.wght[Level]
  
  #dimensions of original data, also helpful for shift parameter
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  maxDimension <- max(nr, nc)
  minDimension <- min(nr, nc)
  
  # Setting up a new, single level LKrig object using the extracted info 
  LKinfoNew <- LKrigSetup(bounds, nlevel = 1, NC = basisNum_big, NC.buffer = buffer, 
                          alpha = alphaNum, a.wght = awght, normalize = FALSE) 
  
  # Setting a default coarse grid size based on the number of basis functions 
  # MINIMUM VALUE is 2 * basisNum - 1
  # NOTE: can play with this for accuracy
  miniGridSize_big <- 4 * basisNum_big
  miniGridSize_small <- 4 * basisNum_small
  
  if (miniGridSize_big >= maxDimension || miniGridSize_small >= minDimension) {
    stop("Warning: Minimum coarse grid based on the number of basis functions is 
         greater than the size of the data. This method is not appropriate here. 
         Either choose less basis functions, or choose a different method, 
         such as exactKronecker or fast. See help file on 
         LKrigNormalizeBasisFFTInterpolate.")
  }
  
  else{
    
    # Creating the actual grid
    # if the first row of basis funcs is larger, the y is the bigger side
    if (gridOrientation == 1){
      gridList<- list( x= seq( bounds[1,1],bounds[2,1],length.out = miniGridSize_big),
                       y= seq( bounds[1,2],bounds[2,2],length.out = miniGridSize_small) )
    }
    # if the second row of basis funcs is larger, the x is the bigger side
    if (gridOrientation == 2){
      gridList<- list( x= seq( bounds[1,1],bounds[2,1],length.out = miniGridSize_small),
                       y= seq( bounds[1,2],bounds[2,2],length.out = miniGridSize_big) )
    }
    
    sGrid<- make.surface.grid(gridList)
    
    # Calling LKrig.cov to evaluate the variance on the coarse grid
    # this is the initial, small variance calculation that we will upsample
    roughMat <- as.surface(sGrid, 
                           LKrig.cov(sGrid, LKinfo = LKinfoNew, marginal=TRUE )
                           )[["z"]]
    
    # FFT step: taking the fft of the small variance
    # reliant on fftwtools package
    fftStep <- fftw2d((roughMat)/length(roughMat))
    
    # Helpful dimensions
    snr <- nrow(fftStep) 
    snc <- ncol(fftStep) 
    
    # Shit parameter
    yShift <- ceiling(((nr/snr) - 1)/2)
    xShift <- ceiling(((nc/snc) - 1)/2)
    
    # Instantiate empty matrix
    temp <- matrix(0, nrow = nr, ncol = nc)
    
    # Helpful for indexing later 
    bigindY <- 1:ceiling(snr/2) 
    bigindX <- 1:ceiling(snc/2) 
    indY <- 1:(snr/2) 
    indX <- 1:(snc/2) 
    # helpful offsets
    bigOffsetY <- (nr - floor(snr/2)) 
    bigOffsetX <- (nc - floor(snc/2))
    smallOffsetY <- (snr - floor(snr/2)) 
    smallOffsetX <- (snc - floor(snc/2))
    
    # Stuffing the small FFT result into the large matrix of zeroes 
    temp[bigindY, bigindX] <- fftStep[bigindY, bigindX] #top left corner
    temp[bigindY, (indX + bigOffsetX)] <- fftStep[bigindY, (indX + smallOffsetX)] #top right corner
    temp[(indY + bigOffsetY), bigindX] <- fftStep[(indY + smallOffsetY), bigindX] #bottom left corner 
    temp[(indY + bigOffsetY), (indX + bigOffsetX)] <- fftStep[(indY + smallOffsetY), (indX + smallOffsetX)] #bottom right corner
    
    # takes the IFFT of the modified big matrix to return our interpolated/upsampled variance
    # again reliant on fftwtools package
    wght <- Re(fftw2d(temp, inverse = 1))
    
    # Shifting due to the periodicity assumed by the FFT
    wght <- LKrig.shift.matrix(wght, shift.row = yShift, shift.col = xShift, periodic = c(TRUE, TRUE))

    #the next steps are done to account for oddly sized/missing data
    # the fft calculation will produce a full grid
    # but only certain observations need a variance associated with them
    
    lat_min <- min(x1[,1]) # minimum latitude
    lat_max <- max(x1[,1]) # maximum latitude
    lon_min <- min(x1[,2]) # minimum longitude
    lon_max <- max(x1[,2]) # maximum longitude
    
    # the grid resolution is known, but if one needs to calculate steps based on nr and nc:
    lat_step <- (lat_max - lat_min) / (nr - 1)
    lon_step <- (lon_max - lon_min) / (nc - 1)
    
    # Map lat-lon to grid indices
    row_indices <- round((x1[,1] - lat_min) / lat_step) + 1
    col_indices <- round((x1[,2] - lon_min) / lon_step) + 1
    
    # Ensure indices are within the bounds and adjust if needed
    row_indices <- pmin(pmax(row_indices, 1), nr)
    col_indices <- pmin(pmax(col_indices, 1), nc)
    
    # Extract the actual values associated with existing data using the indices
    values <- wght[cbind(row_indices, col_indices)]
    
    # vectorize the output matrix to be compatible with LKrig.basis
    wght <- c(values)
  }
  
  return (wght)
  
}
