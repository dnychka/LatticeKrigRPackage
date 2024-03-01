
# This functions evaluates the variance of the basis functions on a coarser grid, 
# then uses 2D interpolation via FFT in order to smooth/interpolate the variance 
# up to the size of the original grid that we were working with. Should provide a 
# significant computational speedup. 


# NOTE: this function is reliant on the library(fftwtools). Needs to be installed 
#      with the LatticeKrig package itself. 




LKrigNormalizeBasisFFTInterpolate <- function(LKinfo, Level, x1){

  
  # Extracting important information from LKinfo 
  bounds <- LKinfo$x
  basisNum <- LKinfo$latticeInfo$mxDomain[Level,2]
  buffer <- LKinfo$NC.buffer
  alphaNum <- LKinfo$alpha[Level]
  awght <- LKinfo$a.wght[Level]
  
  cat("Doing FFT normalization at level", Level, "with", basisNum, "basis functions", fill = TRUE)
  
  #dimensions of original data, also helpful for shift parameter
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  
  # # even number of rows or columns will result in asymmetric error due to fft assuming periodicity
  # if (nr %% 2 == 0 || nc %% 2 == 0) {
  #   warning("Your input has an even number of rows or columns. 
  #         This will result in asymmetric error pattern, yet similar accuracy. 
  #         Considering using odd numbers for symmetric error.
  #           See LKrigNormalizeBasisFFTInterpolate help file.")
  # }
  
  # Setting up a new LKrig object using the extracted info 
  LKinfoNew <- LKrigSetup(bounds, nlevel = 1, NC = basisNum, NC.buffer = buffer, 
                          alpha = alphaNum, a.wght = awght, normalize = FALSE) 
  
  # Setting a default coarse grid size based on the number of basis functions 
  # MINIMUM VALUE is 2 * basisNum - 1
  # NOTE can play with this for accuracy 
  miniGridSize <- 4 * basisNum
  
  if (miniGridSize >= min(c(nr, nc)) ) {
    stop("Warning: Minimum coarse grid based on the number of basis functions is 
         greater than the size of the data. This method is not appropriate here. 
         Either choose less basis functions, or choose a different method, 
         such as exactKronecker or fast. See help file on 
         LKrigNormalizeBasisFFTInterpolate.")
  }
  
  else{
  
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
    
    
    # Takes the IFFT of the modified big matrix to return a interpolated version of the input matrix
    # Again (reliant on fftwtools package)
    wght <- Re(fftw2d(temp, inverse = 1))
    
    # Shifting due to the periodicity assumed by the FFT
    wght <- LKrig.shift.matrix(wght, shift.row = yShift, shift.col = xShift, periodic = c(TRUE, TRUE))
    
    
    # # correcting to have symmetric error for grids with even dimensions
    # if (ncol(wght) %% 2 == 0){
    #   corner_size <- nc/2
    # 
    #   # extract bottom right corner (will show up as top right in imagePlot)
    #   bottom_right_corner <- wght[(corner_size+1):nc, (corner_size+1):nc]
    # 
    #   # create empty matrix to fill in
    #   expanded_matrix <- matrix(0, nrow = 2*corner_size, ncol = 2*corner_size)
    # 
    #   # fill in the new matrix by rotating the bottom corner around (symmetric in both x and y)
    #   # can also be thought of as mirroring it across the x and then mirroring the result of that across y
    #   expanded_matrix[(corner_size+1):(2*corner_size), (corner_size+1):(2*corner_size)] <- bottom_right_corner  # Bottom right
    #   expanded_matrix[1:corner_size, (corner_size+1):(2*corner_size)] <- bottom_right_corner[corner_size:1, ]  # Top right
    #   expanded_matrix[(corner_size+1):(2*corner_size), 1:corner_size] <- bottom_right_corner[, corner_size:1]  # Bottom left
    #   expanded_matrix[1:corner_size, 1:corner_size] <- bottom_right_corner[corner_size:1, corner_size:1]  # Top left
    # 
    #   # reassigning wght as the new matrix
    #   wght <- expanded_matrix
    # }
    # 
    # #error will already be symmetric for odd dimensions
    # if (ncol(wght) %% 2 == 1 || nrow(wght) %% 2 == 1){
    #   wght <- wght
    # }
    
    #vectorize the variance matrix to be compatible with LKrig.basis
  
    lat_min <- min(x1[,1]) # Define your minimum latitude
    lat_max <- max(x1[,1]) # Define your maximum latitude
    lon_min <- min(x1[,2]) # Define your minimum longitude
    lon_max <- max(x1[,2]) # Define your maximum longitude
    
    # The grid resolution is known, but if you need to calculate steps based on nr and nc:
    lat_step <- (lat_max - lat_min) / (nr - 1)
    lon_step <- (lon_max - lon_min) / (nc - 1)
    
    # Map lat-lon to grid indices
    row_indices <- round((x1[,1] - lat_min) / lat_step) + 1
    col_indices <- round((x1[,2] - lon_min) / lon_step) + 1
    
    # Ensure indices are within the bounds of `wght` and adjust if needed
    row_indices <- pmin(pmax(row_indices, 1), nr)
    col_indices <- pmin(pmax(col_indices, 1), nc)
    
    # Extract the values using the indices
    values <- wght[cbind(row_indices, col_indices)]
    
    wght <- c(values)
  }
  
  return (wght)
  
}



