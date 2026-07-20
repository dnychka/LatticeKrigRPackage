

generate_param_field <- function(
  # Big scary function for making param fields 
  # with simple functions/patterns 
    paramtype,
    config,
    sGrid,
    rows,
    columns,
    n,
    #### added because no visible binding DWN
    awghts,
    low_awghts,
    high_awghts
    ){
  
  if (paramtype == "awght"){
    # sampling from predetermined awghts
    # (this can be changed, just a design choice)
    lower_bound <- min(awghts)
    upper_bound <- max(awghts)
    midpoint <- (lower_bound + upper_bound)/2 # 5
    
    param_constant <- sample(awghts, 1)
    param_low <- sample(low_awghts, 1)
    param_high <- sample(high_awghts, 1)
    
    # these variables are dependent on the domain of the param (awght)
    gauss_amps <- runif(2, 0.1, 0.5) * ifelse(param_constant <= midpoint, 1, -1)
    mult_factor <- runif(1, 0.001, 0.1998)
    
  } else if (paramtype == "rho"){
    # rho is constructed simply based on bounds
    lower_bound <- 1
    upper_bound <- 7
    midpoint <- (lower_bound + upper_bound) / 2  # 4
    
    param_constant <- runif(1, lower_bound, upper_bound)
    param_low <- runif(1, lower_bound, param_constant)
    param_high <- runif(1, param_constant, upper_bound)
    
    # dependent on rho domain
    gauss_amps <- runif(2, 0.1, 1.5) * ifelse(param_constant <= midpoint, 1, -1)
    mult_factor <- runif(1, 0.001, 0.7498)
    
  } else if (paramtype == "theta"){
    # making sure to obtain all possible ellipses
    lower_bound <- 0
    upper_bound <- 3
    midpoint <- (lower_bound + upper_bound) / 2  # 1.5 (roughly pi/2)
    
    param_constant <- runif(1, lower_bound, upper_bound)
    param_low <- runif(1, lower_bound, param_constant)
    param_high <- runif(1, param_constant, upper_bound)
    
    # dependent on theta domain
    gauss_amps <- runif(2, 0.1, pi/4) * ifelse(param_constant < midpoint, 1, -1)
    mult_factor <- runif(1, 0.001, 0.9998)
    
  } else {
    stop("Unknown paramtype. Please use either awght, rho, or theta.")
  }
  
  # these quantities dont depend on which param field we are making
  taper_sd <- runif(1, 0.05, 1)
  
  gauss_slopes <- runif(2, 0.2, 0.5)
  # gauss locations need to be within [-1,1],[-1,1] domain
  gauss_locs <- runif(4, sGrid[1], tail(sGrid,1)[1])
  
  coast_sharpnesses <- runif(2, 3, 50)
  coast_bump_scales <- runif(2, 0.1, 0.5)
  coast_freqs <- runif(2, 0.4, 3)
  coast_coefs <- runif(2, -2 , 2)
  # make sure that adding coastlines doesnt go out of bounds
  coast_amp1 <- runif(1,0.1 , 0.9)
  coast_amp2 <- runif(1, 0.1, 1-coast_amp1)
  
  # make sure sin amplitude doesnt go out of bounds
  if (param_constant > midpoint){
    sin_amp <- runif(1, 0, upper_bound - param_constant)
  }
  else{
    sin_amp <- runif(1, 0, param_constant - lower_bound)
  }
  sin_freq <- runif(1, 1.5, 5)
  sin_orientation <- sample(c("horiz","vert"),1)
  
  #number of basis for generating a gp for the param field.
  # a gp will then be generated from this gp.
  #Just like inception...
  num_basis_param <- sample(c(6:32), 1)
  
  if (config == "constant"){
    param_func <- rep(param_constant, n)
  }
  
  else if (config == "taper"){
    taper<- pnorm( sGrid[,1] + sGrid[,1],
                   mean = 0, sd = taper_sd)
    param_func<- param_low*taper +  param_high*(1-taper)
  }
  
  else if (config == "Gaussian"){
    param_func <- param_constant +
      gauss_amps[1] * exp(-((sGrid[,1]^2 + sGrid[,2]^2) / gauss_slopes[1]))
  }
  
  else if (config == "sinwave"){
    if (sin_orientation == "vert"){
      param_func <- param_constant + sin_amp * sin( pi * sGrid[,1] * sin_freq)
    }
    else{
      param_func <- param_constant + sin_amp * cos( pi * sGrid[,2] * sin_freq)
    }
  }
  
  else if (config == "double_Gaussian"){
    
    peak1 <- gauss_amps[1] * exp(-((sGrid[,1] - gauss_locs[1])^2 +
                                     (sGrid[,2] - gauss_locs[2])^2) / gauss_slopes[1])
    peak2 <- gauss_amps[2] * exp(-((sGrid[,1] + gauss_locs[3])^2 +
                                     (sGrid[,2] + gauss_locs[4])^2) / gauss_slopes[2])
    param_func <- param_constant + peak1 + peak2
  }
  
  # finally, reshape back to a matrix and return
  param_field <- matrix(param_func, nrow = rows, ncol = columns)
  return(param_field)
}