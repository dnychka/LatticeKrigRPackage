# Function to create a.wght objects with anisotropy parameters
create_awght_obj <- function(
    sidelen, 
    awght_config,
    awght_override = FALSE,
    awght_override_val,
    theta_config, 
    theta_override = FALSE,
    theta_override_val,
    rho_config,
    rho_override = FALSE,
    rho_override_val,
    plotting = FALSE,
    data_rows,
    #### added because no visible binding DWN
    awghts,
    low_awghts,
    high_awghts
){
  # Setup grid 
  gridList_obj <- list(x= seq( -1,1,length.out= sidelen),
                       y= seq( -1,1,length.out= sidelen))
  sGrid_obj <- make.surface.grid(gridList_obj)
  n_obj <- sidelen^2
  
  # Generate the parameter fields
  kappa2_field <- generate_param_field(paramtype = "awght", 
                                       config = awght_config, 
                                       sGrid_obj, 
                                       sidelen, sidelen, n_obj,
                                       #### added because no visible binding DWN
                                       awghts,
                                       low_awghts,
                                       high_awghts) - 4
  if (awght_override == TRUE){
    kappa2_field <- kappa2_field - kappa2_field + awght_override_val
  }
  
  theta_field <- generate_param_field(paramtype = "theta", 
                                      config = theta_config, 
                                      sGrid_obj, 
                                      sidelen, sidelen, n_obj, #### added because no visible binding DWN
                                      awghts,
                                      low_awghts,
                                      high_awghts)
  if (theta_override == TRUE){
    theta_field <- theta_field - theta_field + theta_override_val
  }
  
  rho_field <- generate_param_field(paramtype = "rho", 
                                    config = rho_config, 
                                    sGrid_obj, 
                                    sidelen, sidelen, n_obj,
                                    #### added because no visible binding DWN
                                    awghts,
                                    low_awghts,
                                    high_awghts)
  if (rho_override == TRUE){
    rho_field <- rho_field - rho_field + rho_override_val
  }
  rhox_field <- sqrt(rho_field)
  rhoy_field <- 1/rhox_field
  
  # Plotting functionality 
  if (plotting == TRUE){
    par(mfrow=c(1,3))
    image.plot(kappa2_field, main = "Kappa2 Surface", col = viridis(256))
    image.plot(theta_field, main = "Theta Surface", col = viridis(256))
    image.plot(rho_field, main = "Rho Surface", col = viridis(256))
    par(mfrow = c(1,1))
  }
  
  # Create H tensor (anisotropy tensor)
  H11_tensor <- ( rhox_field^2 * (cos(theta_field))^2) + 
    ( rhoy_field^2 * (sin(theta_field))^2 ) 
  H12_tensor <- (rhoy_field^2 - rhox_field^2)*(sin(theta_field)*cos(theta_field))
  H21_tensor <- H12_tensor 
  H22_tensor <- (rhox_field^2 * (sin(theta_field))^2) + 
    (rhoy_field^2 * (cos(theta_field))^2)
  
  # Create stencil tensor (9-point stencil for each grid location)
  stencil_tensor_obj <- array( NA, c( sidelen, sidelen, 9))
  stencil_tensor_obj[,,1] <- 0.5 * H12_tensor
  stencil_tensor_obj[,,2] <- -H22_tensor
  stencil_tensor_obj[,,3] <- -0.5 * H12_tensor
  stencil_tensor_obj[,,4] <- -H11_tensor
  stencil_tensor_obj[,,5] <- kappa2_field + 2 * H11_tensor + 2 * H22_tensor
  stencil_tensor_obj[,,6] <- -H11_tensor
  stencil_tensor_obj[,,7] <- -0.5 * H12_tensor
  stencil_tensor_obj[,,8] <- -H22_tensor
  stencil_tensor_obj[,,9] <- 0.5 * H12_tensor
  
  # Plot resulting stencil tensor 
  if (plotting == TRUE){
    par(mfrow = c(3,3))
    image.plot(stencil_tensor_obj[,,1], main = "Top Left Corner", col = viridis(256))
    image.plot(stencil_tensor_obj[,,2], main = "Top Middle", col = viridis(256))
    image.plot(stencil_tensor_obj[,,3], main = "Top Right Corner", col = viridis(256))
    image.plot(stencil_tensor_obj[,,4], main = "Middle Left", col = viridis(256))
    image.plot(stencil_tensor_obj[,,5], main = "Middle (a.wght)", col = viridis(256))
    image.plot(stencil_tensor_obj[,,6], main = "Middle Right", col = viridis(256))
    image.plot(stencil_tensor_obj[,,7], main = "Bottom Left Corner", col = viridis(256))
    image.plot(stencil_tensor_obj[,,8], main = "Bottom Middle", col = viridis(256))
    image.plot(stencil_tensor_obj[,,9], main = "Bottom Right Corner", col = viridis(256))
    par(mfrow = c(1,1))
  }
  
  # Create a.wght object 
  awght_obj_final <- list( x= gridList_obj$x,  y= gridList_obj$y, z=stencil_tensor_obj )
  class( awght_obj_final)<- "multivariateSurfaceGrid"
  
  # Return an object with all components
  result_obj <- list()
  result_obj$awght <- awght_obj_final
  result_obj$kappa2 <- kappa2_field
  result_obj$theta <- theta_field
  result_obj$rhox <- rhox_field
  result_obj$rhoy <- rhoy_field
  
  return(result_obj)
}
