# Function for plotting covariances at a specified location
plot_cov_surface <- function(
    LKinfo, 
    real_field,
    loc_x, 
    loc_y, 
    plotting = TRUE, 
    cov_contour = TRUE,
    cov_lwd = 1.3,
    nrows_advanced,
    sGrid_advanced
){
  cov_surface_calc <- LKrig.cov(sGrid_advanced, rbind(c(loc_x, loc_y)), LKinfo)
  cov_surface_calc <- matrix(cov_surface_calc, nrow = nrows_advanced, ncol = nrows_advanced)
  
  if (plotting == TRUE){
    par(mfrow = c(1,2), mar =  c(5.1, 4.1, 4.1, 2.1))
    image.plot(as.surface(sGrid_advanced, cov_surface_calc), col = viridis(256), 
               main = "Covariance Function")
    if (cov_contour == TRUE){
      contour( 
        as.surface(sGrid_advanced, cov_surface_calc), 
        col = "white", 
        add = TRUE, lwd = cov_lwd
      )
    }
    
    imagePlot(as.surface(sGrid_advanced, real_field), 
              main = "Simulated Spatial Field", col = turbo(256))
    if (cov_contour == TRUE){
      contour( 
        as.surface(sGrid_advanced, cov_surface_calc), 
        col = "black", 
        add = TRUE, lwd = cov_lwd
      )
    }
    par(mfrow = c(1,1))
  }
  invisible(cov_surface_calc)
}