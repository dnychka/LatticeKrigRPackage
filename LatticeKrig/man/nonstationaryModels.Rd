\name{nonstationaryModels}
\alias{nonstationaryModels}
\title{Specifying non-stationary models}
\description{An overview of specifying non-stationary 
and anisotropic models with some specific examples showing how to vary the alpha, sigma2 and a.wght parameters.}
\details{
The Lattice Krig model can be extended in a natural way to have a  non-stationary covariance that has a multi-scale structure. 

The default and  process model has the form 

 \deqn{ g(x)= \sum_{l=1}^L g_l(x)
}
where g_l(x) are processes defined by fixed basis functions and with
coefficients that follow a mean zero multivariate normal distribution
and with dependence described by a spatial autoregression (SAR).

Under the model with approximate stationarity the SAR is parameterized by
a set of a.wght values that are applied in the same way to every
lattice point and its neighbors. When the process is normalized to
have constant marginal variance, the variance of g_l(x) is given by
sigma2* alpha_l. We recommend that the alpha parameters sum to one and so
the marginal variance of g(x) is given by sigma2. This package allows for
the parameters a.wght, sigma2, and alpha to vary over the spatial
domain. In this way the variance of g_l(x) is given by
sigma2(x)alpha_l(x) and g(x) by sigma2(x) provided alpha_l sums to one. The
a.wght parameters can different values at each lattice point and
possibly at each level. Some preliminary results suggest that having
a.wght vary differently for different levels may be too much
flexibility and may not needed, but examples are still provided. One
way in which such flexibility (non-stationarity and anisotropy across 
multiple resolutions) in both alpha and a.wght could be useful is for 
creating models with non-monotonically decaying covariance. 

We describe the formatting for these features in the LKinfo object below. Essentially they involve specifying one or more of the arguments: 
\code{a.wghtObject}, \code{sigma2.object} or \code{alphaObject}, in the call to 
\code{LKrigSetup}. 

\bold{a.wght}  The general form of this object is as a list of matrices. In this case
length( a.wght) = nlevel and a.wght[[l]] is a matrix where the number of rows is equal to the number of lattice points at level $l$ and in the order that they are indexed for the SAR matrix.  The number of columns depends on the particular geometry but we explain how this works for the rectangular spatial domain, LKRectangle. Here  a.wght[[l]]  has either 1 or 9 columns depending if anisotropy is specified. The weights for the LKRectangle anisotropy case are organized as
\preformatted{
    1 4 7
    2 5 8
    3 6 9
}
So an isotropic model with the center value as 4.5  looks like
\preformatted{
     0   -1    0
    -1   4.5  -1
     0   -1    0
}
and to encoded with the anisotropic extension the  row in 
a.wght[[l]] corresponding to this lattice point would be 

\code{ c( 0, -1, 0, -1, 4.5, -1, 0,-1, 0)}

Specifying the matrices of a.wghts directly can be involved because
one needs to known the lattice information. An easier way to
accomplish specifying these models is to use a function on the spatial
domain to define the a.wght values based on the locations of the
lattice points. In this case one would pass an object, a.wghtObject,
that describes the function. This object needs to be defined so that
predict(a.wghtObject, x) will evaluate the function at the locations
x. (See example below.). With this object, LKrigSetupAwght will
evaluate the function at the lattice points and so create the correct
list of matrices.

\bold{alpha}. To specify spatial varying alpha parameters one
specifies \code{alphaObject} as a list of objects where the predict
function works:

\code{predict(alphaObject[[l]], x) }

The result should give the values for alpha_l(x) with \code{x} a
matrix of arbitrary locations.

One should also set the vector of usual alpha parameters equal to all
 ones so only the spatially varying values as used for the variance.
 E.g. alpha= rep( 1, nlevel)
 

\bold{sigma2} The object \code{sigma2.object} is used define a spatially
varying sigma2(x). Like alpha and a.wght this object must work with the
predict function.

\code{predict(sigma2.object, x)}
}

\examples{

#######################################################
# Define some useful predict functions. 
#######################################################
predict.constantValue<- function(object,x){
  n<- length(object$values)
  m<- nrow( x)
  return( matrix( object$values, nrow=m, ncol=n, byrow=TRUE ) )
}

predict.surfaceGrid<- function(object,x){
  interp.surface( object, x)
}

predict.multivariateSurfaceGrid<- function(object,x){
  dimZ<- dim( object$z)
  L<- dimZ[3]
  out<- matrix( NA, nrow= nrow(x), ncol=L)
  for (  l in 1:L){
    out[,l]<- interp.surface( 
      list( x=object$x,y=object$y, z=object$z[,,l]) , x)
  }
  return( out)
}

predict.multivariateSurfaceGridList <- function(object, x, nLevel, level){
  # object is a multivariateSurfaceGridList
  # this is a list of multivariateSurfaceGrids (each has $x, $y, $z[,,1:9])
  # x: lattice nodes for *one* level (data.frame with columns x,y)
  # nLevel is number of levels
  # level is the current level
  
  # sanity and error checks 
  if (!is.list(object) || length(object) == 0L) {
    stop("multivariateSurfaceGridList must be a non-empty list.")
  }
  if (is.null(level)) {
    stop("predict.multivariateSurfaceGridList requires level index.")
  }
  
  if (!is.null(nLevel) && length(object) != nLevel) {
    stop(sprintf("Mismatch: a.wghtObject has %d levels but LKrig expects %d.",
                 length(object), nLevel))
  }
  
  if (level < 1L || level > length(object)) {
    stop(sprintf("Level index level=%d is out of bounds for list length %d.",
                 level, length(object)))
  }
  
  # call prediction function on the single layer
  predict.multivariateSurfaceGrid(object[[level]], x)
}

################################################
##### Non-stationary examples
##### Below are 10, different examples showing 
##### non-stationary, anisotropic, and 
##### multi-resolution capabilities. 
###############################################
# spatial domain    
sDomain<- rbind( c(-1.2,-1.2),
                 c(1,1))

# we will use this coarse grid to define any 
# surfaces of parameters
# (unrelated to the lattice grids and plotting grid!)
# this is larger than the sDomain to accommodate buffer points
# (with larger ranges when NC is small)
gridList<- list( x = seq( -3, 3,,50),
                 y = seq( -3, 3,,75) )
xPoints<- make.surface.grid( gridList)
fineGrid<- make.surface.grid(
  list( x = seq(-1, 1, ,45),
        y = seq(-1, 1, ,60)
  )
)

##################################################
### end of setup 
#################################################
# Example 1: Spatially varying sigma2 parameter
# sigma2 increases across the domain as a function of first coordinate. 
#################################################
sigma2Values<-  .01 +  10* pnorm( xPoints[,1], mean=.25, sd =.3 )
sigma2Object_ex1<- as.surface( xPoints, sigma2Values) 
class( sigma2Object_ex1)<- "surfaceGrid"

LKinfo_sigma2<- LKrigSetup( sDomain, NC= 4, nlevel = 3,
                            a.wght=4.5, nu=1, sigma2.object=sigma2Object_ex1)   
# Simulate a field from this model
set.seed(123)
simField_sigma2<- LKrig.sim( fineGrid, LKinfo_sigma2)
image.plot( as.surface( fineGrid, simField_sigma2), 
            main = "Spatially Varying Sigma2: Increases Left to Right")
xline( .25, col="grey30")
# see also 
# temp<- as.surface( fineGrid, look)
# I<- 20
# matplot(temp$x, temp$z, type="l", xlab="x", ylab="GP slice" )

######################################################
##### Example 2: Spatially varying alpha parameters 
#######################################################

# The alpha surface at each level will just be the result of 
# bi-linear interpolation of values specified on a small grid.
# To keep things identified, the alpha weights at 
# any grid location are normalized to sum to 1. 
#
# Create a 3 column matrix with (proportional) alpha weights
# at each grid point 
#
alphaTaper_ex2<- pnorm( xPoints[,1], mean = .4, sd=.02)
alphaWeights_ex2<- cbind( alphaTaper_ex2,
                          rep( 1, length( alphaTaper_ex2)), 
                          1-alphaTaper_ex2)
# Normalize to sum to one                        
alphaWeights_ex2 <- alphaWeights_ex2/rowSums(alphaWeights_ex2)

# pack as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list 
# (yes this is a "list of lists")
alphaObject_ex2<- list()
for( k in 1:3){
  alphaSurface<- as.surface( xPoints, alphaWeights_ex2[,k]) 
  class( alphaSurface)<- "surfaceGrid"
  alphaObject_ex2<- c( alphaObject_ex2, list( alphaSurface))
}

# Define the 2-d LatticeKrig model
LKinfo_alpha<- LKrigSetup(sDomain, NC = 4, a.wght=4.5,
                          alpha = c(1,1,1), nlevel = 3, 
                          alphaObject =  alphaObject_ex2 )
# Simulate a field 

set.seed(123)
simField_alpha<- LKrig.sim( fineGrid, LKinfo_alpha)
image.plot( as.surface( fineGrid, simField_alpha),
            main = "Spatially Varying Alpha Parameters (3 Levels)")

######################################################
##### Example 3: Spatially varying a.wght parameters 
##### See above comments and setup
##### for steps that are the same 
#######################################################
awghtTaper_ex3<- pnorm( xPoints[,1] + xPoints[,1],
                        mean = 0, sd=.15)
awghtValues_ex3<- 4.001*awghtTaper_ex3 +  10*(1-awghtTaper_ex3)
# Pack up as a list 
# Convert from a vector to the image/list format: $x $y $z
# Give this object a class so that predict.surfaceGrid is called.

awghtObject_ex3 <- as.surface( xPoints, awghtValues_ex3) 
class( awghtObject_ex3)<- "surfaceGrid"

# Define the 2-d LatticeKrig model

LKinfo_awght<- LKrigSetup(sDomain, NC = 10, NC.buffer=0, 
                          alpha = c(1, .5, .125), nlevel = 3, 
                          a.wghtObject =  awghtObject_ex3)

set.seed(123)            
simField_awght<- LKrig.sim( fineGrid, LKinfo_awght)
image.plot( as.surface( fineGrid, simField_awght),
            main = "Spatially Varying a.wght Parameters")
##############################################
###### Example 4: 1-D example
#############################################
xCoarse_1d<- seq( -.5,1.5,, 40)
yCoarse_1d<-  pnorm( xCoarse_1d, mean=.4, sd=.05)*5 + 2.2 
awghtObject_1d<- Tps(xCoarse_1d, yCoarse_1d, lambda=0)
alphaWeights_1d<-c(.5, .3, .2)
LKinfo_1d<- LKrigSetup( rbind(0,1), NC=10,
                        LKGeometry="LKInterval",
                        nlevel=3, alpha=alphaWeights_1d,
                        a.wghtObject = awghtObject_1d,
                        NC.buffer=2
) 
xFine_1d<- cbind(seq( 0,1,length.out= 200))
set.seed( 123)
simField_1d<- LKrig.sim( xFine_1d, LKinfo_1d, M=5)
matplot( xFine_1d, simField_1d, type="l", lty=1, lwd = 1.5,
         main = "1-D LatticeKrig Simulation (5 Realizations)",
         xlab = "x", ylab = "Simulated Values")
##################################################
######## Example 5: Anisotropy in a.wght
##################################################
#### Stationary anisotropic example
awghtMatrix_aniso<- c( rbind( c(  0,   0, -1.5),
                              c(-.5, 4.5,  -.25),
                              c(-1.5,  0,    0)
)
) 

LKinfo_aniso_stat<- LKrigSetup(sDomain, NC = 5, 
                               a.wght= list( awghtMatrix_aniso), 
                               alpha = c(1, .5, .125), nlevel = 3, 
                               a.wghtObject =  NULL, normalize=TRUE )

set.seed(123)            
simField_aniso_stat<- LKrig.sim( fineGrid, LKinfo_aniso_stat)
image.plot( as.surface( fineGrid, simField_aniso_stat),
            main = "Stationary Anisotropic Covariance")

\dontrun{
  
  #### Example 6: Anisotropy varying over space
  #### First check that the constant model can be reproduced
  awghtMatrix_const<- c( rbind( c(  0,   0, -1.5),
                                c(-.5, 4.5,  -.5),
                                c(-1.5,  0,    0)
  )
  )
  
  awghtObject_const<- list( values=awghtMatrix_const)
  class(awghtObject_const )<- "constantValue"
  
  LKinfo_const<- LKrigSetup(sDomain, NC = 5, 
                            alpha = c(1,.5, .125), nlevel = 3, 
                            a.wghtObject =  awghtObject_const, normalize=TRUE )
  set.seed(123)            
  simField_const<- LKrig.sim( fineGrid, LKinfo_const)
  image.plot( as.surface( fineGrid, simField_const),
              main = "Constant Anisotropic a.wght (Verification)" )            
  
  ###### Example 7: Non-stationary anisotropy
  awghtMatrix_A <- c( rbind( c(    0,   0, -2),
                             c( 0, 4.5,  0),
                             c(-2,  0,     0)
  )
  )
  awghtMatrix_B <- c( rbind( c(  -2,   0,     0),
                             c(  0, 4.5,  0),
                             c(    0,    0, -2)
  )
  )
  # Now create multivariate prediction object.
  gridList_aniso<-  attributes( xPoints)$grid.list
  m1_aniso<- length(gridList_aniso$x)
  m2_aniso<- length(gridList_aniso$y) 
  z_aniso<- array( NA, c( m1_aniso, m2_aniso, 9))
  alphaTaper_aniso<- (xPoints[,1] + 1 )/2
  alphaTaper_aniso<- ifelse( alphaTaper_aniso <= 0, 0, alphaTaper_aniso)
  alphaTaper_aniso<- ifelse( alphaTaper_aniso >= 1, 1, alphaTaper_aniso)
  # Loop over the 9 stencil elements   
  for(j in 1:9) {
    # Linear combination of two a.wght matrices
    zTemp_aniso<- awghtMatrix_A[j] * (1-alphaTaper_aniso) +  
                  awghtMatrix_B[j]*(alphaTaper_aniso)
    # Coerce into an image format    
    z_aniso[,,j]<- as.surface( xPoints, zTemp_aniso)$z
  }
  
  awghtObject_nonstat<- list( x= gridList_aniso$x,  y= gridList_aniso$y, z=z_aniso )
  class( awghtObject_nonstat)<- "multivariateSurfaceGrid"
  
  LKinfo_aniso_nonstat<- LKrigSetup(sDomain, NC = 25, NC.buffer=0,
                                    alpha = c(1,.5), 
                                    nlevel = 2, 
                                    a.wghtObject =  awghtObject_nonstat )
  set.seed(122)  
  fineGrid_ex7<- make.surface.grid(
    list( x = seq(-1, 1, ,150),
          y = seq(-1, 1, ,180)
    )
  )
  simField_aniso_nonstat<- LKrig.sim( fineGrid_ex7, LKinfo_aniso_nonstat)
  image.plot( as.surface( fineGrid_ex7, simField_aniso_nonstat), 
              col = terrain.colors(256),
              main = "Non-Stationary Anisotropic Covariance" )
}




\dontrun{
##################################################
######## Some More Advanced Examples
##################################################

  #####################################################
  ######## Nonstationary awght and anisotropy patterns
  ######## Same for both levels 
  #####################################################
  
  # Big scary function for making param fields 
  # with simple functions/patterns 
  generate_param_field <- function(
    paramtype,
    config,
    sGrid,
    rows,
    columns,
    n){
    
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
  
  ### Some setup for params 
  configs <- c("constant", "taper", "Gaussian", "sinwave", "double_Gaussian")
  n_awghts <- 400
  awghts <- 4 + exp(seq(log(0.001),
                        log(10),
                        length.out=n_awghts))
  low_awghts <- awghts[1:(n_awghts/2)]
  high_awghts <- awghts[(n_awghts/2 + 1):n_awghts]
  
  # Setup for data grid (128x128)
  nrows_advanced <- 128
  ncols_advanced <- nrows_advanced
  n_advanced <- nrows_advanced^2
  
  gridList_advanced<- list( x= seq( -1,1,length.out= nrows_advanced),
                            y= seq( -1,1,length.out= nrows_advanced) )
  sGrid_advanced<- make.surface.grid(gridList_advanced)
  
  # Two layers of basis functions: 20x20 and 39x39
  sidelen_bf_1 <- 20
  sidelen_bf_2 <- 39
  
  
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
    data_rows
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
                                         sidelen, sidelen, n_obj) - 4
    if (awght_override == TRUE){
      kappa2_field <- kappa2_field - kappa2_field + awght_override_val
    }
    
    theta_field <- generate_param_field(paramtype = "theta", 
                                        config = theta_config, 
                                        sGrid_obj, 
                                        sidelen, sidelen, n_obj)
    if (theta_override == TRUE){
      theta_field <- theta_field - theta_field + theta_override_val
    }
    
    rho_field <- generate_param_field(paramtype = "rho", 
                                      config = rho_config, 
                                      sGrid_obj, 
                                      sidelen, sidelen, n_obj)
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
  
  #####################################################
  ######## Example 8: Nonstationary a.wght and anisotropy 
  ######## Same for both levels
  #####################################################
  
  # Create the a.wghts for the first level 
  awghtObj_level1 <- create_awght_obj(
    sidelen = sidelen_bf_1, 
    awght_config = "constant",
    awght_override = TRUE, 
    awght_override_val = 0.01,
    theta_config = "constant", 
    theta_override = TRUE,
    theta_override_val = -pi/4,
    rho_config = "constant",
    rho_override = TRUE,
    rho_override_val = 5,
    plotting = TRUE,
    data_rows = nrows_advanced
  )
  
  # The type of object determines the predict function called
  awghtObject_level1 <- awghtObj_level1$awght
  class( awghtObject_level1)<- "multivariateSurfaceGrid"
  
  LKinfo_ex8<- LKrigSetup(sGrid_advanced, NC = sidelen_bf_1, NC.buffer=0,
                          alpha = c(0.3, 0.7), 
                          nlevel = 2, 
                          a.wghtObject =  awghtObject_level1)
  
  simField_ex8 <- LKrig.sim( sGrid_advanced, LKinfo_ex8)
  imagePlot(as.surface(sGrid_advanced, simField_ex8), col = turbo(256),
            main = "Nonstationary Anisotropy (Same for Both Levels)")
  
  
  #####################################################
  ######## Example 9: Nonstationary a.wght and anisotropy 
  ######## Different for each level
  #####################################################
  
  # Create the a.wghts for the second level 
  awghtObj_level2 <- create_awght_obj(
    sidelen = sidelen_bf_2, 
    awght_config = "constant",
    awght_override = TRUE, 
    awght_override_val = 4.0,
    theta_config = "constant", 
    theta_override = TRUE,
    theta_override_val = pi/4,
    rho_config = "constant",
    rho_override = TRUE,
    rho_override_val = 40,
    plotting = TRUE,
    data_rows = nrows_advanced
  )
  
  awghtObject_level2 <- awghtObj_level2$awght
  class( awghtObject_level2)<- "multivariateSurfaceGrid"
  
  # Pack both multivariateSurfaceGrid objects into
  # a multivariateSurfaceGridList object (required for level-specific parameters)
  awghtObjectList_ex9 <- list(awghtObject_level1, awghtObject_level2)
  class(awghtObjectList_ex9) <- "multivariateSurfaceGridList"
  
  LKinfo_ex9<- LKrigSetup(sGrid_advanced, NC = sidelen_bf_1, NC.buffer=0,
                          alpha = c(0.3, 0.7), 
                          nlevel = 2, 
                          a.wghtObject =  awghtObjectList_ex9)
  
  simField_ex9 <- LKrig.sim( sGrid_advanced, LKinfo_ex9)
  imagePlot(as.surface(sGrid_advanced, simField_ex9), col = turbo(256),
            main = "Nonstationary Anisotropy (Different for Each Level)")
  
  #####################################################
  ######## Example 10: Non-monotonically decaying covariance
  ######## Achieved with: 
  ######## Nonstationary a.wght and anisotropy patterns
  ######## Different for each level
  ######## and nonstationary alpha parameter with 
  ######## modifications to create teleconnection
  #####################################################
  alpha_side_ex10 <- 39
  alphaGrid_ex10 <- list(
    x = seq(-1, 1, length.out = alpha_side_ex10),
    y = seq(-1, 1, length.out = alpha_side_ex10)
  )
  alphaPoints_ex10 <- make.surface.grid(alphaGrid_ex10)  
  
  # Define a taper for alpha
  alphaTaper_ex10 <- pnorm(alphaPoints_ex10[, 1], mean = 0.0, sd = 0.1)  # smooth 0->1 horizontally
  alphaWeights_ex10 <- cbind(1 - alphaTaper_ex10, alphaTaper_ex10)
  
  # Forcing some alpha weights to change for the teleconnection effect
  alphaWeights_ex10[1517:1521, ] <- alphaWeights_ex10[1517:1521, ] + 100
  alphaWeights_ex10[1478:1482, ] <- alphaWeights_ex10[1478:1482, ] + 100
  alphaWeights_ex10[1439:1443, ] <- alphaWeights_ex10[1439:1443, ] + 100
  alphaWeights_ex10[1400:1404, ] <- alphaWeights_ex10[1400:1404, ] + 100
  alphaWeights_ex10[1361:1365, ] <- alphaWeights_ex10[1361:1365, ] + 100
  
  alphaWeights_ex10 <- alphaWeights_ex10 / rowSums(alphaWeights_ex10)  # normalize rows to sum to 1
  
  par(mfrow = c(1,2))
  image.plot(
    x = alphaGrid_ex10$x,
    y = alphaGrid_ex10$y,
    z = matrix(alphaWeights_ex10[, 1], nrow = alpha_side_ex10, ncol = alpha_side_ex10),
    col = viridis(256),
    main = "Level-1 Alpha Weight"
  )
  image.plot(
    x = alphaGrid_ex10$x,
    y = alphaGrid_ex10$y,
    z = matrix(alphaWeights_ex10[, 2], nrow = alpha_side_ex10, ncol = alpha_side_ex10),
    col = viridis(256),
    main = "Level-2 Alpha Weight"
  )
  par(mfrow = c(1,1))
  
  # Pack as a list of two "surfaceGrid" objects so predict.surfaceGrid() is used
  alphaObject_ex10 <- list()
  for (k in 1:2) {
    alphaSurface_ex10 <- as.surface(alphaPoints_ex10, alphaWeights_ex10[, k])  # has $x, $y, $z (matrix)
    class(alphaSurface_ex10) <- "surfaceGrid"
    alphaObject_ex10 <- c(alphaObject_ex10, list(alphaSurface_ex10))
  }
  
  LKinfo_ex10<- LKrigSetup(sGrid_advanced, NC = sidelen_bf_1, NC.buffer=0,
                           alpha = c(1,1),
                           alphaObject = alphaObject_ex10, 
                           nlevel = 2, 
                           a.wghtObject =  awghtObjectList_ex9)
  
  simField_ex10 <- LKrig.sim( sGrid_advanced, LKinfo_ex10)
  imagePlot(as.surface(sGrid_advanced, simField_ex10), col = turbo(256), 
            main = "Teleconnection: Non-Monotonic Covariance")
  
  
  
  #####################################################
  ######## Covariance plotting
  ######## For Examples 8, 9, and 10
  #####################################################
  
  # Function for plotting covariances at a specified location
  plot_cov_surface <- function(
    LKinfo, 
    real_field,
    loc_x, 
    loc_y, 
    plotting = TRUE, 
    cov_contour = TRUE,
    cov_lwd = 1.3
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
  
  # Plot covariance for Example 8
  plot_cov_surface(
    LKinfo = LKinfo_ex8,
    real_field = simField_ex8,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = TRUE,
    cov_lwd = 1.3
  )
  
  # Plot covariance for Example 9
  plot_cov_surface(
    LKinfo = LKinfo_ex9,
    real_field = simField_ex9,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = TRUE,
    cov_lwd = 1.0
  )
  
  # Plot covariance for Example 10
  plot_cov_surface(
    LKinfo = LKinfo_ex10,
    real_field = simField_ex10,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = FALSE
  )

}