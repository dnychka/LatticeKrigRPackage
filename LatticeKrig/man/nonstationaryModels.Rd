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

######################################################
##### This is an extended example showing how to define 
##### spatially varying sigma2 parameter 
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
# sigma2 increases across the domain as a function of first coordinate. 
sigma2Temp<-  .01 +  10* pnorm( xPoints[,1], mean=.25, sd =.3 )
sigma2.object<- as.surface( xPoints, sigma2Temp) 
class( sigma2.object)<- "surfaceGrid"

LKinfo<- LKrigSetup( sDomain, NC= 4, nlevel = 3,
                     a.wght=4.5, nu=1, sigma2.object=sigma2.object)   
# simulate a field from this model
set.seed(123)
look<- LKrig.sim( fineGrid, LKinfo)
image.plot( as.surface( fineGrid, look))
xline( .25, col="grey30")
# see also 
# temp<- as.surface( fineGrid, look)
# I<- 20
# matplot(temp$x, temp$z, type="l", xlab="x", ylab="GP slice" )

######################################################
##### spatially varying alpha parameters 
#######################################################

# the alpha surface at each level will just be the result of 
# bi-linear interpolation of values specified on a small grid.
# To keep things identified the alpha weights at 
# any grid location are 
# normalized to sum to 1. 
#
# create a 3 column matrix with  (proportional) alpha weights
# at each grid point 
#
taper<- pnorm( xPoints[,1], mean = .4, sd=.02)
alphaTemp<- cbind( taper,
                   rep( 1, length( taper)), 
                   1-taper)
# normalize to sum to one                        
alphaTemp <- alphaTemp/rowSums(alphaTemp)

# pack as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list 
# (yes this is a "list of lists")
alphaObject<- list()
for( k in 1:3){
  hold<- as.surface( xPoints, alphaTemp[,k]) 
  class( hold)<- "surfaceGrid"
  alphaObject<- c( alphaObject, list( hold))
}

# define the 2-d LatticeKrig model
LKinfo<- LKrigSetup(sDomain, NC = 4, a.wght=4.5,
                    alpha = c(1,1,1), nlevel = 3, 
                    alphaObject =  alphaObject )
# simulate a field 

set.seed(123)
look<- LKrig.sim( fineGrid, LKinfo)
image.plot( as.surface( fineGrid, look))

######################################################
##### spatially varying a.wght parameters 
##### See above comments and setup
##### for steps that are the same 
#######################################################
taper<- pnorm( xPoints[,1] + xPoints[,1],
               mean = 0, sd=.15)
a.wghtTemp<- 4.001*taper +  10*(1-taper)
# pack up as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list (yes this is 
# a "list of lists")

a.wghtObjectA <- as.surface( xPoints, a.wghtTemp) 
class( a.wghtObjectA)<- "surfaceGrid"


# define the 2-d LatticeKrig model

LKinfo2<- LKrigSetup(sDomain, NC = 10, NC.buffer=0, 
                     alpha = c(1, .5, .125), nlevel = 3, 
                     a.wghtObject =  a.wghtObjectA)

set.seed(123)            
look<- LKrig.sim( fineGrid, LKinfo2)
image.plot( as.surface( fineGrid, look))
##############################################
###### 1-d example
#############################################
xCoarse1<- seq( -.5,1.5,, 40)
y<-  pnorm( xCoarse1, mean=.4, sd=.05)*5 + 2.2 
a.wghtObject<- Tps(xCoarse1, y, lambda=0)
alphaTemp<-c(.5, .3, .2)
LKinfoTEST<- LKrigSetup( rbind(0,1), NC=10,
                         LKGeometry="LKInterval",
                         nlevel=3, alpha=alphaTemp,
                         a.wghtObject = a.wghtObject,
                         NC.buffer=2
) 
xFine1<- cbind(seq( 0,1,length.out= 200))
set.seed( 123)
look<- LKrig.sim( xFine1, LKinfoTEST, M=5)
matplot( xFine1, look, type="l", lty=1)
##################################################
######## Anisotropy in a.wght
##################################################
#### stationary example
a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                     c(-.5, 4.5,  -.25),
                     c(-1.5,  0,    0)
)
) 

LKinfo3<- LKrigSetup(sDomain, NC = 5, 
                     a.wght= list( a.wghtM2), 
                     alpha = c(1, .5, .125), nlevel = 3, 
                     a.wghtObject =  NULL, normalize=TRUE )


set.seed(123)            
look<- LKrig.sim( fineGrid, LKinfo3)
image.plot( as.surface( fineGrid, look))

\dontrun{
  
  #### Anisotropy varying over space
  #### First check that the constant model can be reproduced
  a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                       c(-.5, 4.5,  -.5),
                       c(-1.5,  0,    0)
  )
  )
  
  a.wghtObject<- list( values=a.wghtM2)
  class(a.wghtObject )<- "constantValue"
  
  LKinfo4<- LKrigSetup(sDomain, NC = 5, 
                       alpha = c(1,.5, .125), nlevel = 3, 
                       a.wghtObject =  a.wghtObject, normalize=TRUE )
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKinfo4)
  image.plot( as.surface( fineGrid, look) )            
  
  ###### non-stationary anisotropy
  a.wghtA <- c( rbind( c(    0,   0, -2),
                       c( 0, 4.5,  0),
                       c(-2,  0,     0)
  )
  )
  a.wghtB <- c( rbind( c(  -2,   0,     0),
                       c(  0, 4.5,  0),
                       c(    0,    0, -2)
  )
  )
  # Now create multivariate prediction object.
  gridList<-  attributes( xPoints)$grid.list
  m1<- length(gridList$x)
  m2<- length(gridList$y) 
  z<- array( NA, c( m1,m2,9))
  alpha<- (xPoints[,1] + 1 )/2
  alpha<- ifelse( alpha <= 0, 0, alpha)
  alpha<- ifelse( alpha >= 1, 1, alpha)
  # A for loop over the 9 pixels   
  for(j in 1:9) {
    # linear combination of two a.wght matrices
    # 
    zTemp<- a.wghtA[j] * (1-alpha) +  a.wghtB[j]*(alpha)
    # coerce into an image format    
    z[,,j]<- as.surface( xPoints, zTemp)$z
  }
  
  a.wghtObject<- list( x= gridList$x,  y= gridList$y, z=z )
  class( a.wghtObject)<- "multivariateSurfaceGrid"
  
  LKinfo5<- LKrigSetup(sDomain, NC = 25, NC.buffer=0,
                       alpha = c(1,.5), 
                       nlevel = 2, 
                       a.wghtObject =  a.wghtObject )
  set.seed(122)  
  fineGrid<- make.surface.grid(
    list( x = seq(-1, 1, ,150),
          y = seq(-1, 1, ,180)
    )
  )
  look<- LKrig.sim( fineGrid, LKinfo5)
  image.plot( as.surface( fineGrid, look), col = terrain.colors(256) )
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
  
  # Some setup for data (lives on 128x128 grid)
  rows <- 128
  columns <- rows
  n <- rows^2
  
  gridList<- list( x= seq( -1,1,length.out= rows),
                   y= seq( -1,1,length.out= rows) )
  sGrid<- make.surface.grid(gridList)
  
  # two layers of basis functions, 20x20 and 40x40
  sidelen_bf_1 <- 20
  sidelen_bf_2 <- 39
  
  
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
    # setup grid 
    gridList <- list(x= seq( -1,1,length.out= sidelen),
                     y= seq( -1,1,length.out= sidelen))
    sGrid <- make.surface.grid(gridList)
    n <- sidelen^2
    
    # make the actual params 
    kappa2 <- generate_param_field(paramtype = "awght", 
                                   config = awght_config, 
                                   sGrid, 
                                   sidelen, sidelen, n) - 4
    if (awght_override == TRUE){
      kappa2 <- kappa2 - kappa2 + awght_override_val
    }
    
    theta <- generate_param_field(paramtype = "theta", 
                                  config = theta_config, 
                                  sGrid, 
                                  sidelen, sidelen, n)
    if (theta_override == TRUE){
      theta <- theta - theta + theta_override_val
    }
    
    rho <- generate_param_field(paramtype = "rho", 
                                 config = rho_config, 
                                 sGrid, 
                                 sidelen, sidelen, n)
    if (rho_override == TRUE){
      rho <- rho - rho + rho_override_val
    }
    rhox <- sqrt(rho)
    rhoy <- 1/rhox
    
    # plotting functionality 
    if (plotting == TRUE){
      par(mfrow=c(1,3))
      image.plot(kappa2, main = "kappa2 surface", col = viridis(256))
      image.plot(theta, main = "theta surface", col = viridis(256))
      image.plot(rho, main = "rho surface", col = viridis(256))
      par(mfrow = c(1,1))
    }
    
    # create H tensor
    H11 <- ( rhox^2 * (cos(theta))^2) + ( rhoy^2 * (sin(theta))^2 ) 
    H12 <- (rhoy^2 - rhox^2)*(sin(theta)*cos(theta))
    H21 <- H12 
    H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
    
    # create stencil tensor 
    stencil_tensor <- array( NA, c( sidelen,sidelen,9))
    stencil_tensor[,,1] <- 0.5 * H12
    stencil_tensor[,,2] <- -H22
    stencil_tensor[,,3] <- -0.5 * H12
    stencil_tensor[,,4] <- -H11
    stencil_tensor[,,5] <- kappa2 + 2 * H11 + 2 * H22
    stencil_tensor[,,6] <- -H11
    stencil_tensor[,,7] <- -0.5 * H12
    stencil_tensor[,,8] <- -H22
    stencil_tensor[,,9] <- 0.5 * H12
    
    # plot resulting stencil tensor 
    if (plotting == TRUE){
      par(mfrow = c(3,3))
      image.plot(stencil_tensor[,,1], main = "top left corner", col = viridis(256))
      image.plot(stencil_tensor[,,2], main = "top middle", col = viridis(256))
      image.plot(stencil_tensor[,,3], main = "top right corner", col = viridis(256))
      image.plot(stencil_tensor[,,4], main = "middle left", col = viridis(256))
      image.plot(stencil_tensor[,,5], main = "middle (awght)", col = viridis(256))
      image.plot(stencil_tensor[,,6], main = "middle right", col = viridis(256))
      image.plot(stencil_tensor[,,7], main = "bottom left corner", col = viridis(256))
      image.plot(stencil_tensor[,,8], main = "bottom middle", col = viridis(256))
      image.plot(stencil_tensor[,,9], main = "bottom right corner", col = viridis(256))
      par(mfrow = c(1,1))
    }
    
    # create awght object 
    awght_obj_final <- list( x= gridList$x,  y= gridList$y, z=stencil_tensor )
    class( awght_obj_final)<- "multivariateSurfaceGrid"
    
    # return an object with everything in it 
    object <- c()
    object$awght <- awght_obj_final
    object$kappa2 <- kappa2
    object$theta <- theta
    object$rhox <- rhox
    object$rhoy <- rhoy
    
    return(object)
  }
  
  # create the awghts for the first level 
  awght_1 <- create_awght_obj(
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
    data_rows = rows
  )
  
  # the type of object determines the predict function called
  a.wghtObject_1 <- awght_1$awght
  class( a.wghtObject_1)<- "multivariateSurfaceGrid"
  
  LKinfo6<- LKrigSetup(sGrid, NC = sidelen_bf_1, NC.buffer=0,
                       alpha = c(0.3, 0.7), 
                       nlevel = 2, 
                       a.wghtObject =  a.wghtObject_1)
  
  look6 <- LKrig.sim( sGrid, LKinfo6)
  imagePlot(as.surface(sGrid, look6), col = turbo(256))
  
  
  #####################################################
  ######## Nonstationary awght and anisotropy patterns
  ######## Different for each level
  #####################################################
  
  # create the awghts for the second level 
  awght_2 <- create_awght_obj(
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
    data_rows = rows
  )
  
  a.wghtObject_2 <- awght_2$awght
  class( a.wghtObject_2)<- "multivariateSurfaceGrid"
  
  # pack both multivariateSurfaceGrid objects into
  # a multivariateSurfaceGridList object. this is important
  awghtObjectList <- list(a.wghtObject_1, a.wghtObject_2)
  class(awghtObjectList) <- "multivariateSurfaceGridList"
  
  LKinfo7<- LKrigSetup(sGrid, NC = sidelen_bf_1, NC.buffer=0,
                       alpha = c(0.3, 0.7), 
                       nlevel = 2, 
                       a.wghtObject =  awghtObjectList)
  
  look7 <- LKrig.sim( sGrid, LKinfo7)
  imagePlot(as.surface(sGrid, look7), col = turbo(256))
  
  #####################################################
  ######## Non-monotonically decaying covariance
  ######## Achieved with: 
  ######## Nonstationary awght and anisotropy patterns
  ######## Different for each level
  ######## and nonstationary alpha parameter with 
  ######## some minor modifications 
  #####################################################
  alpha_side <- 39
  alpha_grid <- list(
    x = seq(-1, 1, length.out = alpha_side),
    y = seq(-1, 1, length.out = alpha_side)
  )
  alpha_pts <- make.surface.grid(alpha_grid)  
  
  # define a taper for alpha
  taper <- pnorm(alpha_pts[, 1], mean = 0.0, sd = 0.1)  # smooth 0->1 horizontally
  alphaTemp <- cbind(1 - taper, taper)
  
  # forcing some alpha weights to change for the teleconnection
  alphaTemp[1517:1521] <- alphaTemp[1517:1521] + 100
  alphaTemp[1478:1482] <- alphaTemp[1478:1482] + 100
  alphaTemp[1439:1443] <- alphaTemp[1439:1443] + 100
  alphaTemp[1400:1404] <- alphaTemp[1400:1404] + 100
  alphaTemp[1361:1365] <- alphaTemp[1361:1365] + 100
  
  alphaTemp <- alphaTemp / rowSums(alphaTemp)           # normalize rows to sum to 1
  
  par(mfrow = c(1,2))
  image.plot(
    x = alpha_grid$x,
    y = alpha_grid$y,
    z = matrix(alphaTemp[, 1], nrow = alpha_side, ncol = alpha_side),
    col = viridis(256),
    main = "Level-1 alpha weight"
  )
  image.plot(
    x = alpha_grid$x,
    y = alpha_grid$y,
    z = matrix(alphaTemp[, 2], nrow = alpha_side, ncol = alpha_side),
    col = viridis(256),
    main = "Level-2 alpha weight"
  )
  par(mfrow = c(1,1))
  
  # Pack as a list of two "surfaceGrid" objects so predict.surfaceGrid() is used.
  alphaObject <- list()
  for (k in 1:2) {
    surf <- as.surface(alpha_pts, alphaTemp[, k])  # has $x, $y, $z (matrix)
    class(surf) <- "surfaceGrid"
    alphaObject <- c(alphaObject, list(surf))
  }
  
  LKinfo8<- LKrigSetup(sGrid, NC = sidelen_bf_1, NC.buffer=0,
                       alpha = c(1,1),
                       alphaObject = alphaObject, 
                       nlevel = 2, 
                       a.wghtObject =  awghtObjectList)
  
  look8 <- LKrig.sim( sGrid, LKinfo8)
  imagePlot(as.surface(sGrid, look8), col = turbo(256), 
            main = "Top right teleconnection")
  
  
  
  #####################################################
  ######## Covariance plotting
  ######## For the above three examples
  #####################################################
  
  # additional function for plotting covariances at a certain point
  plot_cov_surface <- function(
    LKinfo, 
    real_field,
    loc_x, 
    loc_y, 
    plotting = TRUE, 
    cov_countour = TRUE,
    cov_lwd = 1.3
  ){
    cov_surface <- LKrig.cov(sGrid, rbind(c(loc_x,loc_y)), LKinfo)
    cov_surface <- matrix(cov_surface, nrow = rows, ncol = rows)
    
    if (plotting == TRUE){
      par(mfrow = c(1,2), mar=  c(5.1, 4.1, 4.1, 2.1))
      image.plot(as.surface(sGrid, cov_surface), col = viridis(256), 
                 main = "Covariance")
      if (cov_countour == TRUE){
        contour( 
          as.surface(sGrid, cov_surface), 
          col = "white", 
          add = TRUE, lwd = cov_lwd
        )
      }
      
      imagePlot(as.surface(sGrid, real_field), 
                main = "Spatial Field", col = turbo(256))
      if (cov_countour == TRUE){
        contour( 
          as.surface(sGrid, cov_surface), 
          col = "black", 
          add = TRUE, lwd = cov_lwd
        )
      }
      par(mfrow = c(1,1))
    }
    # return (cov_surface)
  }
  
  plot_cov_surface(
    LKinfo = LKinfo6,
    real_field = look6,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_countour = TRUE,
    cov_lwd = 1.3
  )
  plot_cov_surface(
    LKinfo = LKinfo7,
    real_field = look7,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_countour = TRUE,
    cov_lwd = 1.0
  )
  plot_cov_surface(
    LKinfo = LKinfo8,
    real_field = look8,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_countour = FALSE,
  )

}