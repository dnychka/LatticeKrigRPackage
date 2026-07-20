%%BEGIN HEADER
%
% LatticeKrig is a package for analysis of spatial data written for
% the R software environment.
% Copyright (C) 2026 Colorado School of Mines
% 1500 Illinois St., Golden, CO 80401
% Contact: Douglas Nychka,  douglasnychka@gmail.com,
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% A copy of the GNU General Public License is included
% along with the R software environment if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
% or refer to  http://www.r-project.org/Licenses/GPL-2
%
%%END HEADER

\name{nonstationaryModels}
\alias{nonstationaryModels}
\title{Specifying non-stationary models}
\description{An overview of specifying non-stationary 
and anisotropic models with some specific examples showing how to vary the 
alpha, sigma2 and a.wght parameters.}
\details{
The Lattice Krig model can be extended in a natural way to have a  
non-stationary covariance that has a multi-scale structure. 

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

We describe the formatting for these features in the LKinfo object below.
Essentially they involve specifying one or more of the arguments: 
\code{a.wghtObject}, \code{sigma2.object} or \code{alphaObject}, in the call to 
\code{LKrigSetup}. 

\bold{a.wght}  The general form of this object is as a list of matrices.
In this case
length( a.wght) = nlevel and a.wght[[l]] is a matrix where the number of rows
is equal to the number of lattice points at level $l$ and in the order that
they are indexed for the SAR matrix.  The number of columns depends on the 
particular geometry but we explain how this works for the rectangular spatial
domain, LKRectangle. Here  a.wght[[l]]  has either 1 or 9 columns depending 
if anisotropy is specified. The weights for the LKRectangle anisotropy case 
are organized as

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
and unrolled in the form  of the anisotropic format the  row in 
a.wght[[l]] corresponding to this lattice point would be 

\code{ c( 0, -1, 0, -1, 4.5, -1, 0,-1, 0)}.

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
 ones so only the spatially varying values are used for the variance.
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
    cat("Mismatch: a.wghtObject has", length(object), 
                "levels but LKrig expects", nLevel, fill = TRUE)
    stop()
  }
  
  if (level < 1L || level > length(object)) {
    cat("Level index level=", level, "is out of bounds for list length", 
                length(object), fill = TRUE)
    stop()
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
\dontrun{
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
}
\dontrun{
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
}

\dontrun{
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

}

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
  
  awghtObject_nonstat<- list( x= gridList_aniso$x,  y= gridList_aniso$y,
  z=z_aniso )
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





##################################################
######## Some More Advanced Examples
##################################################

  #####################################################
  ######## Nonstationary awght and anisotropy patterns
  ######## Same for both levels 
  #####################################################
  
  
\dontrun{  
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
    data_rows = nrows_advanced,
    awghts = awghts,
    low_awghts=low_awghts,
    high_awghts=high_awghts
    
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
 }
 
  
  #####################################################
  ######## Example 9: Nonstationary a.wght and anisotropy 
  ######## Different for each level
  #####################################################
  \dontrun{
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
    data_rows = nrows_advanced,
    awghts = awghts,
    low_awghts=low_awghts,
    high_awghts=high_awghts
  )
  
  awghtObject_level2 <- awghtObj_level2$awght
  class( awghtObject_level2)<- "multivariateSurfaceGrid"
  
  # Pack both multivariateSurfaceGrid objects into
  # multivariateSurfaceGridList object (required for level-specific parameters)
  awghtObjectList_ex9 <- list(awghtObject_level1, awghtObject_level2)
  class(awghtObjectList_ex9) <- "multivariateSurfaceGridList"
  
  LKinfo_ex9<- LKrigSetup(sGrid_advanced, NC = sidelen_bf_1, NC.buffer=0,
                          alpha = c(0.3, 0.7), 
                          nlevel = 2, 
                          a.wghtObject =  awghtObjectList_ex9)
  
  simField_ex9 <- LKrig.sim( sGrid_advanced, LKinfo_ex9)
  imagePlot(as.surface(sGrid_advanced, simField_ex9), col = turbo(256),
            main = "Nonstationary Anisotropy (Different for Each Level)")
  }
  
  #####################################################
  ######## Example 10: Non-monotonically decaying covariance
  ######## Achieved with: 
  ######## Nonstationary a.wght and anisotropy patterns
  ######## Different for each level
  ######## and nonstationary alpha parameter with 
  ######## modifications to create teleconnection
  #####################################################
  \dontrun{
  alpha_side_ex10 <- 39
  alphaGrid_ex10 <- list(
    x = seq(-1, 1, length.out = alpha_side_ex10),
    y = seq(-1, 1, length.out = alpha_side_ex10)
  )
  alphaPoints_ex10 <- make.surface.grid(alphaGrid_ex10)  
  
  # Define a taper for alpha
  alphaTaper_ex10 <- pnorm(alphaPoints_ex10[, 1], mean = 0.0, sd = 0.1) 
  # smooth 0->1 horizontally
  alphaWeights_ex10 <- cbind(1 - alphaTaper_ex10, alphaTaper_ex10)
  
  # Forcing some alpha weights to change for the teleconnection effect
  alphaWeights_ex10[1517:1521, ] <- alphaWeights_ex10[1517:1521, ] + 100
  alphaWeights_ex10[1478:1482, ] <- alphaWeights_ex10[1478:1482, ] + 100
  alphaWeights_ex10[1439:1443, ] <- alphaWeights_ex10[1439:1443, ] + 100
  alphaWeights_ex10[1400:1404, ] <- alphaWeights_ex10[1400:1404, ] + 100
  alphaWeights_ex10[1361:1365, ] <- alphaWeights_ex10[1361:1365, ] + 100
  
  alphaWeights_ex10 <- alphaWeights_ex10 / rowSums(alphaWeights_ex10)  
  # normalize rows to sum to 1
  
  par(mfrow = c(1,2))
  image.plot(
    x = alphaGrid_ex10$x,
    y = alphaGrid_ex10$y,
    z = matrix(alphaWeights_ex10[, 1], 
      nrow = alpha_side_ex10,
      ncol = alpha_side_ex10),
    col = viridis(256),
    main = "Level-1 Alpha Weight"
  )
  image.plot(
    x = alphaGrid_ex10$x,
    y = alphaGrid_ex10$y,
    z = matrix(alphaWeights_ex10[, 2], 
      nrow = alpha_side_ex10, 
      ncol = alpha_side_ex10),
    col = viridis(256),
    main = "Level-2 Alpha Weight"
  )
  par(mfrow = c(1,1))
  
  # Pack as a list of two "surfaceGrid" objects so predict.surfaceGrid() is used
  alphaObject_ex10 <- list()
  for (k in 1:2) {
    alphaSurface_ex10 <- as.surface(alphaPoints_ex10, alphaWeights_ex10[, k]) 
    # has $x, $y, $z (matrix)
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
  
  }
  
  #####################################################
  ######## Covariance plotting
  ######## For Examples 8, 9, and 10
  #####################################################
  
  \dontrun{
  # Plot covariance for Example 8
  plot_cov_surface(
    LKinfo = LKinfo_ex8,
    real_field = simField_ex8,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = TRUE,
    cov_lwd = 1.3,
    nrows_advanced,
    sGrid_advanced
  )
  
  # Plot covariance for Example 9
  plot_cov_surface(
    LKinfo = LKinfo_ex9,
    real_field = simField_ex9,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = TRUE,
    cov_lwd = 1.0,
    nrows_advanced,
    sGrid_advanced
  )
  
  # Plot covariance for Example 10
  plot_cov_surface(
    LKinfo = LKinfo_ex10,
    real_field = simField_ex10,
    loc_x = 0,
    loc_y = 0,
    plotting = TRUE,
    cov_contour = FALSE,
    nrows_advanced= nrows_advanced,
    sGrid_advanced= sGrid_advanced
  )
  }

}