# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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

LKrig.basis <- function(x1, LKinfo, verbose = FALSE)
  {
    nlevel        <- LKinfo$nlevel
#    delta         <- LKinfo$latticeInfo$delta
#    overlap       <- LKinfo$basisInfo$overlap
    normalize     <- LKinfo$normalize
    distance.type <- LKinfo$distance.type
    fast          <-  attr( LKinfo$a.wght,"fastNormalize")
    V <- LKinfo$basisInfo$V
# coerce x1 to a matrix    
    x1<- as.matrix( x1)
    # We will transform (scale and rotate) x matrix of locations by   x%*% t(solve(V))
    # 
    # For the radial basis functions this
    # will imply that the Euclidean distance norm used between two vectors X1 and X2
    # is  t(X1-X2) solve(A) (X1-X2)  with A = (V %*%t(V))
    # Here's the deal on the linear transformation V:
    # It should be equivalent to just transforming the locations outside of LatticeKrig.
    # Where ever the observation locations are used
    # they are first transformed with t(solve(V))).
    # Surprisingly this only needs to happen in one place below and in LKrig.setup to
    # determine the grid centers.
    #
    # The RBF centers and the delta scale factor, however, assumed to be
    # in the transformed scale and so a  transformation is not needed.
    # see LKrig.setup for how the centers are determined.
    # The concept is that if LKrig was called with the transformed locations
    # ( x%*%t( solve(V)) instead of x
    # and V set to diag( c(1,1) ) one would get the same results.
    # as passing a non identity V.
    #
    # accumulate matrix column by column in PHI
    PHI <- NULL
    # transform locations if necessary (lattice centers already in 
    # transformed scale) 
    if( !is.null( V[1]) ){     
      x1<- x1 %*% t(solve(V))     
    }
    if( verbose){
          cat("LKrig.basis: Dim x1 ",  dim( x1), fill=TRUE)
    }
    basis.delta <- LKrigLatticeScales(LKinfo)
    for (l in 1:nlevel) {
        # Loop over levels and evaluate basis functions in that level.
        # Note that all the center information based on the regualr grids is
        # taken from the LKinfo object
        #  set the range of basis functions, they are assumed to be zero outside
        #  the radius basis.delta and according to the distance type.
        # 
        # There are two choices for the type of basis functions
        # 
        centers<- LKrigLatticeCenters( LKinfo,Level=l )
        if(LKinfo$basisInfo$BasisType=="Radial" ){ 
        	 t1<- system.time(
        PHItemp <- Radial.basis(  x1, centers, basis.delta[l],
                                max.points = LKinfo$basisInfo$max.points,
                             mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                             BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                             distance.type = LKinfo$distance.type,
                                   verbose = verbose)
                             )
                             }
        if(LKinfo$basisInfo$BasisType=="Tensor" ){  
        	 t1<- system.time(            
        PHItemp <- Tensor.basis(  x1, centers, basis.delta[l],
                                max.points = LKinfo$basisInfo$max.points,
                             mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                             BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                             distance.type = LKinfo$distance.type)
                             )
                             }      	                            
        if( verbose){
          cat(" Dim PHI level", l, dim( PHItemp), fill=TRUE)
          cat("time for basis", fill=TRUE) 
          print( t1)
        }
        if (normalize) { 
        	t2<- system.time( 
        	if( !fast ){
        	# the default choice should work for all models	
               wght<-  LKrigNormalizeBasis( LKinfo,  Level=l,  PHI=PHItemp)               
            }
            else{ 
            # potentially a faster method that is likely very specific
            # to a particular more e.g. LKrectangle with stationary first order GMF
            # NOTE that locations are passed here rather than the basis matrix	
               wght<- LKrigNormalizeBasisFast(LKinfo,  Level=l,  x=x1)
            	}
            	)
                
            	if( verbose){
            		cat("time for normalization", "fast=", fast,  fill=TRUE)
            		print( t2)
            	}
# now normalize the basis functions by the weights treat the case with one point separately
# wghts are in scale of inverse marginal variance of process
# the wght maybe zero if the x location does not overlap with nay basis function (e.g. x outside spatial
# adjust for this case by setting wght =1. This should be valid because the row of PHI will be zero
                indZero <- wght == 0
                if( any( indZero) ){
                    warning( "Some normalization weights are zero")
                }
                wght[ indZero] <- 1.0
           if( nrow( x1)>1){
              PHItemp <- diag.spam( 1/sqrt(wght) ) %*% PHItemp
           }
          else{
              PHItemp@entries <- PHItemp@entries/sqrt(wght)
              }
         if( verbose){
          cat("finding normalized basis functions", fill=TRUE)
        }    
        }   
       
# NOTE: when spam is loaded this is equivalent to just cbind( PHI, PHItemp)
        
# always weight basis functions by scalar alpha
        if (is.null( LKinfo$alphaObject[[l]]) ){
         wght <- LKinfo$alpha[[l]]
        }
        else{
# spatially varying alpha extension         
          wght <- LKinfo$alpha[[l]]*c(predict(LKinfo$alphaObject[[l]], x1))
          }
#  spam does not handle diag of one element so do separately  
        if( length( wght)>1){
          t3<- system.time(
            PHItemp <- diag.spam(sqrt(wght)) %*% PHItemp
            )
          if( verbose){
            cat("time for weights", t3,  fill=TRUE)
            print( t3)
          }
          }
        else{
            PHItemp <-sqrt(wght)*PHItemp
        }
# accumulate this level into growing basis matrix        
        PHI <- spam::cbind.spam(PHI, PHItemp)
    }
 #   
   
    
# finally multiply the basis functions by sqrt(rho) to give the right
# marginal variances. This is either a function of the locations or constant.
# rho may not be provided when estimating a model 
    if (!is.null( LKinfo$rho.object) ) {
      wght <- c(predict(LKinfo$rho.object, x1))
    }
    else{
      wght<- ifelse( is.na( LKinfo$rho), 1.0, LKinfo$rho )
    }
      
    if( length( wght)>1){
      PHI <- diag.spam(sqrt(wght)) %*% PHI
      }
    else{
      PHI <-sqrt(wght)*PHI
      }
#############    attr(PHI, which = "LKinfo") <- LKinfo
    return(PHI)
}


