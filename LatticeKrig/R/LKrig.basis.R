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

LKrig.basis <- function(x1, LKinfo, Level= NULL, 
                        raw = FALSE, verbose = FALSE)
  {
    nlevel        <- LKinfo$nlevel
#    delta         <- LKinfo$latticeInfo$delta
#    overlap       <- LKinfo$basisInfo$overlap

    # don't normalize if raw is TRUE
    normalize     <- LKinfo$normalize & !raw

    normalizeMethod <- LKinfo$normalizeMethod
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
    
    # hook to just evaluate subset  of  levels
    # if level argument in not NA
    # if level is NA then evaluate all levels. 
    if(is.null(Level) ){
      levelRange<-  1:nlevel
    }
    else{
      levelRange<- Level
    }
    for (l in levelRange) {
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
         #  cat("normalization method: ", normalizeMethod, fill=TRUE )
        	t2<- system.time(
# depending on option chosen, the basis functions will either be normalized by the default exact method, 
# the fft interpolation alone, the kronecker alone, or both combined 
# both will select fft for coarser levels, kronecker for levels where the number number of basis functions
# violates the minimum coarse grid size condition
        	  wght <- switch(  
        	    normalizeMethod,  
        	    "exact"= LKrigNormalizeBasis( LKinfo,  Level=l,  PHI=PHItemp),  
        	    "exactKronecker"= LKrigNormalizeBasisFast(LKinfo,  Level=l,  x=x1),  
        	    "fftInterpolation"= LKrigNormalizeBasisFFTInterpolate(LKinfo, Level=l, x1=x1),
        	    "both" = LKrigNormalizeBasisSelector(LKinfo, Level = l, x1 = x1, verbose = verbose)
        	      ) 
            	)

            	if( verbose){
            		cat("time for normalization", "Method=", normalizeMethod,  fill=TRUE)
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
           
        }   
       
# At this point if normalization has been done
# the variance of the process at level i will be 1.0 at the x1 locations
# the next two modifications to the basis happen regardless of
# the normalization. 
# alpha are the weights applied at each level 
# (makes sense that these sum to one but they do not need to )
# sigma2 is an overall weight applied across all levels
# If sum(alpha) = 1 and the basis is normalized
# then the process will have marginal variance sigma2 and each level will 
# have variance alpha[level]*sigma2. This product can vary over space and is 
# motivates  including the objects defining surface of these parameters. 
        
        
wght<- LKFindAlphaVarianceWeights(x1,LKinfo, l)

#  spam does not handle diag of one element so do separately  
        if( length( wght)>1){
          t3<- system.time(
            PHItemp <- diag.spam(sqrt(wght)) %*% PHItemp
            )
          }
        else{
            PHItemp <-sqrt(wght)*PHItemp
        }
# accumulate this level into growing basis matrix        
        PHI <- spam::cbind.spam(PHI, PHItemp)
    }
#   
# finally multiply the basis functions by sqrt(sigma2) to give the right
# marginal variances. This is either a function of the locations or constant.
# sigma2 may not be provided when estimating a model 
     wght<- LKFindSigma2VarianceWeights(x1,LKinfo)
     
    if( length( wght)>1){
      PHI <- diag.spam(sqrt(wght)) %*% PHI
      }
    else{
      PHI <-sqrt(wght)*PHI
      }
#############    attr(PHI, which = "LKinfo") <- LKinfo
    return(PHI)
}
