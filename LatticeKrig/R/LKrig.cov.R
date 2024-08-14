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

LKrig.cov <- function(x1, x2 = NULL, LKinfo, C = NA, 
                      marginal = FALSE, aRange=NA) {
  # ARange is an unused argument to make this compatable with
  # being called from spatialProcess in fields 
  if( marginal){
    if (!is.null(x2)) {
      stop("x2 should not be passed to find marginal variance")
    }
    #
    #  decide if variance needs to be recomputed from the basis functions
    #  NOTE: if LKinfo$normalize = TRUE  and an exact method used
    #  then basis functions
    #  will be already normalized so that the marginal variance 
    #  is 1.0  the alpha weights are multiplied by the basis 
    #  functions even when not normalized 
    #  
    evaluateVariance <-  !(
      LKinfo$normalize == TRUE &
        (
          LKinfo$normalizeMethod == "exact" |
          LKinfo$normalizeMethod == "exactKronecker"
        )
    )
    
    sigma2Variance<- LKFindSigma2VarianceWeights(x1,LKinfo)
    wght <- rep(0, nrow(x1))
    # loop over level of multiresolution
    for (l in 1:(LKinfo$nlevel)) {
      if (evaluateVariance) {
        # raw =TRUE means an approximate method for normalization
        # is not used. Instead the exact method is applied.
        PHILevel <- LKrig.basis(x1, LKinfo, raw = TRUE, Level = l)
        marginal.variance <-
        LKrigNormalizeBasis(LKinfo,  Level = l,  PHI = PHILevel)
        #  square root of alpha weights and sigma2 weights always multiplied by basis
        # function and so included into LKrigNormalizeBasis
      }
      else{
        marginal.variance <-  
          sigma2Variance*LKFindAlphaVarianceWeights(x1,LKinfo,l)
      }
        wght <- wght +  marginal.variance
    }
    return(wght)
  }
  if( !marginal ){
    # find the full covariance matrix 
    # or the covariance matrix times the matrix C
    # this is just for testing or using with  the fields spatialProcess when the
    # number of observations is  <=  1e3
    # this computation works because the basis functions called below include all the
    # variance parameters. 
    PHI1 <- LKrig.basis(x1, LKinfo)
    
    # sparse precision matrix for the basis coeficients	
    Q <- LKrig.precision(LKinfo)
    
    Qc <- chol(Q, memory = LKinfo$choleskyMemory)
     
    if (is.null(x2)) {
      PHI2 <- PHI1
    } else {
      PHI2 <- LKrig.basis(x2, LKinfo)
    }
    
    if (is.na(C[1])) {
      A <- forwardsolve(Qc, transpose = TRUE, t(PHI2), upper.tri = TRUE)
      A <- backsolve(Qc, A)
    } else {
      A <- forwardsolve(Qc, transpose = TRUE, t(PHI2) %*% C, upper.tri = TRUE)
      A <- backsolve(Qc, A)
    }
    return(PHI1 %*% A)
  }
  
}