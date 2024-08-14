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
LKrig.lnPlike<- function(GCholesky, Q, quad.form, nObs, nReps,
                         weights, LKinfo) {
  # sanity check on lambda being ratio of tau^2 and sigma2
  # if tau and sigma2 are passed.
  tau <- LKinfo$tau
  sigma2   <- LKinfo$sigma2
  lambda<- LKinfo$lambda
  if (!is.na(sigma2)) {
    if ((tau^2/sigma2) != lambda) {
      stop(" tau, sigma2 and lambda do not match up")
    }
  }
  #
  n <- nObs
  m <- dim(Q)[1]
  # find log determinant of reg matrix for use in the log likeihood
  lnDetReg <- 2 * sum(log(diag(GCholesky)))
  # log determinant of precision matrix.
  lnDetQ <- 2 * sum(log(diag(chol(Q, memory = LKinfo$choleskyMemory))))
  # now apply a miraculous determinant identity (Sylvester's theorem)
  #  det( I + UV) = det( I + VU)    providing UV is square
  # or more generally
  #  det( A + UV) = det(A) det( I + V A^{-1}U)
  #  or as we use it
  #  ln(det( I + V A^{-1}U)) = ln( det(  A + UV)) - ln( det(A))
  #
  lnDetCov <- lnDetReg - lnDetQ + (n - m) * log(lambda) - sum(log(weights))
  # MLE estimate of sigma2 and tau
  # these are derived by assuming Y is  MN(  Ud, sigma2*M )
  sigma2.MLE <- quad.form/n
  tau.MLE <- sqrt(lambda * sigma2.MLE)
  # the  log profile likehood with  sigma2.MLE  and  dhat substituted
  # leaving a profile for just lambda.
  # note that this is _not_  -2*loglike just the log and
  # includes the constants
  lnProfileLike <- (-(n/2) - log(2 * pi) * (n/2) - (n/2) * log(sigma2.MLE) - 
                      (1/2) * lnDetCov)
  # find log likelihood without profiling if tau and sigma2 have been passed.
  # NOTE: this assumes that  lambda == tau^2/sigma2
  if (!is.na(sigma2)) {
    lnLike <- (-(quad.form)/(2 * sigma2) - log(2 * pi) * (n/2) - (n/2) * 
                 log(sigma2) - (1/2) * lnDetCov)
    lnLike.FULL <- sum(lnLike)
  } else {
    lnLike <- NA
    lnLike.FULL <- NA
  }
  # the full ln profile likelihood  is found 
  # by assuming the replicate fields (columns of y)
  # are independent and putting in a pooled MLE for sigma2.
  sigma2.MLE.FULL <- mean(sigma2.MLE)
  tau.MLE.FULL <- sqrt(lambda * sigma2.MLE.FULL)
  lnProfileLike.FULL <- nReps*(-(n/2) - log(2 * pi) * (n/2) - (n/2) * log(sigma2.MLE.FULL) - 
                                 (1/2) * lnDetCov )
  
  #
  return(list(lnProfileLike = lnProfileLike,
              sigma2.MLE = sigma2.MLE,
              tau.MLE = tau.MLE,
              tau = tau,
              sigma2 = sigma2,
              lnLike = lnLike,
              lnLike.FULL = lnLike.FULL, 
              sigma2.MLE.FULL = sigma2.MLE.FULL,
              tau.MLE.FULL = tau.MLE.FULL,
              lnProfileLike.FULL = lnProfileLike.FULL, 
              quad.form = quad.form,
              lnDetCov = lnDetCov))
}


