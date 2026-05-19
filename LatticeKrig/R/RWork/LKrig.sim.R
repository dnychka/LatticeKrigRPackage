##BEGIN HEADER
#
# LatticeKrig is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
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
# A copy of the GNU General Public License is included
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or refer to  http://www.r-project.org/Licenses/GPL-2
#
##END HEADER

LKrig.sim <- function(x1, LKinfo, M = 1, just.coefficients = FALSE,
                      timing=FALSE) {
	
  t1<- system.time(
    Q <- LKrig.precision(LKinfo)
  )
	#
	#  Q is precision matrix of the coefficients -- not of the field
#  last step does the multiplication to go from coefficients to evaluating
#  values at the field
#  Q = t(H)%*%H = inv((Sigma)
#  So   Sigma= Hi%*% t(Hi)
#  find u= t(Hi) %*% N(0,1)   then cov(u) = t(Hi)%*%Hi
#  Hi is upper triangular
#
# snippet of code to test the algebra ...
#   x<-seq( 0,1,,20); Sigma<- exp(-rdist( x,x)/2.5); Q<- solve( Sigma)
#   Mc <- chol(Q); H<- Mc ; Hi<- solve(H);
#   test.for.zero( Q, t(H)%*%H); test.for.zero(Sigma, Hi%*%t(Hi))
#   E<- rnorm(20);  u1<- Hi%*% E ;   u2<-backsolve(Mc,E)
#   test.for.zero(u1,u2)
#
  t2<- system.time(
    Qc <- chol(Q, memory = LKinfo$choleskyMemory)
  )
	m <- LKinfo$latticeInfo$m
	E <- matrix(rnorm(M * m), nrow = m, ncol = M)
 t3<- system.time(
   randomC <- backsolve(Qc, E)
 )
	if (just.coefficients) {
		return(randomC)
	} 
	else {
	  t4<- system.time(
		PHI1 <- LKrig.basis(x1, LKinfo)
	  )
	  t5<- system.time(
	    simFields<- PHI1 %*% randomC
	  )
	  if( timing){
	    print( rbind(Q= t1, Chol= t2, backsolve= t3, PHI=t4, PhiC= t5))
	  }
	  
		return(simFields)
	}
}	

