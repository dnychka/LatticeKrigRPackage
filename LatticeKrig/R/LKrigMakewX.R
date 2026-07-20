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

LKrigMakewX <- function(object, verbose = FALSE) {
	if(verbose){
	cat( "LKrigMakewX: dim X ", dim( object$X), fill=TRUE)	
	}
	if (!is.null(object$X)) {
		wX <- diag.spam(sqrt(object$weights)) %*% object$X
	} 
	else {
	  tmpX<- LKrig.basis(object$x, 
	                     object$LKinfo, verbose = FALSE)
		wX <- diag.spam(sqrt(object$weights)) %*% tmpX
	}
	if(verbose){
	cat( "LKrigMakewX: dim wX ", dim( wX), fill=TRUE)	
	}
	return(wX)
}



