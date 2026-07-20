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

LKrigMakewU <- function(object, verbose = FALSE) {
	LKinfo<- object$LKinfo
	if (!is.null(object$U)) {
		wU <- sqrt(object$weights) * object$U
	} else {
		if (!is.null(LKinfo$fixedFunction)) {
			wU <- sqrt(object$weights) * do.call(
			LKinfo$fixedFunction, 
			 c(list(x = object$x, 
			     	Z = object$Z,
	    distance.type = LKinfo$distance.type),
	        LKinfo$fixedFunctionArgs))
		}
		else{
			wU<- NULL
			}
	}
	if (verbose) {
    cat("dim wU:", dim(wU),  fill=TRUE)
	}
	return( wU)
}
