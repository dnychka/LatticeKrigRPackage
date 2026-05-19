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

 convertIndexPeriodic <- function(I, nGrid, nPad = NULL) {
	nDim <- length(nGrid)
	if (is.null(nPad)) {
		nPad <- rep(0, nDim)
	}
	J <- rep(0, length(I))
	cI <- cumprod(c(1, nGrid))
	nGridP <- nGrid - 2 * nPad
	cJ <- cumprod(c(1, nGridP))
	L <- I - 1
	#	indP<- NULL
	#	ind<- NULL
for (k in nDim:1) {
		kI <- floor(L/cI[k])
		coordP <- (kI - nPad[k])%%nGridP[k]
		J <- J + (coordP) * cJ[k]
		L <- L - kI * cI[k]
		# 		indP<- cbind( coordP + 1, indP) 	
		# 		ind<- cbind( kI + 1, ind)
}
	J <- J + 1
	return(J)
}
 
convertIndexArray <- function(I, nGrid) {
	nDim <- length(nGrid)
	cI <- cumprod(c(1, nGrid))
	L <- I - 1
	ind <- NULL
	for (k in nDim:1) {
		kI <- (floor(L/cI[k]))
		L <- L - kI * cI[k]
		ind <- cbind((kI + 1), ind)
	}
	return(ind)
}
