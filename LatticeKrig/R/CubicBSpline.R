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

CubicBSpline <- function(d){
  # cubic spline assuming knots at -1,-1/2,0,1/2, 1 and d is the distance from location to 
  # the central knot.
  # convert to unit knot spacing of -2,-1,0,1,2
  d<- d*2 
  # standardized cubic spline basis function
  # only evaluate distances  <= 2 -- rest are zero. 
  ind<- (abs(d) <= 2)
  out<- rep( 0, length( d))
  # only evaluate distances within 2 units of the the central knot (0).
  dInd<- d[ind]
  # piecewise cubic 
  out[ind]<- 
  (1/6) *(pmax(dInd + 2, 0)^3 - 4 * (pmax(dInd + 1, 0)^3) + 6 * (pmax(dInd, 0)^3) -
                   4 * (pmax(dInd - 1, 0)^3)) 
   return(out)
}


