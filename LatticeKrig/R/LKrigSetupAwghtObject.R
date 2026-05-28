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

LKrigSetupAwghtObject<- function( object){
  # object is an LKinfo object
  # should already have the lattice node information. 
  nLevel<- object$nlevel
  a.wght<- list()
  for( k in 1: nLevel ){
    latticeGrid<-  make.surface.grid((object$latticeInfo$grid )[[k]])
    # if there is only one object use this for the a.wghts at
    # every level
    if (length(class(object$a.wghtObject)) == 1 &&
        inherits(object$a.wghtObject, "multivariateSurfaceGridList") ){
      aTemp <- predict(object$a.wghtObject, latticeGrid, nLevel, k) 
    }else{
      aTemp <- predict(object$a.wghtObject, latticeGrid)
    }
    # check for missing values -- most likely outside range of spatial domain.    
    if( any( is.na(aTemp)) ){
      print(cbind( latticeGrid, aTemp))
      stop("Some a.wghts set to NAs from predict")
    }
    # NOTE to check this for a 2-d geometry
    # image.plot( as.surface( latticeGrid, aTemp))
    # coerce to a matrix, if aTemp is a vector this will be a one column 
    # matrix
    aTemp<- as.matrix( aTemp )
    a.wght<- c( a.wght, list( aTemp)  )
  }
  return(a.wght)
}