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

LKDefaultVarNames<- function(A, tag=NULL){
  if( is.null(A)){
    return(NULL)
  }
  if( ! is.matrix(A)){
    stop("Not a matrix! ")
  }
  colA<- colnames(A)
  if( is.null(colA)){
    if( is.null(tag)){ 
      stop("need a tag for column names")
    }
    colA<- paste0(tag,1:ncol(A))
  }
  return(colA)
}
