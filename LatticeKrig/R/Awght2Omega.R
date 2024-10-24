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



omega2Awght<- function( omega, LKinfo){
  # NOTE: typically in 2D the paremeter used is
  #  kappa where  a.wght=  4 + kappa^2
  #  omega = log(kappa)
  #
  # for a rectangle this should be:
  #  Awght <- 4 +  exp( omega)^2
  # in general  a.wght =  2*dim + exp( dim*omega)
  #
  xDimension<- dim(LKinfo$x)[2]
  Awght<- LKinfo$floorAwght + exp( omega * xDimension )
  return( Awght )
}

Awght2Omega<- function( Awght, LKinfo){
  # for a rectangle this should be:
  #  Awght <- 4 +  exp( omega)^2
  # omega  <- log(  Awght - 4)/2
  if( any(Awght <= LKinfo$floorAwght) ){
    stop(paste("Awght is less than or equal to (lower limit) " ,
               LKinfo$floorAwght) )
  }
  xDimension<- dim(LKinfo$x)[2]
  omega<-  log(Awght - LKinfo$floorAwght)/xDimension
  return( omega )
}