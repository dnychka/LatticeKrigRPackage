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

LatticeKrigEasyDefaults<- function( argList,nlevel,x){
	  if( is.null(argList$LKGeometry) ){
                xDimension<- ncol( x)
                argList$LKGeometry<- c( "LKInterval", "LKRectangle", "LKBox")[xDimension]
                NC<- argList$NC      
                if( is.null(NC)){
                  N<- nrow(x)
                  a<-  sum(  2^(xDimension*(0:(nlevel-1)))  )
                  NCtest<-  (N/a)^( 1/xDimension)            

                  argList$NC<- round(max(4, NCtest ))
              }
#                  
# This is a convenient function to hard wire good  default choices for 
# some of the more common models and where the 
# parameters that are not explicitly given -- basically to save typing
#
                
                # NC chosen so that with  d= xDimensison   NCtest^d * ( 1 + 2^d  + 2^(2d) + ...)
                #  gives a basis size that is at least the number of observations.
                # Note that if NC.buffer is not zero this can still add quite few extra basis function 
                # outside the domain. 

                if( is.null(argList$a.wght)){
              	  # a thinplate spline-like a.wght
                  argList$a.wght<- 2*xDimension +.01
                }
                if( is.null(argList$nu)){
                	argList$nu<-1
                	}
	  }
  #
  # make sure overlap is 2.0 for cubic spline tensor basis
  # but don't mess with this  related arguments if they are passed 
              EuclideanModel<- is.element( argList$LKGeometry,
                                   c("LKInterval", "LKRectangle", "LKBox" )) 
              if(EuclideanModel ){
              #
              if( !is.null(argList$BasisFunction) ) {
                if( argList$BasisFunction == "CubicBSpline" &
                    is.null(argList$overlap ))
              argList$overlap<- 2.0
              #
              if( is.null(argList$BasisType) ){
                argList$BasisType<-"tensor"
              }
              #
              }
              }
              
              #
              return( argList)
}
