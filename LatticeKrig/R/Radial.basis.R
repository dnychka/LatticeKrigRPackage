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

Radial.basis <- function(x1, centers, basis.delta,
                     max.points = NULL, 
                  mean.neighbor = 50,
                  BasisFunction = "WendlandFunction", 
                  distance.type = "Euclidean",
                        verbose = FALSE,
                         timing = FALSE)
{    
    d <- ncol(x1)
    n1 <- nrow(x1)
    if (is.null(max.points)) {
        Nmax <- n1 * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
    if( verbose){
        cat('Radial.basis Nmax',         Nmax       , fill=TRUE)
    	cat('Radial.basis basis.delta' , basis.delta, fill=TRUE)
    }
    # evalute RBF basis functions at the x1 locations with
    # nodes given by centers. Returned is a sparse matrix in
    # a simpler format than spam  (row/column indices and value)
    # returned values ind and rd below in 'out' object.
    out <- LKrigDistance( x1, centers,
                             delta = basis.delta,
                        max.points = max.points,
                     mean.neighbor = mean.neighbor,
                     distance.type = distance.type,
                        components = FALSE)
    # evaluate distance  with RBF ---  usually Wendland2.2 
    out$ra <- do.call(BasisFunction, list(d = out$ra/basis.delta) )
    out <- spind2spam(out)
    
    #
    # following code takes much longer
    # t3<- system.time(
    #   out <- spam(out[c("ind", "ra")],
    #               nrow=out$da[1],
    #               ncol= out$da[2] )
    # )
    
    return(out)
}

