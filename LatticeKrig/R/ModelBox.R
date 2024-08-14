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
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

## LKrig model for 3-d data in a box
#
setDefaultsLKinfo.LKBox <- function(object, ...) {
	# object == LKinfo
  object$floorAwght<- 6
  
	if (is.null(object$NC.buffer)) {
		object$NC.buffer <- 2
	}
	#a lazy default: Set alpha to 1 if only one level.
	if (object$nlevel == 1 & is.na(object$alpha[1])) {
		object$alpha <- list(1)
	}
	if (is.null(object$setupArgs$a.wght)) {
		object$setupArgs$a.wght <- 6.01
	}
	return(object)
}

LKrigSAR.LKBox <- function(object, Level, ...) {
	m <- object$latticeInfo$mLevel[Level]
	a.wght <- (object$a.wght)[[Level]]
	if (length(a.wght) > 1) {
		stop("a.wght must be constant")
	}
	da <- c(m, m)
	# INTIALLY create all arrays for indices ignoring boudaries
	#  e.g. an edge really only has 2 or 3 neighbors not 4.
ra <- c(rep(a.wght, m), rep(-1, m * 6))
	Bi <- c(rep(1:m, 7))
	Bindex <- array(1:m, object$latticeInfo$mx[Level, ])
	# indexing is East, West, South, North, Down, Up.
	Bj <- c(1:m, c(LKArrayShift(Bindex, c(-1, 0, 0))), c(LKArrayShift(Bindex, 
		c(1, 0, 0))), c(LKArrayShift(Bindex, c(0, -1, 0))), c(LKArrayShift(Bindex, 
		c(0, 1, 0))), c(LKArrayShift(Bindex, c(0, 0, -1))), c(LKArrayShift(Bindex, 
		c(0, 0, 1))))
	inRange <- !is.na(Bj)
	Bi <- Bi[inRange]
	Bj <- Bj[inRange]
	ra <- ra[inRange]
	return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}

LKrigLatticeCenters.LKBox <- function(object, Level, ...) {
	gridl <- object$latticeInfo$grid[[Level]]
	# return the grid describing the centers -- not the center themselves
	class(gridl) <- "gridList"
	return(gridl)
}


LKrigSetupLattice.LKBox <- function(object, verbose=FALSE,  
	...) {
  ######
  if (ncol(object$x) != 3) {
    stop("x is not 3-d !")
  }
  LKinfo <- object
  ##### some checks
  if (class(LKinfo)[1] != "LKinfo") {
    stop("object needs to be an LKinfo object")
  }
  if (is.null(LKinfo$NC)& is.null(LKinfo$delta)) {
    stop("Need to specify NC  or delta for grid size")
  }
  if (is.null(LKinfo$NC.buffer)) {
    stop("Need to specify NC.buffer for lattice buffer region")
  }
  if(!is.null(LKinfo$delta )){
    if( length(LKinfo$delta)!= LKinfo$nlevel)
    {
      stop("length of delta is different from number of levels")
    }
  } 
  # create lattice
  out<- LKBoxCreateLattice(LKinfo, verbose=verbose)
  
  return(out)
}

LKBoxCreateLattice<- function(LKinfo, verbose=FALSE){
  # these two items used throughout by passed or created 
  # lattice info
  
  nlevel<- LKinfo$nlevel
  NC.buffer<- LKinfo$NC.buffer
  
  # find range of scaled locations
  if (is.null(LKinfo$basisInfo$V)) {
    Vinv <- diag(1, 3)
  } else {
    Vinv <- solve(LKinfo$basisInfo$V)
  }
  
  if( is.null(LKinfo$rangeDomain)){
    rangeDomain <- apply(LKinfo$x %*% t(Vinv), 2, "range")
  }
  else{
    rangeDomain<- LKinfo$rangeDomain
  }
  
  if (verbose) {
    cat("Domain range", rangeDomain, fill = TRUE)
  }
  
  grid.info <- list( range = rangeDomain)
  # set the coarsest spacing of centers           
  
  # create basis function spacings if delta is not passed
  if( is.null(LKinfo$delta)){
    NC<- LKinfo$NC
    dMax<- max( c(rangeDomain[2,]- rangeDomain[1,] ))
    delta.level1 <- dMax/(NC - 1)
    deltaVector <- rep(NA, nlevel)
    for (j in 1:nlevel) {
      deltaVector[j] <- delta.level1/(2^(j - 1))
    }
  }
  else{
    deltaVector<- LKinfo$delta
  }
  #
  #  if NC.buffer is greater than zero the number of grid points in
  # both dimensions will be increased by this buffer around the edges:
  # a total of  NC + 2* NC.buffer grid points along the longer dimension
  #
  mx <- mxDomain <- matrix(NA, nrow = nlevel, ncol = 3)
  #
  # build up series of nlevel multi-resolution grids
  # and accumulate these in a master list -- creatively called grid
  grid.all.levels <- NULL
  # loop through multiresolution levels decreasing delta by factor of 2
  # and compute number of grid points.
  # build in a hook for buffer regions to differ in x and y
  # currently this distinction is not supported.
  NC.buffer.x1 <- NC.buffer
  NC.buffer.x2 <- NC.buffer
  NC.buffer.x3 <- NC.buffer
  
  for (j in 1:nlevel) {
    delta <- deltaVector[j]
    # the width in the spatial coordinates for NC.buffer grid points at this level.
    buffer.width.x1 <- NC.buffer.x1 * delta
    buffer.width.x2 <- NC.buffer.x2 * delta
    buffer.width.x3 <- NC.buffer.x3 * delta
    # rectangular case
    grid.list <- list(
      x1 = seq(rangeDomain[1,1] - buffer.width.x1,
               rangeDomain[2,1] + buffer.width.x1, delta),
      x2 = seq(rangeDomain[1,2] - buffer.width.x2,
               rangeDomain[2,2] + buffer.width.x2, delta),
      x3 = seq(rangeDomain[1,3] - buffer.width.x3,
               rangeDomain[2,3] + buffer.width.x3, delta)
    )
    class(grid.list) <- "gridList"
    mx[j, 1] <- length(grid.list$x1)
    mx[j, 2] <- length(grid.list$x2)
    mx[j, 3] <- length(grid.list$x3)
    mxDomain[j, ] <- mx[j, ] - 2 * NC.buffer            
    grid.all.levels <- c(grid.all.levels, list(grid.list))
  }
  # end multiresolution level loop
  #	
  # create a useful index that indicates where each level starts in a
  # stacked vector of the basis function coefficients.
  mLevel <- mx[, 1] * mx[, 2] * mx[, 3]
  offset <- as.integer(c(0, cumsum(mLevel)))
  m <- sum(mLevel)
  
  mLevelDomain <- (mx[, 1] - 2 * NC.buffer.x1) * 
    (mx[, 2] - 2 * NC.buffer.x2) *
    (mx[, 3] - 2 * NC.buffer.x3)
  
  out <- list(     m = m, 
                   offset = offset, 
                   mLevel = mLevel, 
                   delta = deltaVector,                 
                   rangeDomain = rangeDomain
  )
  # specific arguments for LKBox 
  out <- c(out, list(mLevelDomain = mLevelDomain,
                     mx = mx, 
                     mxDomain = mxDomain, 
                     NC.buffer = NC.buffer,
                     grid = grid.all.levels, 
                     grid.info = grid.info)
  )
  
  return(out)
  
}








