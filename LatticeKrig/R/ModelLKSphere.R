
# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

###################################################################
## LKrig model for data on Sphere using icosohedral grid.
##  Zach Thomas and Doug Nychka authors
####################################################################

###### Geometry class of this method is LKSphere

## These are obvious defaults for this model and saves
## specifying them in LKrigSetup
# This function is called in LKrigSetup before
# creating the lattice

setDefaultsLKinfo.LKSphere <- function(object, ...) {
  object$floorAwght<- 1.0 
# Definitely do not want Euclidean by default!
# set to Great Circle 
  if( object$distance.type == "Euclidean"){
      object$distance.type <-"GreatCircle"
      attr( object$distance.type, "Radius" )<- 1.0
  }
  
# If the radius attribute is missing set to  1.0
  if( is.null( attr( object$distance.type, "Radius" )) ){
      attr( object$distance.type, "Radius" )<- 1.0
  }
  
# A lazy default: set alpha to 1.0 if only one level.
  if (object$nlevel == 1 & is.na(object$alpha[1])) {
    object$alpha <- list(1.0)
  }
#  
# hard wire the fixed part to just fit a constant function (m =1)
  if( !is.null( object$fixedFunction)){
    object$fixedFunction <- "LKrigDefaultFixedFunction"	
    object$fixedFunctionArgs$m <- 1
  }
# lazy default: set a.wght close to 6 giving a thin plate spline-like 
# model 
# (For the 12 points with 5 neighbors an adjustment is made
  if (is.na(object$a.wght)) {
    object$a.wght <- 1.1
  }
  # delta cutoffs for support of the basis functions found empirically ...
  # Nearest neighobors are within delta great circle distance (and second order
  # neighbors are excluded) 
 
  if( is.null( object$delta)){
          object$delta <- 1.408
  }
  
  #print( object$delta)
  
  if(length( object$delta) > 1 ){
    stop("delta must be just one value")
  }
          
  
  return(object)
  }

# setup the lattice based on subdividing the faces of
# an icosohedron. There are 12 points at first
# level 
# Note that spatial coordiante passed in are assumed as 
# lon/lat in degrees with lon being [-180,180] to 
# be consistent with R maps package
 LKrigSetupLattice.LKSphere<- function(object, x=NULL, verbose,                                         
                                      ... ){
  if( is.null(x)){
    x<- object$x
  }
  startingLevel<- object$setupArgs$startingLevel
  ###### some common setup opertations to all geometries
  LKinfo<- object
  if(  class( LKinfo)[1] != "LKinfo") {
    stop("object needs to an LKinfo object")
  }
  
  rangeLocations<- apply( x,2, "range")
  nlevel<- LKinfo$nlevel
###### end common operations  
# if x not passed then  assume full Sphere
  if( is.null(x)){
    x<- cbind( c(-90,90), c( -180,180))
  }
  if (is.null(startingLevel)) {
    stop("Need to specify startingLevel initial geodesic grid level")
  }
# Here startingLevel is used as setting the initial or coarsest resolution geodesic grid. We may want to change the
# Right now, only allow use of
# up to 7th resolution grid due to memory issues
# if startingLevel=3 and nlevel =5 then centers at resolutions 3,4, and 5 will
# be generated
  if (startingLevel+nlevel-1 > 8){
    stop("startingLevel+nlevel-1 cannot exceed 8")
  }else{
    ##This vector allows us to subset the full list of grids and just keep the ones we need.
    R<- seq(startingLevel,startingLevel+nlevel-1,1) 
    Rmax<- max(R)
  }
# delta cutoffs found empirically ...
# Nearest neighobors are within delta great circle distance (and second order
# neighbors are excluded) 
  #       maybe allow this to be passed as 
  #       object$setupArgs$deltaBaseValue
  # deltaBaseValue  <- object$delta
  deltaBaseValue  <- 1.408
  delta<- deltaBaseValue* 2^( -(0:(Rmax-1)) )
  delta.save<- delta[R]
  
##Build and subset geodesic grid up to level startingLevel+nlevel-1; 
## returns each level in a list but in 
# 3-d coordinates (a.k.a. direction cosines)
  MultiGrid<- IcosahedronGrid(Rmax) ##Get full list of geodesic grids stopping at level Rmax
  grid.all.levels<- list()
  grid3d.all.levels<- list()
  mLevel<- rep(NA,nlevel)
  for(l in (1:nlevel) ){
    # to Sphere converts to lon/lat coordinates
    grid3d<- MultiGrid[[ l + (startingLevel -1) ]]
    gridTemp<- toSphere( grid3d )
      # trim to range of locations (in lon/lat coordinates)
    ind<-  gridTemp[,1] >= rangeLocations[1,1] &
           gridTemp[,1] <= rangeLocations[2,1] &
           gridTemp[,2] >= rangeLocations[1,2] &
           gridTemp[,2] <= rangeLocations[2,2] 
    ind2<-  (1:length( ind)) <= 12
# first 12 coordinates are always the initial isocosohedron points    
    grid.all.levels[[l]]<-  gridTemp[ ind, ] 
    grid3d.all.levels[[l]]<-    grid3d[ ind, ]
    numberNodes<- sum( ind)
    if( numberNodes <= 3){
      stop("must have at least 4 nodes at lowest 
            resolution in the  spatial domain.
            Try increasing startLevels or increasing the 
           boundaries of the spatial domain.")
       }
# A possible logical attribute to indicate which are the
# initial vertices icosohedral
# within the subset determined by ranges     
#    attr(grid.all.levels[[l]],"pentagon")<- ind2[ind]
    mLevel[l]<- nrow( grid.all.levels[[l]] )             
  }
  m<- sum(mLevel)
  offset <- as.integer(c(0, cumsum(mLevel)))
  out<- list( m=m,
              offset=offset,
              mLevel=mLevel, 
              delta=delta.save, 
              rangeLocations=rangeLocations,
# specific arguments for LKSphere              
              startingLevel=startingLevel,
              grid=grid.all.levels,
            grid3d= grid3d.all.levels)
  return(out)
}


LKrigSAR.LKSphere = function(object, Level, ...) {
  if( Level>7){
    stop("can not handle more than 7 levels")}
# delta is the cutoff to find nearest neighbors -- tuned to this
# particular lattice.
  grid <- object$latticeInfo$grid[[Level]]
  grid3d <- object$latticeInfo$grid3d[[Level]]
  delta <- object$latticeInfo$delta[[Level]]
  a.wght <- object$a.wght[[Level]]
## Find nearest 5 or 6 neighbors
  dType<- object$distance.type
  n <-  nrow(grid)
#  logical for pentagon points
# sparse distance matrix in spind format.
# $ind are the indices that are nonzero
# Since original 12 iocosahedral vertices are stored as the first 12 points in each grid, we can easily
#  set their diagonal entries differently if this is needed..
# The B matrix is already in spind format and the entries  (ra) are 
# converted to be the SAR matrix.
  B = LKDist(  grid[,1:2], grid[,1:2], delta = delta,
                  distance.type= dType)
# find the diagonal elements of the sparse matrix
# compute weights for slightly unequal distributions
  ind1<- B$ind[,1]
  ind2<- B$ind[,2]
  raTemp<- rep( NA, length( ind1))
  Diagonal<- ind1==ind2
  if( sum( Diagonal)!= nrow(grid)) {
    stop( "Number of diagonal elements in B different from grid")}
# fill diagonal elements  
    raTemp[ Diagonal ] <- a.wght
    for (I in 1:n ){
# find positions in matrix entries that correspond to the Ith node point
# and is not the diagonal one.
      indNeighbors <- (ind1 == I) & !Diagonal 
# indices of nearest neighbors    
      J<- ind2[ indNeighbors ]
      nJ<- length(J) 
      # this should be either 5 or 6  nearest neighbors for the full sphere
      # but may  be less if the region is a lon/lat rectangle bounding the data.
      x1<- grid3d[J,]
      x0<- grid3d[I,]
      u<- projectionSphere( x0,x1) 
    # u are local 2 d coordinates on tangent plane to sphere at x0
    # x0 projects to (0,0)
      X<- cbind( rep( 1,nJ), u )
    # find weights to predict a linear function
    # at node from the nearest neighbors. 
      if( nJ>=3){
      c2<- c( (X)%*%(solve( t(X)%*%X, c( 1,0,0) )  ))
    #   NOTE: c2 sum to 1 by properties of unbiasedness
    #   for constant function.
      }
      else{
        c2<- 1/nJ
      }
      raTemp[ indNeighbors ]<- -1*c2
    }
#    
  B$ra<- raTemp
  return(B)
# NOTE: B is converted to spam format in LKrig.precision   
}



# return the lattice centers at a given Level
LKrigLatticeCenters.LKSphere<- function(object, Level, ... ){
  return( object$latticeInfo$grid[[Level]] )
} 

