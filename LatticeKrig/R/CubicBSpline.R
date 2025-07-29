

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


