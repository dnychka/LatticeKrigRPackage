
# simulation example for estimating covariance
# LatticeKrig
# Copyright 2004-2016, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

suppressMessages(library(LatticeKrig))
#options(echo = FALSE)
test.for.zero.flag <- 1

# NOTE: these tests also exercise the LKinfoUpdate
# function which is complex and 
# has many potential places for creative bugs!

# most of these test are predicated on the replications (M) overwhelming sampling error to 
# give the expected estiamtes ....  

M<-1
N<- 40 # number of obs
set.seed(222)
x<- matrix( runif(2*N ), N,2)
lambdaTrue<- .1^2
#NOTE: true sigma2 is 1.0 dont add fixed function so likelihood is precise.
LKinfo1<- LKrigSetup(x,NC=4, nlevel=3, a.wght= 4.2, nu=1,
                         NC.buffer=0)

look<- LKRectangleCreateLattice(LKinfo1)
cat(" Testing LKRectangleCreateLattice default case 1", fill=TRUE)
for( k in 1:3){
test.for.zero(
  (LKinfo1$latticeInfo$grid[[k]])$x, (look$grid[[k]])$x)
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$y, (look$grid[[k]])$y)
}


LKinfo1<- LKrigSetup(cbind( c(0,1),
                            c(0,1.2)),
                    NC=4, nlevel=3, a.wght= 4.2, nu=1,
                     NC.buffer=0)

look<- LKRectangleCreateLattice(LKinfo1)
delta<- look$delta
cat("  Testing LKRectangleCreateLattice default case 2", fill=TRUE)
for( k in 1:3){
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$x, (look$grid[[k]])$x)
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$y, (look$grid[[k]])$y)
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$x, seq( 0,1, delta[k]))
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$y, seq( 0,1.2, delta[k]))
}
  


# test with buffer

LKinfo1<- LKrigSetup(cbind( c(0,1),
                            c(0,1.2)),
                     NC=4, nlevel=3, a.wght= 4.2, nu=1,
                     NC.buffer=5)
LKinfo2<- LKinfo1
LKinfo2$delta<- LKinfo1$latticeInfo$delta
look<- LKRectangleCreateLattice(LKinfo2)

cat("testing LKRectangleCreateLattice with NC.buffer = 5 ", fill=TRUE)
delta<- LKinfo2$delta
for( k in 1:3){
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$x, (look$grid[[k]])$x)
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$y, (look$grid[[k]])$y)
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$x,
    seq( 0-delta[k]*5, 1+delta[k]*5, delta[k]) )
  test.for.zero(
    (LKinfo1$latticeInfo$grid[[k]])$y,
    seq( 0-delta[k]*5, 1.2+delta[k]*5, delta[k])
  )
}


look<- LKRectangleCreateLattice(LKinfo2)

cat("testing LKRectangleCreateLattice with arbitrary delta NC.buffer = 5 ", fill=TRUE)

delta<-  c( .1,.05,.025)

LKinfo3<- LKrigSetup(cbind( c(0,1),
                            c(0,1.2)),
                     nlevel=3, a.wght= 4.2, nu=1,
                     NC=4,
                     NC.buffer=5, 
                     delta= delta)

LKinfo3$delta<- delta
look<- LKRectangleCreateLattice(LKinfo3)

for( k in 1:3){
  test.for.zero(
    (look$grid[[k]])$x,
    seq( 0-delta[k]*5, 1+delta[k]*5, delta[k]) )
  test.for.zero(
    (look$grid[[k]])$y,
    seq( 0-delta[k]*5, 1.2+delta[k]*5, delta[k])
  )
}





