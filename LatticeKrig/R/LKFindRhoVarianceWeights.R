 LKFindRhoVarianceWeights<- function (x1,LKinfo){
   if (!is.null( LKinfo$rho.object) ) {
     wght <- c(predict(LKinfo$rho.object, x1))
   }
   else{
     wght<- ifelse( is.na( LKinfo$rho), 1.0, LKinfo$rho )
   }
   return(wght)
 }