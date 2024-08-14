 LKFindSigma2VarianceWeights<- function (x1,LKinfo){
   if (!is.null( LKinfo$sigma2.object) ) {
     wght <- c(predict(LKinfo$sigma2.object, x1))
   }
   else{
     wght<- ifelse( is.na( LKinfo$sigma2), 1.0, LKinfo$sigma2 )
   }
   return(wght)
 }