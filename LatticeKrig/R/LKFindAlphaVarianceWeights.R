LKFindAlphaVarianceWeights <- function(x1, LKinfo, level) {
  # always weight basis functions by scalar alpha
  if (is.null(LKinfo$alphaObject[[level]])) {
    wght <- LKinfo$alpha[[level]]
  }
  else{
    #  compute spatially varying alpha extension if alphaObject
    # is included.
    wght <- LKinfo$alpha[[level]] * c(predict(LKinfo$alphaObject[[level]], x1) )
  }
  return(wght)
}