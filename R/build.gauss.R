#' Build Guassian
#'
#' The function allows you to build a sum of Gaussian functions from peak-wise parameters
#'
#' @param x a dataframe of the Gaussian parameters
#' @param tbins an integer indicating the number of time bins (or x-values) to build along

#' @export
build.gauss <- function(x, tbins=500){
  tbins = 1:tbins
  sumyi = 0
  for(i in seq(nrow(x))){
    sumyi = sumyi + x[i,2] * exp(-abs(tbins - x[i,3])**x[i,5]/(2*x[i,4]**2))
  }
  return(sumyi)
}