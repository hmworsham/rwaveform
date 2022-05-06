#' Fix waveform peaks
#' 
#' This function finds problematic peaks in a waveform and introduces a minor correction to improve Gaussian curve fitting. After deconvolution, some first peaks in the waveform jump from 0 to the maximum value with no ramp-up. A Gaussian fit can't accommodate this, so when it happens we introduce a value in the prior time bin at 0.99x the first non-zero return value. This creates a shape that a Gaussian can be fit to without invalidating the data.
#' 
#' @param y is a vector of return intensity values
#'
#' @return the peak-fixed vector
#' @import data.table
#' @keywords waveform, deconvolution, peaks
#' @export 

peakfix <- function(y){
  # Find the first nonzero
  firstnonzero = which(y!=0)[1]
  if (is.na(firstnonzero)) {
    #ParallelLogger::logError(paste('No peaks found in vector', y[1]))
    return(NULL)
  } else if (firstnonzero == 1) {
    #ParallelLogger::logError(paste('First peak invalid', y))
    return(NULL)
  } else if (y[[firstnonzero]] >=y [[firstnonzero+1]]){
        y[[firstnonzero-1]] <- 0.99 * y[[firstnonzero]]
        }
  return(y)
  }