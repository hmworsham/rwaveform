#' Clean Error-Generating Waveforms for Decomposition
#' 
#' This function finds returns in a set of raw or deconvolved waveforms that are likely to generate errors in decomposition and removes them from the set. 
#' 
#' @param x is a list of waveform arrays, such as that resulting from use of the `rwaveform::ingest` function
#'
#' @return a list of cleaned/trimmed waveform arrays
#' @import data.table
#' @examples
#' 
#' @keywords waveform, LiDAR, Gaussian decomposition

#' @export 

# Function to find waveforms that will break Gaussian fit procedure
gfitesc <- function(x){
  xre <- x$return
  apply(xre, 1, function(xwf){
    tre <- try(decom.adaptive(
        xwf, 
        smooth=T,
        peakfix=T,
        thres=0.2,
        width=3),
      silent = T)
    if(class(tre) == 'try-error'){
      tre <- NULL
    }
    return(tre)
  })
}

#' @export
# Function to find indices of problematic waveforms
gfitindxs <- function(breaks){
  okindx = which(lengths(breaks)>0)
  return(okindx)
}

#' @export
# Function to remove those waveforms
rmbreaks <- function(indxs, wf){
  wfre = wf$return
  wfout = wf$outgoing
  wfgeo = wf$geolocation
  newre = wfre[indxs]
  newout = wfout[indxs]
  newgeo = wfgeo[indxs]

  return(list('return' = newre, 'outgoing' = newout, 'geolocation' = newgeo))
}