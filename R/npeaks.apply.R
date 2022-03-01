#' Apply Peak Counter to N Waveforms
#' 
#' This function allows you to quickly count the number of waveform peaks in an array of N waveforms
#' @param wfarray is the array of waveform intensities.
#' @param drop is a tuple indicating first and last elements of a range of values we should ignore (e.g., index or non-intensity info). Default is c(0,0), which means use the all input data.
#' @param smooth is tell whether you want to smooth the waveform to reduce the effect of some obvious noise. Default is TRUE.
#' @param threshold is to determine if the detected peak is the real peak whose intensity should be higher than threshold*maximum intensity. Default is 0.2.
#' 
#' @return A list of the number of peaks in each waveform vector
#' @import data.table
#' @import parallel
#' @examples
#' np <- deconv.apply(wfarrays, cliparrays, method = 'Gold', np=2, rescale = F)
#' 

#' @keywords waveform, deconvolution, peaks

#' @export 
#' 
npeaks.apply <- function(wfarray, drop=c(0,0), smooth=F, thres=0.2) {
  
  wfl = lapply(as.list(as.data.frame(t(wfarray))), as.numeric)
  
  np = mclapply(
    wfl,
    rwaveform::npeaks,
    drop=drop,
    smooth=smooth,
    thres=thres,
    mc.cores=getOption("mc.cores", ceiling(detectCores()*.25)))
    
  return(np)
}