#' Deconvolve, Decompose, and Geolocate Waveforms
#' 
#' This function combines other functions in the rwaveform package to create a full processing workflow
#' @param fp is string pointing to the filepath of a directory containing one waveform binary package. The package should typically contain, at minimum, header (.hdr) and binary data files describing:
#' - observation parameters 
#' - geolocation
#' - outgoing pulse
#' - return pulse
#' - impulse response
#' - outgoing impulse response
#' @param clip specifies whether to clip the waveform to a particular geometry before processing. Default is not to clip (False)
#' @param aoi is list describing the extent of an area of interest in the form `[xmin, xmax, ymin, ymax]`
#'
#' @return a data.frame of points with return parameters and geolocation information
#' @import data.table
#' @import rPeaks
#' @import parallel
#' @examples
#' 
#' @keywords waveform, LiDAR

#' @export 
#' 
#' 
process_wf <- function(fp, clip=FALSE, aoi){
  
  # Ingest waveforms
  wfarrays = ingest(fp)
  print(paste('flightpath', fp, 'ingested'))
  
  # Clip waveforms to geometry if specified
  if (clip) {
    
    arrays = rwaveform::clipwf(wfarrays, aoi, buff=20)
    
    if(dim(arrays[[1]])[1]==0) {
      print('clip failed: flightpath does not intersect plot')
    } else {
      print('clip succeeded')
    }
  } else {
    arrays = wfarrays
  }
  
  # deconvolve
  decon <- rwaveform::deconv.apply(
    wfarrays, 
    arrays, 
    method='Gold', 
    rescale=F,
    small_paras = list(c(30,2,1.2,30,2,2)),
    large_paras = list(c(40,4,1.8,30,2,2)))
  
  # Check for NaNs and extreme values
  #decon <- subset(decon, select = -index)
  nanrows = which(rowSums(is.na(decon))>0)
  bigrows = which(rowSums(decon[,2:length(decon)])>10^5)
  
  # Clean NaNs and extreme values
  if(length(nanrows) | length(bigrows)) {
    decon <- deconv.clean(decon)
  }
  
  # Find npeaks
  decon <- data.table(t(apply(decon, 1, peakfix)))
  np <- apply(decon, 1, npeaks, smooth=F, threshold=0)
  
  # Store indices of returns with potentially unreasonable number of peaks or 0 peaks
  #unreasonable <- decon[np>12]$index
  nopeaks <- which(np==0)
  
  # Filter out 0 and unreasonable peak vectors
  decon <- decon[-nopeaks,]
  #decon <- decon[-unreasonable,]
  
  # Decompose waveforms
  decomp <- mcmapply(
                  rwaveform::decom.adaptive,
                  decon,
                  smooth=T,
                  peakfix=T,
                  thres=0.2,
                  width=3,
                  mc.cores=getOption("mc.cores", ceiling(detectCores()/2)))
  
  # geotransform waveforms to points
  wfpts = geotransform(decomp = decom$repars, decom$geolocation)
  
  return(wfpts)
}
