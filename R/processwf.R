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
process_wf <- function(fp){
  
  # Ingest waveforms
  wfarrays = ingest(fp)
  print(paste('flightpath', fp, 'ingested'))
  
  # Subset to a manageable set of waveforms for testing
  out <- wfarrays$out
  re <- wfarrays$re
  geol <- wfarrays$geol
  outir <- wfarrays$outir
  sysir <- wfarrays$sysir
  
  set.seed(123)
  sub <- sample(seq(1:nrow(re)), 100000, replace=F)
  #sub <- seq(1:5000)
  #sub <- seq(1:nrow(re))
  out_sub <- out[sub]
  re_sub <- re[sub]
  geol_sub <- geol[sub]
  sub_arrays = list('out'=out_sub, 're'=re_sub, 'geol'=geol_sub)
  
  # deconvolve
  decon <- rwaveform::deconv.apply(
    wfarrays, 
    sub_arrays, 
    method='Gold', 
    rescale=F,
    small_paras = list(c(30,2,1.2,30,2,2)),
    large_paras = list(c(40,4,1.8,30,2,2)))
  
  # Restore index
  decon$index <- re_sub$index
  
  # Check for NaNs and extreme values
  #decon <- subset(decon, select = -index)
  print(paste('Deconvolved. Result is of dim:', list(dim(decon))))
  nanrows = which(rowSums(is.na(decon))>0)
  bigrows = which(rowSums(decon[,2:length(decon)])>10^5)
  
  # Clean NaNs and extreme values
  if(length(nanrows) | length(bigrows)) {
    decon <- deconv.clean(decon)
    print(paste('Cleaned. Result is of dim:', list(dim(decon))))
  }
  
  # Find npeaks
  decon <- data.table(t(apply(decon, 1, peakfix)))
  np <- apply(decon, 1, npeaks, smooth=F, threshold=0)
  
  # Store indices of returns with potentially unreasonable number of peaks or 0 peaks
  unreasonable <- decon[np>18]$index
  nopeaks <- decon[np==0]$index
  
  # Filter out 0 and unreasonable peak vectors
  if (length(nopeaks) | length(unreasonable)) {
    decon <- decon[!decon$index %in% nopeaks,]
    decon <- decon[!decon$index %in% unreasonable,]
  }
  
  # Decompose waveforms
  decomp <- rwaveform::decom.apply(
    sub_arrays,
    deconvolved=T,
    decon,
    smooth=T,
    peakfix=F,
    thres=0.25,
    window=3)
  
  print(paste('Decomposed. Result is of dim:', list(dim(decomp$repars))))
  
  # geotransform waveforms to points
  wfpts = geotransform(decomp = decomp$repars, decomp$geolocation)
  
  # write geotransformed results to csv
  outname = tail(str_split(fp, '/')[[1]],1)
  write.csv(wfpts, paste0('/global/scratch/users/worsham/geolocated_returns/', outname, '_returnpoints.csv'))
  return(wfpts)
}
