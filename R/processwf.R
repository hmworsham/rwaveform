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
#' @importFrom ParallelLogger addDefaultFileLogger clearLoggers logInfo logTrace logWarn logError logFatal
#' @examples
#' 
#' @keywords waveform, LiDAR

#' @export
#' 
process_wf_deprec <- function(fp, logpath=logpath, outdir=outdir){
  
  # Define logging filepath
  clearLoggers()
  addDefaultFileLogger(logpath)
  
  # Get flightpath name from path
  filenam = tail(unlist(strsplit(fp, '/')), 1)
  
  # Ingest waveforms
  wfarrays = ingest(fp)
  
  # Check ingest results and log errors 
  if(is.null(wfarrays)) {
    ParallelLogger::logError(
      paste(
        filenam, 
        ': Ingest failed. A data or header file may be missing.'))
    return(NULL)
  } else {
    ParallelLogger::logTrace(
      paste(
        filenam, 
        ': Successfully ingested.'))
  }
  
  # Separate waveform arrays into distinct variables
  out <- wfarrays$out
  re <- wfarrays$re
  geol <- wfarrays$geol
  outir <- wfarrays$outir
  sysir <- wfarrays$sysir
  
  ParallelLogger::logTrace(
      paste(
        filenam, 
        ':',
        nrow(re),
        'returns in array.'))
  
  #set.seed(99)
  #sub <- sample(seq(1:nrow(re)), 10000, replace=F)
  #sub <- seq(1:5000)
  sub <- seq(1:nrow(re))
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
    small_paras = list(c(20,2,1.2,20,2,1.2)),
    large_paras = list(c(40,4,1.8,30,2,2)))
  
  # Check results of deconvolution and log errors
  if(any(unlist(apply(decon, 1, is))=='try-error')) {
    ParallelLogger::logError(
      paste(
        filenam, 
        ': Deconvolution failed for index(es):', 
        which(unlist(apply(decon, 1, is))=='try-error')))
  } else {
    ParallelLogger::logTrace(
      paste(
        filenam, 
        ': All deconvolutions succeeded.'))
  }
  
  # Restore index
  decon$index <- re_sub$index
  
  # Check for NaNs and extreme values
  #print(paste('Deconvolved. Result is of dim:', list(dim(decon))))
  nanrows = which(rowSums(is.na(decon))>0)
  bigrows = which(rowSums(decon[,2:length(decon)])>10^5)

  # Clean NaNs and extreme values
  if(length(nanrows) | length(bigrows)) {
    decon = deconv.clean(decon)
    #print(paste('Cleaned. Result is of dim:', list(dim(decon))))
    ParallelLogger::logWarn(paste(filenam, ':', length(nanrows), 'vectors had NaNs and were cleaned:', paste(nanrows, collapse=',')))
    ParallelLogger::logWarn(paste(filenam, ':', length(bigrows), 'vectors had unrealistically large values and were cleaned:', paste(bigrows, collapse=',')))
  }

  # Find npeaks
  decon = data.table(t(apply(decon, 1, peakfix)))
  np = rwaveform::npeaks.apply(
    decon,
    drop=c(1,1),
    smooth=F,
    thres=0
  )

  # Store indices of returns with potentially unreasonable number of peaks or 0 peaks
  unreasonable = as.integer(which(np>24))
  nopeaks = as.integer(which(np==0))
  
  # Filter out 0 and unreasonable peak vectors
  if (length(nopeaks)) {
    decon = decon[!decon$index %in% nopeaks,]
    ParallelLogger::logWarn(paste(filenam, ':', length(nopeaks), 'vectors had no peaks and were removed:', paste(nopeaks, collapse=',')))
  }
  if (length(unreasonable)) {
    decon = decon[!decon$index %in% unreasonable,]
    ParallelLogger::logWarn(paste(filenam, ':', length(unreasonable), 'vectors had unrealistically many peaks and were removed:', paste(unreasonable, collapse=',')))
  }

  # Decompose waveforms
  decomp = rwaveform::decom.apply(
    sub_arrays,
    deconvolved=T,
    decon,
    peakfix=F,
    smooth=T,
    thres=0.01,
    window=5)
  
  #print(paste('Decomposed. Result is of dim:', list(dim(decomp$repars))))
  if (length(decomp$problems)) {
    ParallelLogger::logWarn(
      paste(
        filenam,
        ':',
        length(decomp$problems), 'vectors failed in decomposition and were removed:', paste(decomp$problems, collapse=',')))
  }
  
  # geotransform waveforms to points
  wfpts = geotransform(decomp = decomp$repars, decomp$geolocation)

  # write geotransformed results to csv
  outname = paste0(filenam, '_returnpoints.csv')
  write.csv(wfpts, file.path(outdir, outname))
  
  return(print('Done'))
}
