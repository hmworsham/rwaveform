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


process_wf <- function(fp, logpath, outdir){
  
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
  
  # Repeat the system and outgoing impulse responses to nrow of return pulse
  outir = outir[rep(seq_len(nrow(outir)), nrow(re))]
  sysir = sysir[rep(seq_len(nrow(sysir)), nrow(re))]
  
  rm(wfarrays)
  
  # Convert the arrays to lists for batch deconvolution
  re = lapply(as.list(as.data.frame(t(re))), as.numeric)
  out = lapply(as.list(as.data.frame(t(out))), as.numeric)
  geol = lapply(as.list(as.data.frame(t(geol))), as.numeric)
  sysir = lapply(as.list(as.data.frame(t(sysir))), as.numeric)
  outir = lapply(as.list(as.data.frame(t(outir))), as.numeric)
  
  ParallelLogger::logTrace(
    paste(
      filenam, 
      ':',
      length(re),
      'returns in array.'))
  
  processvector <- function(out.v, re.v, geol.v, outir.v, sysir.v) {
    
    # Install rPeaks if not already installed
    if (!requireNamespace("rPeaks", quietly = TRUE)) {
      devtools::install_github("jrminter/rPeaks")
    }
    
    # deconvolve
    idx = as.integer(re.v[1])
    re.v = as.vector(as.numeric(re.v)[-1])
    out.v = as.vector(as.numeric(out.v)[-1])
    geol.v = as.vector(as.numeric(geol.v)) # Keep geol indices
    outir.v = as.vector(as.numeric(outir.v)[-1])
    sysir.v = as.vector(as.numeric(sysir.v)[-1])
    
    decon = deconvolution(
      re.v, 
      out.v, 
      sysir.v, 
      outir.v,
      method='Gold',
      np=2,
      rescale=T,
      small_paras = c(20,2,1.2,20,2,1.2),
      large_paras = c(40,4,1.8,30,2,2))
    
    # Clean if NAs or unreasonably large values
    if(sum(is.na(decon))>0 || sum(decon, na.rm=T)>10^5){
      decon = deconv.clean(decon)
    }
    
    # Restore index
    decon = c(idx, decon)
    
    # Fix problematic peaks and find n peaks
    decon = peakfix(decon)
    np = npeaks(decon, drop=c(0,0), smooth=F, thres=0)
    
    # Exit if no peaks or unreasonable peaks
    if(np==0 || np>24){
      return(NA)
    }
    
    decon = as.numeric(t(decon))
    
    # Define safe decomp function
    safedecomp = function(x, thres, width){
      
      # Use error handling to identify erroneous or un-decomposable returns
      tryCatch({
        decom.adaptive(x, peakfix=F, smooth=T, thres=0.2, width=3)
      },
      error = function(e) {
        return(NA)
      })
    }
    
    # Decompose with very low threshold
    th = 0.01
    decom = safedecomp(decon, thres=th, width=3)
    
    # If it fails, retry with incrementally stricter thresholds
    thset = seq(th, 0.99, length.out=5)
    ldc = length(decom)
    for(i in thset[-1]){
      if(ldc<3){
        decom = safedecomp(decon, thres=i, width=3)
        ldc = length(decom)
      } else {
        decom = decom
      }
    }
    
    if(length(decom)<3) return(NA)
    
    gpar = decom[[3]]

    return(list('rfit' = idx, 'gpar'=gpar, 'geol'=geol.v))
  }
  
  results = mcmapply(
    processvector,
    out,
    re, 
    geol,
    outir,
    sysir,
    mc.preschedule = T, 
    mc.cores = getOption("mc.cores", ceiling(detectCores()*0.6))
  )
  
  assemble <- function(results, filenam, outdir){
  
    # Pull Gaussian parameters
    rfit = do.call('rbind', lapply(lapply(results, '[', 1), '[[', 1)) # Indices of correctly processed waveforms
    gpars = do.call('rbind', lapply(lapply(results, '[', 2), '[[', 1)) # The Gaussian parameters
    geols = do.call('rbind', lapply(lapply(results, '[', 3), '[[', 1)) # The geolocations
    
    repars = gpars[!is.na(gpars[,1]),]
    repars = data.table(repars)
    colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
  
    geolcols <- c(1:9,16)
    geols = geols[!is.na(geols[1,]),]
    geols = data.table(geols)
    colnames(geols)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
    
    # geotransform waveforms to points
    wfpts = geotransform(decomp = repars, geols)
    
    # write geotransformed results to csv
    outname = paste0(filenam, '_returnpoints.csv')
    write.csv(wfpts, file.path(outdir, outname))
    
    return(ParallelLogger::logTrace(paste(filenam, ': succeeded.')))
  }
  
  tryCatch({
    assemble(results, filenam, outdir)
    }, error = function(e) {
    return(ParallelLogger::logError(paste(filenam, ': failed.')))    
   })
  
  gc()
  Sys.sleep(5)
}


#'@export
process_wf_clip <- function(fp, plt, datadir, shapedir, buff, logpath, outdir){
  
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
  
  # Clip to plot bounds
  pltext = rwaveform::aoiextent(plt, shapedir)
  wfarraysc = rwaveform::clipwf(wfarrays, pltext, buff=3)
  
  # Separate waveform arrays into distinct variables
  out <- wfarraysc$out
  re <- wfarraysc$re
  geol <- wfarraysc$geol
  outir <- wfarrays$outir
  sysir <- wfarrays$sysir
  
  # Repeat the system and outgoing impulse responses to nrow of return pulse
  outir = outir[rep(seq_len(nrow(outir)), nrow(re))]
  sysir = sysir[rep(seq_len(nrow(sysir)), nrow(re))]
  
  rm(wfarrays)
  
  # Convert the arrays to lists for batch deconvolution
  re = lapply(as.list(as.data.frame(t(re))), as.numeric)
  out = lapply(as.list(as.data.frame(t(out))), as.numeric)
  geol = lapply(as.list(as.data.frame(t(geol))), as.numeric)
  sysir = lapply(as.list(as.data.frame(t(sysir))), as.numeric)
  outir = lapply(as.list(as.data.frame(t(outir))), as.numeric)
  
  ParallelLogger::logTrace(
    paste(
      filenam, 
      ':',
      length(re),
      'returns in array.'))
  
  processvector <- function(out.v, re.v, geol.v, outir.v, sysir.v) {
    
    # Install rPeaks if not already installed
    if (!requireNamespace("rPeaks", quietly = TRUE)) {
      devtools::install_github("jrminter/rPeaks")
    }
    
    # deconvolve
    idx = as.integer(re.v[1])
    re.v = as.vector(as.numeric(re.v)[-1])
    out.v = as.vector(as.numeric(out.v)[-1])
    geol.v = as.vector(as.numeric(geol.v)) # Keep geol indices
    outir.v = as.vector(as.numeric(outir.v)[-1])
    sysir.v = as.vector(as.numeric(sysir.v)[-1])
    
    decon = deconvolution(
      re.v, 
      out.v, 
      sysir.v, 
      outir.v,
      method='Gold',
      np=2,
      rescale=T,
      small_paras = c(20,2,1.2,20,2,1.2),
      large_paras = c(40,4,1.8,30,2,2))
    
    # Clean if NAs or unreasonably large values
    if(sum(is.na(decon))>0 || sum(decon, na.rm=T)>10^5){
      decon = deconv.clean(decon)
    }
    
    # Restore index
    decon = c(idx, decon)
    
    # Fix problematic peaks and find n peaks
    decon = peakfix(decon)
    np = npeaks(decon, drop=c(0,0), smooth=F, thres=0)
    
    # Exit if no peaks or unreasonable peaks
    if(np==0 || np>24){
      return(NA)
    }
    
    decon = as.numeric(t(decon))
    
    # Define safe decomp function
    safedecomp = function(x, thres, width){
      
      # Use error handling to identify erroneous or un-decomposable returns
      tryCatch({
        decom.adaptive(x, peakfix=F, smooth=T, thres=0.2, width=3)
      },
      error = function(e) {
        return(NA)
      })
    }
    
    # Decompose with very low threshold
    th = 0.01
    decom = safedecomp(decon, thres=th, width=3)
    
    # If it fails, retry with incrementally stricter thresholds
    thset = seq(th, 0.99, length.out=5)
    ldc = length(decom)
    for(i in thset[-1]){
      if(ldc<3){
        decom = safedecomp(decon, thres=i, width=3)
        ldc = length(decom)
      } else {
        decom = decom
      }
    }
    
    if(length(decom)<3) return(NA)
    
    gpar = decom[[3]]
    
    return(list('rfit' = idx, 'gpar'=gpar, 'geol'=geol.v))
  }
  
  results = mcmapply(
    processvector,
    out,
    re, 
    geol,
    outir,
    sysir,
    mc.preschedule = T, 
    mc.cores = getOption("mc.cores", ceiling(detectCores()*0.6))
  )
  
  assemble <- function(results, filenam, outdir){
    
    # Pull Gaussian parameters
    rfit = do.call('rbind', lapply(lapply(results, '[', 1), '[[', 1)) # Indices of correctly processed waveforms
    gpars = do.call('rbind', lapply(lapply(results, '[', 2), '[[', 1)) # The Gaussian parameters
    geols = do.call('rbind', lapply(lapply(results, '[', 3), '[[', 1)) # The geolocations
    
    repars = gpars[!is.na(gpars[,1]),]
    repars = data.table(repars)
    colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
    
    geolcols <- c(1:9,16)
    geols = geols[!is.na(geols[1,]),]
    geols = data.table(geols)
    colnames(geols)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
    
    # geotransform waveforms to points
    wfpts = geotransform(decomp = repars, geols)
    
    # write geotransformed results to csv
    outname = paste0(plt, '_', filenam, '_returnpoints.csv')
    write.csv(wfpts, file.path(outdir, outname))
    
    return(ParallelLogger::logTrace(paste(filenam, ': succeeded.')))
  }
  
  tryCatch({
    assemble(results, filenam, outdir)
  }, error = function(e) {
    return(ParallelLogger::logError(paste(filenam, ': failed.')))    
  })
  
  gc()
  Sys.sleep(5)
}

