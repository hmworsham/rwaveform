#' Deconvolve, Decompose, and Geolocate Waveforms for a specified geographic area
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

#'@export
process.wf.clip <- function(fp, plt, datadir, shapedir, buff, logpath, outdir){
  
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
  
  # Log the number of returns in the array
  ParallelLogger::logTrace(
    paste(
      filenam, 
      ':',
      length(re),
      'returns in array.'))
  
  # Exit if no returns
  if(length(re)==0) return(NULL)
  
  # Function to process single waveform vector
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
      ParallelLogger::logTrace(
        paste(
          filenam,
          ': Deconvolved wf cleaned.'))
    }
    
    # Return to list
    decon = as.numeric(t(decon))
    
    # Fix problematic peaks
    try.pf = tryCatch({
      decon = peakfix(decon)
    },
    error=function(e) {
      ParallelLogger::logTrace(
        paste(
          filenam,
          ':',
          idx,
          ': Error in cleaning. Vector removed.'
        )
      )
      return(NULL)
    })
    
    # If result of try.pf wrapper is NA, return null result and exit
    if(is.null(try.pf) | is.null(decon)) return(list(idx, NULL, NULL))
    
    # Smooth with window size 3
    decon <- runmean(decon, 3, "C")
    
    # Get peak metrics to start out
    try.pm = tryCatch({
      peak.metrics = peak.info(decon, idx, 0.01)
    },
    error=function(e) {
      ParallelLogger::logTrace(
        paste(
          filenam,
          ':',
          idx,
          ': Failed to derive peak metrics.'
        )
      )
      return(NULL)
      })

    # If result of try.pm wrapper is NA, return null result and exit
    if(is.null(try.pm)) return(list(idx, NULL, NULL))

    # Restore index
    decon = c(idx, decon)

    # Find n peaks
    np = nrow(peak.metrics)

    # If no peaks, return null result and exit
    if(is.null(np) | np==0){
      ParallelLogger::logTrace(
        paste(
          filenam,
          ':',
          idx,
          ': Deconvolution yielded 0 peaks.'))
      return(list(idx, NULL, NULL))
    }

    # ParallelLogger::logTrace(
    #   paste(
    #     filenam,
    #     ':',
    #     idx,
    #     'Deconvolution steps complete'))

    # Define safe decomp function
    safedecomp = function(x, thres){

      # Use error handling to identify erroneous or un-decomposable returns
      tryCatch({
        decom(x, peakmetrics=peak.metrics, peakfix=F, smooth=F, thres=0.1)
      },
      error = function(e) {
        return(list(idx, NULL, NULL))
      })
    }
  
    # If unreasonable peaks, take largest and resmooth
    if(np>16){
      peak.metrics = peak.metrics[order(peak.metrics$A, decreasing=T),][1:16,]
      decon = runmean(decon[-1], 7, "C")
      decon = c(idx, decon)
      
      # Decompose with a conservative threshold for safety
      decom = safedecomp(decon, thres=0.1)
      
    } else {
      
      # Decompose with very low threshold
      th = 0.001
      decom = safedecomp(decon, thres=th)
  
      # If it fails, retry with incrementally stricter thresholds
      # thset = seq(0.1, 0.5, length.out=4)
      # ldc = length(decom)
      # for(th in thset){
      if(is.null(decom[[3]])){
        # ParallelLogger::logTrace(
        #   paste(filenam, ':', idx, ldc,': Decomp failed. Trying with t =', i))
        th = 0.1
        decom = safedecomp(decon, thres=th)
        ldc = length(decom)
      # }
      } # if
    } # else

    # Finally, if the result is null, substitute deconvolution estimates
    if(is.null(decom[[3]])) {

      pm <- matrix(NA, nrow(peak.metrics), 9)

      pm[,1] <- peak.metrics$A
      pm[,4] <- NA

      pm[,2] <- peak.metrics$mu
      pm[,5] <- NA

      pm[,3] <- peak.metrics$sigma
      pm[,6] <- NA

      pm[,7] <- peak.metrics$w
      pm[,8] <- peak.metrics$fs
      pm[,9] <- peak.metrics$ha

      pmi = cbind(idx,pm)
      colnames(pmi) = c('index', 'A', 'u', 'sigma', 'A_std','u_std', 'sigma_std', 'pw', 'fs', 'ha')
      ga = NA

      decom = list(idx, ga, pmi)

      ParallelLogger::logTrace(
        paste(filenam, ':', idx, ': Decomposition failed. Estimated params from deconvolved peaks.'))
    }

    if(is.null(decom[[3]])) {
      return(list(idx, NULL, NULL))
      } else {
        gpar = decom[[3]]
        return(list(idx, gpar, geol.v))
      }
  }
  
  # Run `decon.decom` processing function in parallel over all vectors
  results = mcmapply(
    processvector,
    out,
    re, 
    geol,
    outir,
    sysir,
    mc.preschedule = T, 
    mc.cores = getOption("mc.cores", ceiling(detectCores()-2))
  )
  
  # writedecon <- function(results, filenam, outdir){
  #   deconout <- do.call('rbind', lapply(lapply(results),'[', 4), '[[', 1) # The deconvolved return
  #   outname = paste0(plt, '_', filenam, '_deconvolved.csv')
  #   write.csv(deconout, file.path(outdir, outname))
  # }
  # 
  # tryCatch({
  #   writedecon(results, filenam, outdir)
  # }, error = function(e) {
  #   return(ParallelLogger::logError(paste(filenam, ': failed.')))    
  # })
  
  # Function to assemble decon.decom results, geotransform, and write out to csv
  assemble <- function(results, filenam, outdir){
    
    # Pull Gaussian parameters
    #rfit = do.call('rbind', lapply(lapply(results, '[', 1), '[[', 1)) # Indices of correctly processed waveforms
    #gpars = do.call('rbind', lapply(lapply(results, '[', 2), '[[', 1)) # The Gaussian parameters
    #geols = do.call('rbind', lapply(lapply(results, '[', 3), '[[', 1)) # The geolocations
    indices = do.call('rbind', results[seq(1, length(results), 3)])
    gpars = do.call('rbind', results[seq(2, length(results), 3)])
    geols = do.call('rbind', results[seq(3, length(results), 3)])
    
    # Filter NAs and coerce parameters to data.table
    repars = gpars[!is.null(gpars[,1]),]
    repars = data.table(repars)
    colnames(repars) = c('index', 'A', 'u', 'sigma', 'A_std','u_std', 'sigma_std', 'pw', 'fs', 'ha')
    
    # Filter NAs and coerce geolocation data to data.table
    geolcols <- c(1:9,16)
    geols = geols[!is.null(geols[1,]),]
    geols = data.table(geols)
    colnames(geols)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
    
    # Geotransform waveforms to points
    wfpts = geotransform(decomp = repars, geols)
    
    # Write geotransformed results to csv
    outname = paste0(plt, '_', filenam, '_returnpoints.csv')
    write.csv(wfpts, file.path(outdir, outname), row.names = F)
    
    return(ParallelLogger::logTrace(paste(filenam, ': succeeded.')))
  }
  
  # Run `assemble` with error handling
  tryCatch({
    assemble(results, filenam, outdir)
  }, error = function(e) {
    return(paste(ParallelLogger::logError(paste(filenam, ': failed.', e))))    
  })
  
  # Cleanup
  clearLoggers()
  gc()
  Sys.sleep(0.5)
  
}

