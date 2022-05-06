# Function to process single waveform vector

decon.decom <- function(out.v, re.v, geol.v, outir.v, sysir.v) {
  
  # Install rPeaks if not already installed
  if (!requireNamespace("rPeaks", quietly = TRUE)) {
    devtools::install_github("jrminter/rPeaks")
  }
  
  # Prep for deconvolution
  idx = as.integer(re.v[1])
  re.v = as.vector(as.numeric(re.v)[-1])
  out.v = as.vector(as.numeric(out.v)[-1])
  geol.v = as.vector(as.numeric(geol.v)) # Keep geol indices
  outir.v = as.vector(as.numeric(outir.v)[-1])
  sysir.v = as.vector(as.numeric(sysir.v)[-1])
  
  # Deconvolve
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
  
  # Clean if NAs or unreasonably large values are present
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
  decon = peakfix(decon)
  
  # Smooth with window size 3
  decon <- runmean(decon, 3, "C")
  
  # Get peak metrics for subsequent fitting
  try.pm = tryCatch({
    peak.metrics = peak.info(decon, idx, 0.01)
  },
  error=function(e) {
    ParallelLogger::logTrace(
      paste(
        filenam,
        ': ', 
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
        ': ',
        idx,
        ': Deconvolution yielded 0 peaks.'))
    return(list(idx, NULL, NULL))
  }
  
  # If unreasonable peaks, take the 16 largest
  if(np>16){
    peak.metrics = peak.metrics[order(peak.metrics$A, decreasing=T),][1:16,]
  }
  
  # ParallelLogger::logTrace(
  #   paste(
  #     filenam, 
  #     ': ',
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
  
  # Decompose with very low threshold
  th = 0.01
  decom = safedecomp(decon, thres=th)
  
  # If it fails, retry with k incrementally stricter thresholds
  thset = seq(0.1, 0.5, length.out=4)
  ldc = length(decom)
  for(i in thset){
    if(is.null(decom[[3]]))
      # ParallelLogger::logTrace(
      #   paste(filenam, ': ', idx, ldc,': Decomp failed. Trying with t =', i))
      decom = safedecomp(decon, thres=i)
    ldc = length(decom)
  }
  
  # If decomposition fails after k attempts, substitute estimates from deconvolution as parameters
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
      paste(filenam, ': ', idx, ': Decomposition failed. Estimated params from deconvolved peaks.'))
  }
  
  if(is.null(decom[[3]])) {
    return(list(idx, NULL, NULL))
  } else {
    gpar = decom[[3]]
    return(list(idx, gpar, geol.v))
  }
}