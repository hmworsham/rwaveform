#' Apply Waveform Decomposition to N Waveforms
#' 
#' This function applies an iterative Gaussian decomposition to a set of waveforms of dim (MxN).
#' 
#' @param rawarray is a list of the full set of waveform output from the ingest function. Includes return, outgoing, system impulse response and outgoing impulse response arrays, if all available.
#' @param deconarray is the set of deconvolved waveform returns
#' @param method is the deconvolution function using either Gold or Richardson-Lucy algorithm. method=c("Gold","RL"). Default is method=c("Gold"). 
#' @param np is a threshold value which specifies when small_paras (fewer iterations and repetitions) or large_paras (more iterations and repetitions) should be applied in the deconvolution, given the number of estimated peaks in the waveform. More complex waveforms require more iterations and more repetitions to converge. Default is 2 peaks.
#' @param rescale determines whether to rescale the waveform intensity. Returns are rescaled to the minimum waveform intensity to conduct rescaling. Default is TRUE.
#'
#' @return The decomposed return waveform
#' @import data.table
#' @import rPeaks
#' @import pbmcapply
#' @examples
#' decon <- deconv.apply(wfarrays, cliparrays, method = 'Gold', np=2, rescale = F)
#' 
#' # Check dims of result
#' dim(decon.dt)
#' 
#' # Count number of weird returns
#' length(decon.dt[rowSums(decon.dt!=0)])
#' 
#' # Plot some deconvolved returns to check
#' plot(seq(500), decon[1:500,16], type = 'l')
#' lines(seq(500), return1[16][[1]], type = 'l')

#' @keywords waveform, deconvolution, Gold, Richardson-Lucy

reprocess_ <- function(probwfs, re, ge, ge2, rf, gp, th, wd){
  
  # Set problematic waveforms to list
  decon.prob = lapply(as.list(as.data.frame(t(re[re$index %in% probwfs]))), as.numeric)
  geol.prob = ge[ge$index %in% probwfs,]
  
  # Use error handling to identify erroneous or un-decomposable returns
  safedecomp_ = function(x, peakfix=peakfix, smooth = smooth, thres, width){
    tryCatch(decom.adaptive(x, peakfix=peakfix, smooth = smooth, thres, width), 
             error = function(e){NA})}
  
  # Attempt decomposition with stricter threshold parameters
  decomp.prob = pbmclapply(
    decon.prob,
    safedecomp_,
    peakfix=peakfix,
    smooth=smooth,
    thres=th,
    width=wd,
    mc.cores=getOption("mc.cores", ceiling(detectCores()-2))
  )
  
  # Filter out returns that threw exceptions
  successes2 = which(!is.na(decomp.prob) & lengths(decomp.prob)>=3)
  idx3 = probwfs[successes2]
  
  if (length(successes2)) {
    decomp3 = decomp.prob[successes2]
    geol3 = geol.prob[geol.prob$index %in% idx3,]
    
    # Pull indices and Gaussian parameters
    rfit.prob = do.call('rbind', lapply(decomp3, '[[', 1)) # Indices of correctly processed waveforms
    gpars.prob = do.call('rbind', lapply(decomp3, '[[', 3)) # The Gaussian parameters
    
    # Bind re-processed problematic waveform indices and Gaussian params to first successful set
    rfit2 = rbind(rf, rfit.prob)
    gpars2 = rbind(gp, gpars.prob)
    geol2 = rbind(ge2, geol3)
    print(c('Lengths of rfit, gpars & geol2:', nrow(rfit2), nrow(gpars2), nrow(geol2)))
  
  } else {
      rfit2 = rf
      gpars2 = gp
      geol2 = ge2
    }
  
  return(list('rfit'=rfit2, 'gpars' = gpars2, 'geol2'=geol2, 'successidx'=idx3))
}

#' @export 
# Function to run waveform decomposition
decom.apply <- function(rawarray, deconvolved=T, deconarray, peakfix=F, smooth=T, thres=0.22, window=3) {
  
  if (deconvolved){
    re <- deconarray
  } else {
    re <- rawarray$re
  }
  
  # Define return and geo arrays
  geol = rawarray$geol
  
  ## Convert the arrays to lists for batch deconvolution
  idx = re$index
  decon2 = lapply(as.list(as.data.frame(t(re))), as.numeric)
  
  # Run adaptive decomposition algorithm on clipped returns
  # Use error handling to identify erroneous or un-decomposable returns
  safedecomp_ = function(x, peakfix=peakfix, smooth = smooth, thres = thres, width = window){
    tryCatch(decom.adaptive(x, peakfix=peakfix, smooth = smooth, thres = thres, width = window), 
             error = function(e){NA})}
  
  # Apply safe decomposition to the set
  decomp = pbmclapply(
    decon2,
    safedecomp_,
    peakfix=peakfix,
    smooth=smooth,
    thres=thres,
    width=window,
    mc.cores=getOption("mc.cores", ceiling(detectCores()-2))
  )
  
  # Filter out returns that threw exceptions
  successes = which(!is.na(decomp) & lengths(decomp)>=3)
  decomp2 = decomp[successes]
  idx2 = idx[successes]
  geol2 = geol[geol$index %in% idx2,]
  
  # Pull Gaussian parameters
  rfit = do.call('rbind', lapply(decomp2, '[[', 1)) # Indices of correctly processed waveforms
  gpars = do.call('rbind', lapply(decomp2, '[[', 3)) # The Gaussian parameters
  
  # Get indices of waveforms that need to be reprocessed
  problem_wfs = setdiff(as.numeric(idx), as.numeric(idx2))
  
  # Retry with more strict params
  if (length(problem_wfs)) {
    print(paste('There were', length(problem_wfs), 'problems. Retrying with stricter parameters.'))
    retry1 = reprocess_(
      probwfs=problem_wfs, 
      re=re, 
      ge=geol, 
      ge2=geol2, 
      rf=rfit, 
      gp=gpars, 
      th=0.25, 
      wd=3)
    problem_wfs = setdiff(as.numeric(problem_wfs), as.numeric(retry1$successidx))
    rfit = retry1$rfit
    gpars = retry1$gpars
    geol2 = retry1$geol2
    
    # Retry with extremely strict params
    if(length(problem_wfs)) {
      print(paste('There were', length(problem_wfs), 'problems. Retrying with strictest parameters.'))
      retry2 = reprocess_(
        problem_wfs, 
        re, 
        geol, 
        geol2, 
        rfit, 
        gpars, 
        0.99, 
        9)
      problem_wfs = setdiff(as.numeric(problem_wfs), as.numeric(retry2$successidx))
      rfit = retry2$rfit
      gpars = retry2$gpars
      geol2 = retry2$geol2
    }
  }
  
  # Preserve parameters that resulted from successful decomposition
  repars = gpars[!is.na(gpars[,1]),]
  colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
  geolcols <- c(1:9,16)
  colnames(geol2)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
  
  return(list('repars' = repars, 'geolocation' = geol2))
}