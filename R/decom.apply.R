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

#' @export 
# Function to run waveform decomposition
decom.waveforms <- function(rawarray, deconarray, peakfix=F, smooth=T, thres=0.22, window=3) {
  
  # Define return and geo arrays
  #re = wfarray$return
  geol = rawarray$geol
  
  # Remove the index columns -- we'll replace them later
  decon = subset(decon, select = -index)
  #geol = subset(geol, select = -index)
    
  ## Convert the arrays to lists for batch deconvolution
  decon2 = lapply(as.list(as.data.frame(t(decon))), as.numeric)
  #geol2 = lapply(as.list(as.data.frame(t(geol))), as.numeric)
  
  # Run adaptive decomposition algorithm on clipped returns
  # Use error handling to identify erroneous or un-decomposable returns
  safe_decomp = function(x){
    tryCatch(decom.adaptive(x, peakfix=peakfix, smooth = smooth, thres = thres, width = window), 
             error = function(e){NA})}
  
  # Apply safe decomposition to the set
  decomp = pbmcmapply(
    safe_decomp,
    decon2,
    mc.cores=getOption("mc.cores", ceiling(detectCores()/2))
  )
  
  # Filter out returns that threw exceptions
  successes = which(!is.na(decomp))
  decomp = decomp[successes]
  geol = geol[successes]
  
  # Pull Gaussian parameters
  rfit = do.call('rbind', lapply(decomp, '[[', 1)) # Indices of correctly processed waveforms
  gpars = do.call('rbind', lapply(decomp, '[[', 3)) # The Gaussian parameters
  
  # Get indices of waveforms that need to be reprocessed
  problem_wfs = setdiff(as.numeric(re[,1]$index), rfit[!is.na(rfit)])
  problem_index = which(lapply(decomp, 'is.null')==T)
  
  # Preserve parameters that resulted from successful decomposition
  repars = gpars[!is.na(gpars[,1]),]
  colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
  geol = geol[!problem_index]
  geolcols <- c(1:9,16)
  colnames(geol)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
  
  return(list('repars' = repars, 'geolocation' = geol))
}