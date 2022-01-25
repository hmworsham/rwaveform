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
#' @return The deconvolved return waveform
#' @import data.table
#' @import rPeaks
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
decom.waveforms <- function(rawarray, deconarray, smooth=T, thres=0.22, window=3) {
  
  # Define return and geo arrays
  #re = wfarray$return
  re = deconarray
  geol = wfarray$geolocation
  
  # Run adaptive decomposition algorithm on clipped returns
  # Use error handling to identify erroneous or un-decomposable returns
  safe_decomp = function(x){
    tryCatch(decom.adaptive(x, smooth = smooth, thres = thres, width = window), 
             error = function(e){NA})}
  
  # Apply the safe decomposition to the set
  decom = apply(re, 1, safe_decomp)
  
  # Filter out returns that threw exceptions
  successes = which(!is.na(decom)) 
  decom = decom[successes]
  geol = geol[successes]
  
  # Pull Gaussian parameters
  rfit = do.call('rbind', lapply(decom, '[[', 1)) # Indices of correctly processed waveforms
  gpars = do.call('rbind', lapply(decom, '[[', 3)) # The Gaussian parameters
  
  # Get indices of waveforms that need to be reprocessed
  problem_wfs = setdiff(as.numeric(re[,1]$index), rfit[!is.na(rfit)])
  problem_index = which(lapply(decom, 'is.null')==T)
  
  # Preserve parameters that resulted from successful decomposition
  repars = gpars[!is.na(gpars[,1]),]
  colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
  geol = geol[!problem_index]
  geolcols <- c(1:9,16)
  colnames(geol)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
  
  return(list('repars' = repars, 'geolocation' = geol))
}
