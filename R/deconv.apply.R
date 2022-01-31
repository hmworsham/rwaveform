#' Apply Waveform Deconvolution to N Waveforms
#' 
#' This function applies a choice of Gold or Richardson-Lucy deconvolution algorithm to a set of waveforms of dim(M x N).
#' @param rawarray is a list of the full set of waveform output from the ingest function. Includes return, outgoing, system impulse response and outgoing impulse response arrays, if all available.
#' @param subarray is either the full set of waveforms or a subset, the output of the clipwf function
#' @param method is the deconvolution function using either Gold or Richardson-Lucy algorithm. method=c("Gold","RL"). Default is method=c("Gold"). 
#' @param np is a threshold value which specifies when small_paras (fewer iterations and repetitions) or large_paras (more iterations and repetitions) should be applied in the deconvolution, given the number of estimated peaks in the waveform. More complex waveforms require more iterations and more repetitions to converge. Default is 2 peaks.
#' @param rescale determines whether to rescale the waveform intensity. Returns are rescaled to the minimum waveform intensity to conduct rescaling. Default is TRUE.
#'
#' @return The deconvolved return waveform
#' @import data.table
#' @import rPeaks
#' @import parallel
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
#' 
#' 
deconv.apply <- function(rawarray, subarray, method = 'Gold', np = 2, rescale = T, small_paras=list(c(30,2,1.8,30,1.2,2)), large_paras=list(c(40,4,1.8,40,4,1.8))){
  
  # Install rPeaks if not already installed
  if (!requireNamespace("rPeaks", quietly = TRUE)) {
    devtools::install_github("jrminter/rPeaks")
  }
  
  # Repeat the system and outgoing impulse responses to nrow of return pulse
  outir_rep = rawarray$outir[rep(seq_len(nrow(rawarray$outir)), nrow(subarray$re))]
  sysir_rep = rawarray$sysir[rep(seq_len(nrow(rawarray$sysir)), nrow(subarray$re))]
  
  # Remove the index columns -- we'll replace them later
  out = subset(subarray$out, select = -index)
  re = subset(subarray$re, select = -index)
  outir_rep = subset(outir_rep, select = -index)
  sysir_rep = subset(sysir_rep, select = -index)
  
  # Convert the arrays to lists for batch deconvolution
  re1 <- lapply(as.list(as.data.frame(t(re))), as.numeric)
  out1 <- lapply(as.list(as.data.frame(t(out))), as.numeric)
  sysir2 <- lapply(as.list(as.data.frame(t(sysir_rep))), as.numeric)
  outir2 <- lapply(as.list(as.data.frame(t(outir_rep))), as.numeric)
  
  # Run deconvolution
  decon = pbmcmapply(rwaveform::deconvolution,
                   re1,
                   out = out1, 
                   imp = sysir2,
                   imp_out = outir2,
                   method = method,
                   np = np,
                   rescale = rescale,
                   small_paras = small_paras,
                   large_paras = large_paras,
                   mc.cores=getOption("mc.cores", ceiling(detectCores()/2)))
  
  # Transpose deconvolution result
  tdecon <- t(decon)
  
  # Create a data table of dim nrow tdecon, ncol decon 
  decon.dt <- data.table(index=1:nrow(tdecon), tdecon)
  
  return(decon.dt)
}