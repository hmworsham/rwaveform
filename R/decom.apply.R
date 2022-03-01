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
#' @import future
#' @import rPeaks
#' @import doParallel
#' @import parallel
#' @import plyr
#' 
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

# Define function to reprocess
reprocess_ <- function(retries, re, ge, ge2, rf, gp, th, wd){
  
  # Set problematic waveforms to list
  decon.rt = re[re$index %in% retries]
  #decon.rt = lapply(as.list(as.data.frame(t(decon.rt))), as.numeric)
  geol.rt = ge[ge$index %in% retries,]
  
  safedecomp.rt_ = function(x, thres, width){
    
    # Use error handling to identify erroneous or un-decomposable returns
    tryCatch({
      decom.adaptive(x, peakfix=F, smooth=T, thres=thres, width=width)
    },
    error = function(e) {
      return(NA)
    })
  }
  
  # # Attempt decomposition with stricter threshold parameters
  
  # foreach
  # decomp.rt = foreach(i=1:length(decon.rt)) %dopar% {
  #   safedecomp.rt_(decon.rt[i], thres=th, wid=wd)
  # }
  
  # parLapply
  # decomp.rt = parLapply(parallelCluster, decon.rt, safedecomp.rt_, thres=thres, width=window)
  
  # mclapply  
  decomp.rt = apply(
    decon.rt,
    1,
    safedecomp.rt_,
    thres=th,
    width=wd
    #mc.preschedule=F,
    #mc.cores=getOption("mc.cores", ceiling(detectCores()*0.8))
  )
  
  # Filter out returns that threw exceptions
  successes2 = which(!is.na(decomp.rt) & lengths(decomp.rt)>=3)
  idx3 = retries[successes2]
  
  if (length(successes2)) {
    decomp3 = decomp.rt[successes2]
    geol3 = geol.rt[geol.rt$index %in% idx3,]
    
    # Pull indices and Gaussian parameters
    rfit.rt = do.call('rbind', lapply(decomp3, '[[', 1)) # Indices of correctly processed waveforms
    gpars.rt = do.call('rbind', lapply(decomp3, '[[', 3)) # The Gaussian parameters
    
    # Bind re-processed problematic waveform indices and Gaussian params to first successful set
    rfit2 = rbind(rf, rfit.rt)
    gpars2 = rbind(gp, gpars.rt)
    geol2 = rbind(ge2, geol3)
    #print(c('Lengths of rfit, gpars & geol2:', nrow(rfit2), nrow(gpars2), nrow(geol2)))
    
  } else {
    rfit2 = rf
    gpars2 = gp
    geol2 = ge2
  }
  
  return(list('rfit'=rfit2, 'gpars' = gpars2, 'geol2'=geol2, 'successidx'=idx3))
}

#' @export 
decom.apply <- function(rawarray, deconvolved=T, deconarray, peakfix=F, smooth=T, thres=0.2, window=3) {
  
  if (deconvolved){
    re = deconarray
  } else {
    re = rawarray$re
  }
  
  # Define return and geo arrays
  geol = rawarray$geol
  
  # Convert the arrays to lists for batch deconvolution
  idx = re$index
  #decon2 = lapply(as.list(as.data.frame(t(re))), as.numeric)
  
  # Set data.table to single thread for multicore processing
  #setDTthreads(1)
  
  # Set up some parallel backend stuff
  # parallelCluster <- makeCluster(
  #   ceiling(detectCores()*0.8),
  #   type = "FORK",
  #   methods = FALSE)
  # setDefaultCluster(parallelCluster)
  # registerDoParallel(parallelCluster)
  # 
  # # # Set threading to 1 for MKL and data.table
  # clusterEvalQ(cl = parallelCluster, {
  #   library(data.table)
  #   library(devtools)
  #   library(rPeaks)
  #   load_all('~/Repos/rwaveform')
  #   setDTthreads(1)
  # })
  
  # future
  # nworkers <- as.numeric(ceiling(detectCores()*0.8))
  # plan(multicore, workers = nworkers)
  
  # Define safe decomposition function
  safedecomp_ = function(x, thres, width){
    
    # Use error handling to identify erroneous or un-decomposable returns
    tryCatch({
      rwaveform::decom.adaptive(x, peakfix=F, smooth=T, thres, width)
    },
    error = function(e) {
      return(NA)
    })
  }
  
  # Run adaptive decomposition algorithm safely on listed returns
  # parLapply
  # decomp = parLapply(parallelCluster, decon2, safedecomp_, thres=thres, width=window)
  
  # foreach
  # decomp = foreach(i=1:length(decon2)) %dopar% {
  #   decom.adaptive(unlist(decon2[i]), thres=thres, width=window)
  # }
  
  # future
  # decomp = future.apply::future_apply(
  #   re,
  #   1,
  #   safedecomp_,
  #   thres=thres,
  #   width=window
  # )
  
  # llply
  # decomp = llply(
  #   decon2,
  #   safedecomp_,
  #   thres=thres,
  #   width=window,
  #   .parallel=T
  # )
  
  # mclapply
  decomp = apply(
    re,
    1,
    safedecomp_,
    thres=thres,
    width=window
    #mc.preschedule=F,
    #mc.cores=getOption("mc.cores", ceiling(detectCores()*0.8))
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
  
  #Retry, iterating through increasingly strict params
  thres.set = seq(thres+0.1, 0.99, length.out=5)

  for (i in seq(1, length(thres.set))) {
    if (length(problem_wfs)) {
      message(paste('There were', length(problem_wfs), 'problems. Retrying with stricter parameters.'))
      retry1 = reprocess_(
        retries=problem_wfs,
        re=re,
        ge=geol,
        ge2=geol2,
        rf=rfit,
        gp=gpars,
        th=thres.set[i],
        wd=3)
      problem_wfs = setdiff(as.numeric(problem_wfs), as.numeric(retry1$successidx))
      rfit = retry1$rfit
      gpars = retry1$gpars
      geol2 = retry1$geol2
    }
  }
  
  #   # Retry with extremely strict params
  #   if(length(problem_wfs)) {
  #     print(paste('There were', length(problem_wfs), 'problems. Retrying with strictest parameters.'))
  #     retry2 = reprocess_(
  #       retries=problem_wfs,
  #       re=re,
  #       ge=geol,
  #       ge2=geol2,
  #       rfit,
  #       gpars,
  #       0.95,
  #       11)
  #     problem_wfs = setdiff(as.numeric(problem_wfs), as.numeric(retry2$successidx))
  #     rfit = retry2$rfit
  #     gpars = retry2$gpars
  #     geol2 = retry2$geol2
  #   }
  # }
  
  # Preserve parameters that resulted from successful decomposition
  repars = gpars[!is.na(gpars[,1]),]
  colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
  geolcols <- c(1:9,16)
  colnames(geol2)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')
  
  # Stop cluster for doParallel approaches
  # on.exit({
  #   try({
  #     cat("Attempting to stop cluster\n")
  #     stopImplicitCluster()        # package: `doParallel`
  #     stopCluster(parallelCluster) # package: `parallel`
  #   })
  # })
  
  return(list('repars' = repars, 'geolocation' = geol2, 'problems' = problem_wfs))
}