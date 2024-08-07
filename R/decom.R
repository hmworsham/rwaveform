#' decom
#'
#' The function allows you to estimate parameters charcterizing waveforms and to pave the way for generating waveform-based point cloud.
#'
#' @param x is a waveform with a index at the beginning and followed with intensities.
#' @param smooth is tell whether you want to smooth the waveform to remove some obvious outliers. Default is TRUE.
#' @param thres is to determine if the detected peak is the real peak whose intensity should be higher than threshold*maximum intensity. Default is 0.22.
#' @param width width of moving window.Default is 3, must be odd integer between 1 and n.This parameter ONLY work when the smooth is TRUE.

#' @return A list contains estimates of A, u, sig after decomposition.
#' @import caTools
#' @import minpack.lm
#' @importFrom stats na.omit
#' @export
#' @references
#'   Zhou, Tan*, Sorin C. Popescu, Keith Krause, Ryan D. Sheridan, and Eric Putman, 2017. Gold-A novel deconvolution algorithm with
#'   optimization for waveform LiDAR processing. ISPRS Journal of Photogrammetry and Remote Sensing 129 (2017):
#'   131-150. https://doi.org/10.1016/j.isprsjprs.2017.04.021
#'
#' @examples
#'
#' ##import return waveform data
#' library(data.table)
#' data(return)
#' return<-data.table(index=c(1:nrow(return)),return)
#' x<-return[1,] ###must be a dataset including intensity with index at the beginning.
#' r1<-decom(x)
#' r2<-decom(x,smooth=TRUE,width=5) ###you can assign different smooth width for the data
#' ###when it comes very noisy waveforms, it may give you some problems
#' xx<-return[182,]
#' r3<-decom(xx)  ##this one returns NULL which means the function didn't work for the
#'                ##complex waveform or too noisy waveform,we should try to reprocess
#'                ##these unsucessful waveforms using larger width to smooth the waveforms.
#' r4<-decom(xx,smooth=TRUE,width=5) ##when you change to a larger width, it can work,
#'                                   #but give you some unreasonable estimates, return NA
#'
#' ###original result from this decom is (you will not see it, the function filter this result
#' ###and put NA for the estimation since they maybe not right results)
#' #Nonlinear regression model
#' #model:y~A1*exp(-(x-u1)^2/(2*sigma1^2))+A2*exp(-(x-u2)^2/(2*sigma2^2)) n\
#' #+A3*exp(-(x-u3)^2/(2*sigma3^2))
#' #data: df
#' #A1      A2      A3      u1      u2      u3  sigma1  sigma2  sigma3
#' #228.709 -30.883  81.869  41.640  42.131  71.680  14.613   3.522   8.073
#' ##A (ampilitude should not be negative)
#'
#' r5<-decom(xx,width=10) ##this will work by smoothing the waveform
#' r6<-decom(xx,thres=0.1,width=5)  ##by adjusting width and thres of real peak, you may
#'                                  ##get a reasonable results
#' \donttest{
#' # for the whole dataset
#' dr<-apply(return,1,decom)
#' }

#'


decom <- function(x, peakmetrics, smooth=TRUE, peakfix=FALSE, width=3, thres=0.1){
  
  waveform <- as.numeric(x)
  index <- waveform[1]
  y <- waveform[-1]
  y[y==0] <- NA
  
  ### Direct decomposition
  y <- y-min(y,na.rm = T) + 1
  
  # Smooth waveform with running mean and window of size width, using 'C' algorithm; 'fast' can't handle na
  if (smooth ==TRUE) {
    y <- runmean(y,width,'C')
    } 
  
  # Restore NAs to 0
  y[is.na(y)] <- 0
  
  # Fix problematic peaks if necessary
  if (peakfix == T) {
    y <- peakfix(y)
  }

  # Identify peaks
  #peakrecord <- lpeak(y, 3)
  
  # Find return time of peak
  #peaktime <- which(peakrecord == T)
  peaktime = peakmetrics$mu
  npeaks = length(peaktime)
  
  # If no peaks in the deconvolved waveform, take the max of the vector as only peak
  if (npeaks == 0){
    peaktime = which.max(y)
  }
  
  # Filter out noisy peaks (those less than threshold*max intensity in the return vector)
  ymax = max(y, na.rm=T)
  peak.idxs = y[peaktime] >= thres*ymax
  
  # Get the time of true peaks
  realpeak.idxs = peaktime[peak.idxs]
  
  # Get the number of true peaks
  real.npeaks <- length(realpeak.idxs)
  
  # Get the intensity of real peaks
  newpeak = peakmetrics[peakmetrics$mu %in% realpeak.idxs,]$A
  
  # Get sigma of true peaks
  sigmas <- peakmetrics[peakmetrics$mu %in% realpeak.idxs,]$sigma
  
  #then we fliter peak we have in the waveform
  #you must define newpeak as a list or a 1D vector; otherwise it's just a scalar

  ### Initialize parameters
  ### For normal Gaussian
  gu = realpeak.idxs
  gi = newpeak*0.8
  gsd = sigmas
  # gsd<-realpeak.idxs[1]/5
  # if (real.npeaks>1){
  #   gsd[2:real.npeaks]<-diff(realpeak.idxs)/4
  # }
  
  # Fit gaussians using the auto generate formula

  init0 <- gennls(A=gi, u=gu, sig=gsd)
  df<-data.frame(x=seq_along(y),y)

  log<-tryCatch(
    fit<-nlsLM(
      init0$formula,
      data=df,
      start=init0$start,
      algorithm='LM',
      control=nls.lm.control(
        factor=100,
        maxiter=1024,
        ftol = .Machine$double.eps,
        ptol = .Machine$double.eps),
      na.action=na.omit),
    error=function(e) NULL)
  
  ### Determine whether nls fit was successful
  if (!is.null(log)){
    result = summary(fit)$parameters
    pn = sum(result[,1]>0)
    rownum = nrow(result)
    npeak = rownum/3
    
    # Record the shot number of unsuccessful fits
    goodfit.idx = index
    ga = cbind(index,result)
    pmi = NULL
    
    if (pn==rownum){
      goodfit.idx = index

      # Assemble parameters into matrix
      # Make a matrix
      pm = matrix(NA,npeak,9)
      
      # Populate
      # A
      pm[,1]<-result[1:npeak,1]
      pm[,4]<-result[1:npeak,2]
      
      # mu
      s2<-npeak+1;e2<-2*npeak
      pm[,2]<-result[s2:e2,1]
      pm[,5]<-result[s2:e2,2]
      
      # sigma
      s3<-2*npeak+1;e3<-3*npeak
      pm[,3]<-result[s3:e3,1]
      pm[,6]<-result[s3:e3,2]

      # Estimate full width at half-maximum
      FWAHM = pm[,3]*2*sqrt(2*log(2))
      pm[,7] = FWAHM
      
      # Estimate front slope
      FS = peakmetrics[peakmetrics$mu %in% realpeak.idxs,]$fs
      pm[,8] = FS
      
      # Find first half-intensity bin
      HA <- peakmetrics[peakmetrics$mu %in% realpeak.idxs,]$ha
      pm[,9] = HA
     
      # Bind index to parameters
      pmi<-cbind(index,pm)
      colnames(pmi) = c('index','A','u','sigma','A_std','u_std','sig_std', 'w', 'fs', 'ha')
    }
    return (list(goodfit.idx, ga, pmi))
  }
}




