#' @export
decom.adaptive<-function(x,smooth=TRUE, peakfix=FALSE, thres=0.22, width=3){
  
  y0<-as.numeric(x)
  index<-as.integer(y0[1])
  y<-y0[-1]
  y[y==0]<-NA
  
  ### Prep for direct decomposition
  y <- y-min(y,na.rm = T)+1
  
  # Smooth waveform with running mean and window of size width, using 'C' algorithm; 'fast' can't handle na
  if (smooth==TRUE) {
    y <- runmean(y,width,"C")
  }
  
  # Fix problematic peaks if necessary
  # if (peakfix == TRUE) {
  #   y <- rwaveform::peakfix(y)
  # }
  
  # Restore NAs to 0
  y[is.na(y)] <- 0
  
  # Identify peaks
  peakrecord <- lpeak(y, 3)
  
  # Find return time of peak
  peaktime <- which(peakrecord == T)
  n.peaks = length(peaktime)
  
  # Catch errors where the deconvolved waveform finds no peaks
  if (n.peaks == 0){
    peaktime <- which.max(y)
  }
  
  # Filter out noisy peaks (those less than threshold*max intensity in the return vector)
  imax <- max(y, na.rm=T)
  ind <- y[peaktime] >= thres*imax
  
  # Get the time of true peaks
  realind<-peaktime[ind]
  
  # Get the intensity of real peaks
  newpeak<-y[realind]
  
  # Get the number of true peaks
  z <- length(realind)

  ### Initialize parameters for adaptive Gaussian fitting
  
  # First try to run the standard general nls algorithm to attempt Gaussian fit
  tre<-decom(x)
  
  # If that's unsuccessful, try the adaptive approach
  if (is.null(tre[[1]])){

    agu<-0.9*realind
    agi<-newpeak*2/3
    agsd<-realind[1]/6
    if (z>1){
      agsd[2:z]<-diff(realind)/5
    }
    ari<- rep(2,z)
    
    # Fit gaussians using the auto generate formula
    init0 <- agennls(agi, agu, agsd, ari)
    sv<-as.numeric(init0$start)
    ad1<-sv*c(rep(0.4,z),rep(0.35,z),rep(0.35,z),rep(0.3,z))
    #ad1<-c(rep(60,z),rep(12,z),rep(8,z),rep(0.5,z))

    up<-sv+ad1
    low<-sv-ad1
  } else if (is.na(tre[[1]])) {
    agu<-0.9*realind
    agi<-newpeak*2/3
    agsd<-realind[1]/6
    if (z>1){
      agsd[2:z]<-diff(realind)/5
    }
    ari<- rep (2,z)

    init0 <- agennls(agi, agu, agsd, ari)
    sv<-as.numeric(init0$start);
    ad1<-sv*c(rep(0.3,z),rep(0.3,z),rep(0.3,z),rep(0.25,z))
    #ad1<-c(rep(60,z),rep(12,z),rep(8,z),rep(0.5,z))

    up<-sv+ad1
    low<-sv-ad1
  } else {
    pars<-tre[[3]]

    agi<-pars[,2]
    agu<- pars[,3]
    agsd<-pars[,4]

    ari<- rep(2,z)  ###for adaptive Gaussian function

    init0 <- agennls(agi, agu, agsd, ari)
    sv<-as.numeric(init0$start);
    ad1<-sv*c(rep(0.2,z),rep(0.2,z),rep(0.25,z),rep(0.25,z))

    up<-sv+ad1
    low<-sv-ad1
  }

  #init$formula
  #init$start
  df<-data.frame(x=seq_along(y),y)
  log<-tryCatch(
    fit<-nlsLM(
      init0$formula,
      data=df,
      start=init0$start,
      algorithm='LM',
      lower=low,
      upper=up,
      control=nls.lm.control(
        factor=100,
        maxiter=1024,
        ftol = .Machine$double.eps, 
        ptol = .Machine$double.eps),
      na.action=na.omit),
    error=function(e) NULL)
  
  ###then you need to determine if this nls is sucessful or not?
  if (!is.null(log)){
    result=summary(fit)$parameters
    pn<-sum(result[,1]>0)
    rownum<-nrow(result)
    npeak<-rownum/4
    
    #record the shot number of not good fit
    rightfit<-NA;
    ga<-matrix(NA,rownum,5);#pmi<-matrix(NA,npeak,9)
    ga<-cbind(index,result)
    pmi<-NULL
    if (pn==rownum){
      rightfit<-index

      ####directly get the parameters
      ###make a matrix
      pm<-matrix(NA,npeak,8)
      pm[,1]<-result[1:npeak,1];pm[,5]<-result[1:npeak,2]
      s2<-npeak+1;e2<-2*npeak
      pm[,2]<-result[s2:e2,1];pm[,6]<-result[s2:e2,2]
      s3<-2*npeak+1;e3<-3*npeak
      pm[,3]<-result[s3:e3,1];pm[,7]<-result[s3:e3,2]
      s4<-3*npeak+1;e4<-4*npeak
      pm[,4]<-result[s4:e4,1];pm[,8]<-result[s4:e4,2]
      pmi<-cbind(index,pm)
      colnames(pmi) = c("index","A","u","sigma","r","A_se","u_se","sigma_se","r_se")

    }
    return (list(rightfit,ga,pmi))
  }
}