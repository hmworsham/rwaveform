#'
#'@export
peak.info <- function(y, index, thr, minPH=2, minPW=1, step=0.49){
  
  x = seq_along(y)
  min.y <- min(y, na.rm=T)
  max.y <- max(y, na.rm=T)
  
  # Check curve is ascending
  yy <- c(min.y, y)
  xx <- c(x[1], x)
  
  # If no threshold specified, take min y
  if (missing(thr)) thr <- min(yy, na.rm=TRUE)
  
  # If step>0.5, exit
  # Smaller step=longer calculation and more accurate answer
  if (step>=0.5) stop('step must be smaller than 0.5')
  
  peak.x <- peak.y <- peak.w <- peak.sig <- peak.fs <- peak.ha <- NULL
  lev <- thr - step*minPH
  len.yy <- length(yy)
  
  repeat { ## find maximum above level lev
    lev <- lev + step*minPH
    if (lev >= max.y) break
    hi <- yy>lev
    start <- which(diff(c(F,hi))>0)
    end <- which(diff(c(hi,F))<0)
    len <- length(start) 
    if(len==0) next
    
    for (ii in 1:len){
      ## find maximum in each continuous part above lev
      x <- xx[start[ii]:end[ii]]
      y <- yy[start[ii]:end[ii]]
      i <- which.max(y)
      
      ## is local maximum already defined?
      if (is.element(x[i], peak.x)) next
      
      miny <- min(yy[max(1,start[ii]-1):min(end[ii]+1,len.yy)],
                  na.rm=TRUE)
      
      ## calculate height of Peak
      PH <- y[i]-miny
      
      ## calculate width at half peak height
      PW <- paste(as.numeric(y > (miny+PH/2)),collapse='')
      PW <- max(attr(gregexpr('1+', PW)[[1]],'match.length'))
      
      ## PW is calculated correctly only when the positions
      ## within the observed window is constant
      ## in rare cases PW might stay -1
      PW <-  if(PW>=0) abs(x[PW]-x[1])
      
      ## estimate sigma
      SIG <- PW/(2*sqrt(2*log(2)))
      
      # estimate front slope
      slope.x1 <- max(1,start[ii]-1)
      slope.x2 <- x[i]
      slope.y1 <- miny
      slope.y2 <- y[i]
      FS <- (slope.y2 - slope.y1) / (slope.x2 - slope.x1)
      
      # find half-intensity bin
      HA <- x[min(which(y>miny+PH/2))]
      
      
       if (PH >= minPH && PW >= minPW) {
        peak.x <- c(peak.x,x[i])
        peak.y <- c(peak.y,y[i])
        peak.w <- c(peak.w,PW)
        peak.sig <- c(peak.sig,SIG)
        peak.fs <- c(peak.fs, FS)
        peak.ha <- c(peak.ha, HA)
       }
       
    } ## for
  } ## repeat
  
  
  ####send parameters to a matrix
  # pm<-matrix(NA, length(peak.x), 9)
  # pm[,1] <- peak.y
  # pm[,4] <- NA
  # 
  # pm[,2] <- peak.x
  # pm[,5] <- NA
  # 
  # pm[,3] <- peak.sig
  # pm[,6] <- NA
  # 
  # pm[,7] <- peak.w
  # pm[,8] <- peak.fs
  # pm[,9] <- peak.ha
  # 
  # pmi<-cbind(index,pm)
  # colnames(pmi) = c('index', 'A', 'u', 'sigma', 'A_std','u_std', 'sigma_std', 'pw', 'fs', 'ha')
  # ga <- NA
  # return (list(index,ga,pmi))
  # 
  
  res <- data.frame(mu=peak.x, A=peak.y, w=peak.w, sigma=peak.sig, fs=peak.fs, ha=peak.ha)
  res <- res[order(res$mu),]
  return(res)

}

