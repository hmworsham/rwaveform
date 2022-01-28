#' QA and Clean Results of Deconvolution
#' 
#' This function finds NaN returns in a set of deconvolved waveforms and infills with 0 or mean of surrounding values
#' @param de is a data.table of deconvolved waveform results with popential errors
#'
#' @return The cleaned, deconvolved return waveform
#' @import data.table
#' @keywords waveform, deconvolution, Gold, Richardson-Lucy
#' @export 

deconv.clean <- function(de) {
  
  # Clean up NaNs
  nanrows = which(rowSums(is.na(de[,2:length(de)]))>0)
  for (nr in nanrows){
    hasnan = de[nr]
    if(is.nan(sum(hasnan))){
      nan.times = list(which(apply(hasnan, 1, is.nan)))
      nanline = paste('Waveform', nr, 'returned NaN at time:', nan.times)
      write(nanline, file="deconvlog.txt", append=TRUE)
      for (n in nan.times) {
        if (n>2 & n < 500) {
          hasnan[[n]] = mean(c(hasnan[[n+1]], hasnan[[n-1]]))
        } else if (n < 500) {
          hasnan[[n]] = hasnan[[n+1]]
        } else {
          hasnan[[n]] = hasnan[[n-1]]
        }
      }
    }
    de[nr] = hasnan
  }
  
  # Clean up extremely unrealistic values
  bigrows = which(rowSums(de[,2:length(de)])>10^5)
  for (br in bigrows){
    hasbig = de[br]
    rescale <- function(x){
      if (x!=0){
        basis = ceiling(log10(x))
        rs = x*10^-(basis+3)
      } else {
        rs = 0
      }
      return(rs)
    }
    hasbig = data.table(t(apply(hasbig, 1, rescale)))
    de[br] = hasbig
  }
  
  return(de)
  }
  