#' Clip Binary Waveform Data
#' 
#' This function allows you to clip waveform data ingested from the ENVI binary format to the extent of a local area of interest.
#' 
#' @param wf An array of waveform data, from one of: return pulse, outgoing pulse, system impulse response, or outgoing impulse response
#' @param geol An array of waveform geolocation data
#' @param in_arrays A list of input arrays, including the waveform and geolocation
#' @param aoi A list describing the extent of an area of interest in the form `[xmin, xmax, ymin, ymax]`

#' 
#' @return A list with three components, `outgoing`, `return`, and `geolocation`, each corresponding to a clipped data array

#'#' @examples
#' wf <- clipwf(wfarrays, 'SG-NES2')
#' for w in wf:
#'     print(w[1:10])

#' @keywords waveform


# Generic function to clip any waveform
clipfun_ <- function(wf, geol, aoi, buff=20) {
  
  # Clip the waveforms that intersect the aoi
  waveform1 = wf[,-1]
  colnames(geol)[2:9] = c('x', 'y', 'z', 'dx', 'dy', 'dz', 'or', 'fr')
  ll = apply(waveform1, 1, wavelen)
  x = geol$x + geol$dx*(round(ll/2)-geol$fr) # use the middle point to represent the waveform position
  y = geol$y + geol$dy*(round(ll/2)-geol$fr)
  ind = which (x >= aoi[1] & x<= aoi[2] & y >= aoi[3] & y<= aoi[4])
  swaveform = wf[ind,]
  
  return(swaveform)
}

safeclipfun_ <- function(wf, geol, aoi, buff=20) {
  swaveform = tryCatch(
    {
      clipfun_(wf, geol, aoi, buff)
      print(paste('Array clipped:', names(wf)))
    },
    error = function(cond){
      message(paste('Clip failed for array:', names(wf)))
      message('This is the original error message:')
      message(cond)
      # Specify return value in case of error
      return(NULL)
    },
    warning = function(cond){
      message(paste('Clipping raised a warning:', names(wf)))
      message('This is the original error message:')
      message(cond)
      # Specify return value in case of error
      return(NULL)
    }
  )
}

#' @export
clipwf <- function(in_arrays, aoi, buff=20){
  
  # specify waveform arrays
  out = in_arrays$out
  re = in_arrays$re
  geol = in_arrays$geol
  #colnames(geol)[2:9] = c('x', 'y', 'z', 'dx', 'dy', 'dz', 'or', 'fr')
  
  # Apply the clipwf function to subset waveforms
  out_sub_tmp = clipfun_(out, geol, aoi, buff)
  re_sub_tmp = clipfun_(re, geol, aoi, buff)
  geol_sub_tmp = clipfun_(geol, geol, aoi, buff)
  
  # Trim the results to the same indexes to make sure they match
  out_sub  = out_sub_tmp[out_sub_tmp$index %in% geol_sub_tmp$index]
  out_sub = out_sub[out_sub$index %in% re_sub_tmp$index]
  re_sub = re_sub_tmp[re_sub_tmp$index %in% geol_sub_tmp$index]
  re_sub = re_sub[re_sub$index %in% out_sub$index]
  geol_sub = geol_sub_tmp[geol_sub_tmp$index %in% re_sub$index]
  geol_sub = geol_sub[geol_sub$index %in% out_sub$index]
  
  # Return a list of the clipped arrays
  clips = list('out' = out_sub, 're' = re_sub, 'geol' = geol_sub)
  
  return(clips)
}