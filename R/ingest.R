#' Ingest Binary Waveform Data
#' 
#' This function allows you to ingest waveform data in the ENVI binary format provided by NEON Airborne Observing Platform acquisitions.
#' @param flightpath A string pointing to the filepath of a directory containing one waveform binary package. The package should typically contain, at minimum, header (.hdr) and binary data files describing:
#' - observation parameters 
#' - geolocation
#' - outgoing pulse
#' - return pulse
#' - impulse response
#' - outgoing impulse response
#'
#' @return A list with five components, `out`, `re`, `geol`, `sysir`, and `outir`, each representing the arrays contained in the binary package

#'#' @examples
#' wf <- ingest('./Data/2018_CRBU_1_2018061214_FL002-001')
#' for w in wf:
#'     print(w[1:10])

#' @keywords waveform

readbinary_ <- function(datafile, headerfile) {
  out = tryCatch(
    {
      read.ENVI(datafile, headerfile = headerfile)
    },
    error = function(cond) {
      message(paste('Data or header file seems to be missing:', datafile[1]))
      message('This is the original error message:')
      message(cond)
      # Specify return value in case of error
      return(NULL)
    },
    warning = function(cond) {
      message(paste('Ingesting file', datafile[1], 'caused a warning:'))
      message('This is the original warning message:')
      message(cond)
      # Specify return value in case of warning
      return(NULL)
    }
  )
}

#' @export
ingest <- function(flightpath){
  
  # Name the waveform files to ingest
  obs_bin = grep(list.files(flightpath, full.names = T), # Observation
                 pattern = 'observation', 
                 value = T)
  geo_bin = grep(list.files(flightpath, full.names = T), # Geolocation array
                 pattern = 'geolocation', 
                 value = T)
  out_bin = grep(list.files(flightpath, full.names = T), # Outgoing pulse
                 pattern = 'outgoing', 
                 value = T)
  re_bin = grep(list.files(flightpath, full.names = T), # Return pulse
                pattern = 'return_pulse', 
                value = T)
  imp_re_bin = grep(list.files(flightpath, full.names = T), # System impulse response (for deconvolution)
                    pattern = 'impulse_response_', 
                    value = T)
  imp_out_bin = grep(list.files(flightpath, full.names = T), # Outgoing impulse response (for deconvolution)
                     pattern = 'impulse_response_T0', 
                     value = T)
  
  # Read the binary files as arrays
  obs_array = readbinary_(obs_bin[1], headerfile = obs_bin[2])
  out_array = readbinary_(out_bin[1], headerfile = out_bin[2])
  geo_array = readbinary_(geo_bin[1], headerfile = geo_bin[2])
  re_array = readbinary_(re_bin[1], headerfile = re_bin[2])
  imp_re_array = readbinary_(imp_re_bin[1], headerfile = imp_re_bin[2])
  imp_out_array = readbinary_(imp_out_bin[1], headerfile = imp_out_bin[2])
  
  # Load return as reshaped data table
  out = data.table(index=c(1:nrow(out_array)), out_array)
  re = data.table(index=c(1:nrow(re_array)), re_array)
  geol = data.table(index=c(1:nrow(geo_array)), geo_array)
  sysir = data.table(index=c(1:nrow(imp_re_array)), imp_re_array)
  outir = data.table(index=c(1:nrow(imp_re_array)), imp_re_array)
  
  # Assign new geo column names to work with downstream functions
  geoindex = c(1:9,16)
  colnames(geol)[geoindex] = c(
    'index', 
    'orix', 
    'oriy', 
    'oriz', 
    'dx', 
    'dy', 
    'dz', 
    'outref', 
    'refbin', 
    'outpeak')
  
  # Return values
  wf_arrays = list(
    'out' = out, 
    're' = re, 
    'geol' = geol, 
    'sysir' = sysir, 
    'outir' = outir)
  
  return(wf_arrays)
  
}
