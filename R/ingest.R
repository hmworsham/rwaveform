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
#' @import progress
#' @import caTools
#'
#' @examples
#' wf <- ingest('./Data/2018_CRBU_1_2018061214_FL002-001')
#' for w in wf:
#'     print(w[1:10])

#' @keywords waveform

readbinary_ <- function(file) {
  datafile = file[1]
  headerfile = file[2]
  
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
  
  wf_paths = list(
    obs_bin, 
    geo_bin,
    out_bin,
    re_bin,
    imp_re_bin,
    imp_out_bin)
  
  # Read the binary files as arrays
  # pb$tick()
  # obs_array = readbinary_(obs_bin[1], headerfile = obs_bin[2])
  # pb$tick()
  # out_array = readbinary_(out_bin[1], headerfile = out_bin[2])
  # pb$tick()
  # geo_array = readbinary_(geo_bin[1], headerfile = geo_bin[2])
  # pb$tick()
  # re_array = readbinary_(re_bin[1], headerfile = re_bin[2])
  # pb$tick()
  # imp_re_array = readbinary_(imp_re_bin[1], headerfile = imp_re_bin[2])
  # pb$tick()
  # imp_out_array = readbinary_(imp_out_bin[1], headerfile = imp_out_bin[2])
  wf_arrays = pbmclapply(wf_paths, readbinary_)
  
  #Load return as reshaped data table
  obs = data.table(index=c(1:nrow(wf_arrays[[1]])), wf_arrays[[1]])
  geol = data.table(index=c(1:nrow(wf_arrays[[2]])), wf_arrays[[2]])
  out = data.table(index=c(1:nrow(wf_arrays[[3]])), wf_arrays[[3]])
  re = data.table(index=c(1:nrow(wf_arrays[[4]])), wf_arrays[[4]])
  outir = data.table(index=c(1:nrow(wf_arrays[[5]])), wf_arrays[[5]])
  sysir = data.table(index=c(1:nrow(wf_arrays[[6]])), wf_arrays[[6]])

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
    'outir' = outir, 
    'obs' = obs)

  return(wf_arrays)
  
}
