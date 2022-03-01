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
#' @import caTools
#'
#' @examples
#' wf <- ingest('./Data/2018_CRBU_1_2018061214_FL002-001')
#' for w in wf:
#'     print(w[1:10])

#' @keywords waveform

readbinary_ <- function(file) {
  
  # Define header file and binary data file
  datafile = file[1]
  headerfile = file[2]
  out = read.ENVI(datafile, headerfile = headerfile)
  
  return(out)
}

cleanup_ <- function(wf_arrays) {
  
  # Define waveform data as reshaped data table for output
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
  
  # Return values as arrays in named list
  wf_arrays = list(
    'out' = out,
    're' = re,
    'geol' = geol,
    'sysir' = sysir,
    'outir' = outir,
    'obs' = obs)
  
  return(wf_arrays)
}

#' @export
ingest <- function(flightpath){
  
  filenam = tail(unlist(strsplit(flightpath, '/')), 1)
  
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
  
  # Concatenate all of those files to a list
  wf_paths = list(
    obs_bin, 
    geo_bin,
    out_bin,
    re_bin,
    imp_re_bin,
    imp_out_bin)
  
  # Ingest all binary/hdr file pairs in list
  wf_arrays = mclapply(
    wf_paths,
    readbinary_)

  # Clean the arrays and adjust column names where needed
  wf_arrays <- tryCatch(
    {
      cleanup_(wf_arrays)
    },

    error = function(cond) {
      message(
        paste(
          filenam,
          ': Cleaning flightpath threw an error. The original message is:',
          cond))

      # Specify return value in case of error
      return(NULL)
    },

    warning = function(cond) {
      message(
        paste(
          filenam,
          'Cleaning flightpath threw a warning. The original message is:',
          cond))

      # Specify return value in case of warning
      return(NULL)
    }
  )

  return(wf_arrays)

}
