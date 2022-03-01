#' Apply Functions in Parallel with Safe Exit and Error Handling
#' 
#' This function allows for generic parallel processing using other functions in the package
#' @param X is a list of iterables
#' @param FUN is an object pointing to a function to be applied in parallel
#' @param clip specifies whether to clip the waveform to a particular geometry before processing. Default is not to clip (False)
#' @param aoi is list describing the extent of an area of interest in the form `[xmin, xmax, ymin, ymax]`
#'
#' @return a data.frame of points with return parameters and geolocation information
#' @import parallel
#' @examples
#' 
#' @keywords parallelization

message_parallel_ <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

#' @export 
safe_mclapply <- function(X, FUN, mc.cores, stop.on.error=T, ...) {
  fun <- function(x){
    inner_result <- tryCatch({
      withCallingHandlers(
        expr = {
          FUN(x, ...)
        },
        warning = function(e) {
          message_parallel_(trimws(paste('[WARNING]', ':', e)))
          invokeRestart('muffleWarning')
          return(NULL)
        },
        error = function(e) {
          message_parallel_(trimws(paste0('[ERROR]', x, ':', e)))
          return(NULL)
        }
      )},
      error = function(e){
        # error is returned gracefully; other results of this core won't be affected
        return(e)
      }
      )
    return(inner_result)
  }
  
  res = mclapply(X, fun, mc.cores=mc.cores)
  failed = sapply(res, inherits, what="error")
  if (any(failed == T)){
    error_indices <- paste0(which(failed == T), collapse=", ")
    error_traces <- paste0(lapply(res[which(failed == T)], function(x) x$message), collapse="\n\n")
    error_message <- sprintf("Elements with following indices failed with an error: %s. Error messages: \n\n%s", 
                             error_indices,
                             error_traces)
    if (stop.on.error)
      stop(error_message)
    else
      warning(error_message, "\n\n### Errors will be ignored ###")
  }
  return(res[!failed])
}




#' @export
safe_mcmapply <- function(X, FUN, mc.cores, stop.on.error=T, ...) {
  X = X
}
  