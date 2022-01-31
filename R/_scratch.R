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
    message(paste('File caused a warning:', datafile[1]))
    message('This is the original warning message:')
    message(cond)
    # Specify return value in case of warning
    return(NULL)
  }
)





# tdecon <- t(decon)
# 
# split_df <- function(d, var) {
#   base::split(d, get(var, as.environment(d)))
# }
# 
# safe_decomp2 <- function(dt) {
#   safe_decomp = function(x){
#     tryCatch(decom.adaptive(x, smooth = T, peakfix=F, thres = 0.2, width = 3), 
#              error = function(e){NA})}
#   
#   dtnew = apply(
#     dt,
#     1,
#     safe_decomp
#   )
#   
#   return(dtnew)
# }
# 
# decom.parapply <- function(dt) {
#   dt2 <- split_df(dt, 'index')
#   
#   require(parallel)
#   cl <- parallel::makeCluster(min(nrow(dt), parallel::detectCores()))
#   clusterExport(cl, varlist= 'safe_decomp2')
#   clusterExport(cl, varlist= "dt2", envir = environment())
#   clusterEvalQ(cl, library("data.table"))
#   
#   dt2 <- parallel::parLapply(cl, X=dt2, fun=safe_decomp2)
#   
#   parallel::stopCluster(cl)
#   return(dt2)
# }
# 
# decomp = decom.parapply(decon)