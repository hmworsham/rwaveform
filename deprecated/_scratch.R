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


# DECONVOLUTION
# Remove the index columns -- we'll replace them later
#decon = subset(decon, select = -index) #****
#geol = subset(geol, select = -index)

## Convert the arrays to lists for batch deconvolution
decon2 = lapply(as.list(as.data.frame(t(decon))), as.numeric)
#geol2 = lapply(as.list(as.data.frame(t(geol))), as.numeric)

# Run adaptive decomposition algorithm on clipped returns
# Use error handling to identify erroneous or un-decomposable returns
safe_decomp = function(x){
  tryCatch(decom.adaptive(x, peakfix=F, smooth = T, thres = 0.75, width = 3), 
           error = function(e){NA})}

# Apply safe decomposition to the set
decomp = pbmcmapply(
  safe_decomp,
  decon2,
  mc.cores=getOption("mc.cores", ceiling(detectCores()/2))
)

# Filter out returns that threw exceptions
successes = which(!is.na(decomp) & lengths(decomp)>0) #****
decomp2 = decomp[successes] #****
geol2 = geol[successes] #****


# Pull Gaussian parameters
rfit = do.call('rbind', lapply(decomp2, '[[', 1)) # Indices of correctly processed waveforms
gpars = do.call('rbind', lapply(decomp2, '[[', 3)) # The Gaussian parameters

# Get indices of waveforms that need to be reprocessed
problem_wfs = setdiff(as.numeric(decon$index), successes)
problem_wfs = c(81,99,205)
decon3 = lapply(as.list(as.data.frame(t(decon[problem_wfs]))), as.numeric)

decomp.prob = pbmcmapply(
  safe_decomp,
  decon3,
  mc.cores=getOption("mc.cores", ceiling(detectCores()/2))
)

dcp.ind <- do.call('rbind', lapply(decomp.prob, '[[', 1))
dcp.rp <- do.call('rbind', lapply(decomp.prob, '[[', 3))

rbind(rfit, dcp.ind)
rbind(gpars, dcp.rp)

# Preserve parameters that resulted from successful decomposition
repars = gpars[!is.na(gpars[,1]),]
colnames(repars) = c('index', 'A', 'u', 'sigma', 'r', 'A_std', 'u_std', 'sig_std', 'r_std')
geolcols <- c(1:9,16)
colnames(geol2)[geolcols] <- c('index', 'orix', 'oriy', 'oriz', 'dx', 'dy', 'dz', 'outref', 'refbin', 'outpeak')

geol

return(list('repars' = repars, 'geolocation' = geol))





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