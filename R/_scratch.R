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