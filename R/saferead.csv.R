#'@export
saferead.csv <- function(infile){
  tryCatch(read.csv(infile, header=T), 
           error = function(cond) {
             message(paste('Reading csv failed'))
           })
}