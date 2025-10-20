#'@export
saferead.csv <- function(infile){
  tryCatch(fread(infile), 
           error = function(cond) {
             message(paste('Reading csv failed'))
           })
}