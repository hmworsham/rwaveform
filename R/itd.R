#' Individual Tree Detection
#' 
#' This function applies a given tree detection algorithm to find individual trees in a LAS/non-LAS point cloud
#' 
#' @param nlas is a ground-normalized las dataset
#' @param aoi is a string indicating a plot ID for subset processing
#' 
#' @return a dataframe of predicted trees with ID, X, Y, and Z values
#' @import rlas
#' @import lidR
#' @import parallel
#' 
#' @examples
#' sgtrees <- itcdelineate(testlas, aoi)
#'
#' @keywords lidar, point cloud, las, individual tree detection

saferead_ <- function(x){
  tryCatch(read.csv(x, header=T), 
           error = function(cond) {
             message(paste('Reading csv failed'))
           })
}

#' @export

itd <- function(nlas, aoi){
  
  # Get plot boundary for aoi
  plotpath = list.files(shapedir, 
                        pattern = glob2rx(paste0(aoi,"*shp")),
                        full.names = T)
  plotsf = readOGR(plotpath, verbose = F)
  geoextent = as.list(extent(plotsf))
  
  # Run locate_trees
  f = function(x) {
    y <- 2.2 * (-(exp(-0.08*(x-2)) - 1)) + 3
    y[x < 2] <- 3
    y[x > 20] <- 7
    return(y)
  }
  
  predicted.trees = locate_trees(nlas, lmf(f))
  
  return(predicted.trees)
  
}

