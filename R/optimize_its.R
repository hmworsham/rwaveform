#' Optimize Individual Tree Segmentation
#' 
#' This function selects the optimal of multiple tree segmentation algorithms to find individual trees in a LAS/non-LAS point cloud
#' It selects algorithm and parameters that minimize the loss between field-identified and LiDAR-modeled trees
#' 
#' @param nlas is TK TK
#' @param aoi is TK TK
#' 
#' @return a dataframe of predicted trees with ID, X, Y, and Z values
#'
#' @import lidR
#' @import lidRplugins
#' @import parallel
#' 
#' @examples
#' sgtrees <- itcdelineate(testlas, aoi)
#'
#' @keywords lidar, point cloud, las, individual tree detection

#' @export



optimize.its <- function(){
  
  # Ingest point cloud at plots with a buffer
  
  
}
