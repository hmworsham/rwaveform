#' Tree Detection Loss
#'
#' This function compares predicted to field trees for individual-scale loss computation.
#'
#' @param predictrees is a dataframe of predicted trees with ID, X, Y, and Z values
#' @param invtrees is a dataframe of field-inventory trees including ID, X, Y, and Z values
#' @param plot.sf is an sf object delineating a field inventory plot
#' @param aoi is a string indicating a plot ID for subset processing
#'
#' @return a value describing total loss in matching field to predicted trees over an area of interest
#'
#' @examples
#' sgtrees <- itcdelineate(testlas, aoi)
#'
#' @keywords lidar, point cloud, individual tree detection, loss, RMSE

#' @export
predictloss <- function(predictrees, invtrees, plot.sf, aoi) {

  plot(st_geometry(invtrees), pch='+', col='lightblue4')
  plot(st_geometry(predictrees), pch='+', add=TRUE, col='firebrick2')
  legend('topright',
         legend=c('Field trees', 'Modeled trees'),
         col=c('lightblue4', 'firebrick2'),
         pch='+',
         cex=0.6)

  nn = st_nn(predictrees, y, k=1, returnDist=T)
  dists = unlist(nn$dist)
  nns = unlist(nn$nn)
  delta.ht = c()
  delta.dist = c()
  misses = c()

  for(el in seq(length(dists))){
    if(dists[el]<=2){
      dh = df$Z[nns[el]] - predictrees$Z[el]
      delta.ht <- c(delta.ht, dh)
      delta.dist <- c(delta.dist, dists[el])
    } else {
      misses = c(misses, el)
    }
  }

  loss = sqrt(sum(delta.ht^2, delta.dist^2))*(1-length(misses)/length(dists))

  return(loss)
}
