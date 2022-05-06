#' Individual Tree Detection
#' 
#' This function applies a given tree detection algorithm to find individual trees in a LAS/non-LAS point cloud
#' 
#' @param predictrees is a a dataframe of predicted trees with ID, X, Y, and Z values
#' @param aoi is a string indicating a plot ID for subset processing
#' 
#' @return a value describing total loss in matching field to predicted trees over an area of interest
#' 
#' @examples
#' sgtrees <- itcdelineate(testlas, aoi)
#'
#' @keywords lidar, point cloud, individual tree detection, loss, RMSE

#' @export
predictloss <- function(inventorydir, shapedir, predictrees, aoi) {
  
  # Get plot boundary for aoi
  plotpath = list.files(shapedir, 
                        pattern = glob2rx(paste0(aoi,"*shp")),
                        full.names = T)
  plotsf = readOGR(plotpath, verbose = F)
  geoextent = as.list(extent(plotsf))
  
  
  invfiles = list.files(fidir, pattern = paste0(aoi,'_inventory_data_202'), recursive=T, full.names = T)
  inv = read_excel(invfiles[1], sheet='inventory_data')
  df = data.frame('Z'=as.numeric(inv$Height_Avg_M), 'X'=as.numeric(inv$Longitude), 'Y'=as.numeric(inv$Latitude))
  df = na.omit(df)
  
  st_crs(predictrees) <- '+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' 
  y = st_as_sf(df, coords = c('X', 'Y'), crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
  y = st_transform(y, '+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
  
  plot(plotsf, axes=T)
  plot(st_geometry(y), pch='+', col='lightblue4', add=T)
  plot(st_geometry(predictrees), pch='+', add=TRUE, col='firebrick2')
  legend('topright', legend=c('Field trees', 'Modeled trees'), col=c('lightblue4', 'firebrick2'), pch='+', cex=0.6)
  
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