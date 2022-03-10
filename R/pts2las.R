#' Convert Return Intensity Points to LAS
#' 
#' This function applies an iterative Gaussian decomposition to a set of waveforms of dim (MxN).
#' 
#' @param aoi is a string indicating a plot ID for subset processing
#'
#' @return a LAS point cloud dataset
#' @import rlas
#' @import lidR
#' @import parallel
#' 
#' @examples
#' testlas <- pts2las(aoi)
#'
#' @keywords lidar, point cloud, las

saferead_ <- function(x){
  tryCatch(read.csv(x, header=T), 
           error = function(cond) {
             message(paste('Reading csv failed'))
           })
}

#' @export
pts2las <- function(datadir, aoi){
  
  # Get point cloud for aoi
  pc_csv = list.files(datadir, pattern = aoi, full.names = T)
  pc <- mclapply(pc_csv, saferead_, mc.cores=getOption('mc.cores', 16))
  pc <- do.call(rbind.fill, pc)
  pc <- pc[order(pc$px, pc$py),]
  
  # Filter outliers and NAs
  pc = pc[which(pc$pi < quantile(pc$pi, 0.95)),]
  pc = pc[which(pc$pz < quantile(pc$pz, 0.99)),]
  pc = pc[which(pc$pz > quantile(pc$pz, 0.01)),]
  pc = na.omit(pc)
  
  # Convert to las and classify ground points before normalizing
  lasdata = data.frame(X = pc$px,
                       Y = pc$py,
                       Z = pc$pz,
                       gpstime = 0.0,
                       Intensity = as.integer(pc$t),
                       ReturnNumber=1L,
                       scale=1L,
                       offset=0,
                       ScanDirectionFlag = 0L,
                       EdgeOfFlightline = 0L,
                       Classification = 0L,
                       ScanAngleRank = 0L,
                       UserData = 0L,
                       PointSourceID = 0L)
  
  # Create the las header file
  lasheader = header_create(lasdata)
  
  # Update some default las header attributes with correct values
  header_set_epsg(lasheader, 32613)
  lasheader$`X scale factor` = 0.1
  lasheader$`Y scale factor` = 1
  lasheader$`Z scale factor` = 1
  lasheader$`X offset` = 0
  lasheader$`Y offset` = 0
  lasheader$`Z offset` = 0
  
  # Create a temp file for output
  lasfile = file.path(tempdir(), "temp.las")
  
  # Write las file out
  write.las(lasfile, lasheader, lasdata)
  
  # Read las file in
  newlas = readLAS(lasfile)
  
  # Classify ground points to create a normalization surface
  ws = seq(3, 9, 3)
  th = seq(0.1, 3, length.out = length(ws))
  gclas = classify_ground(newlas, algorithm = pmf(ws = ws, th = th))

  # Normalize heights to surface
  nlas = normalize_height(gclas, kriging())
  
  return(nlas)
}
