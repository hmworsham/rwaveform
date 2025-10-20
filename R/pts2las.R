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



#' @export
pts2las <- function(infile, outpath=NULL){
  
  # Get point cloud
  pc = saferead.csv(infile)
  pc = pc[order(pc$px, pc$py),]
  
  # Filter outliers and NAs
  #pc = pc[which(pc$pi < quantile(pc$pi, 0.95, na.rm=T)),]
  #pc = pc[which(pc$pz < quantile(pc$pz, 0.99, na.rm=T)),]
  #pc = pc[which(pc$pz > quantile(pc$pz, 0.01, na.rm=T)),]
  #pc = na.omit(pc)
  
  # Count
  pc = pc %>%
    group_by(index) %>%
    dplyr::mutate(nreturns=n()) %>%
    # dplyr::mutate(rn=row_number())
    ungroup()
  
  # Convert to las and classify ground points before normalizing
  lasdata = data.frame(X = pc$px,
                       Y = pc$py,
                       Z = pc$pz,
                       gpstime = 0,
                       Intensity = as.integer(pc$pi),
                       ReturnNumber=as.integer(pc$rn),
                       NumberOfReturns=as.integer(pc$nreturns),
                       scale=1L,
                       offset=0,
                       ScanDirectionFlag = 0L,
                       EdgeOfFlightline = 0L,
                       Classification = 0L,
                       ScanAngleRank = 0L,
                       UserData = 0L,
                       PointSourceID = 0L, 
                       t = pc$t,
                       sd = pc$sd, 
                       pise = pc$pise, 
                       tse = pc$tse, 
                       sdse = pc$sdse, 
                       peakwidth = pc$pw,
                       frontslope = pc$fs,
                       hat = pc$ha, 
                       lowx = pc$lowx,
                       lowy = pc$lowy,
                       lowz = pc$lowz,
                       upx = pc$upx,
                       upy = pc$upy,
                       upz = pc$upz,
                       uncerUXpeak = pc$uncerUXpeak,
                       uncerUYpeak = pc$uncerUYpeak,
                       uncerUZpeak = pc$uncerUZpeak,
                       uncerLXpeak = pc$uncerLXpeak,
                       uncerLYpeak = pc$uncerLYpeak,
                       uncerLZpeak = pc$uncerLZpeak,
                       uncerUXleading = pc$uncerUXleading,
                       uncerUYleading = pc$uncerUYleading,
                       uncerUZleading = pc$uncerUZleading,
                       uncerLXleading = pc$uncerLXleading,
                       uncerLYleading = pc$uncerLYleading,
                       uncerLZleading = pc$uncerLZleading,
                       rn = pc$rn,
                       nreturns = pc$nreturns
                       )
  
  # Create the las header file
  lasheader = header_create(lasdata)
  
  # Update some default las header attributes with correct values
  # wktcrs = 'PROJCS["WGS 84 / UTM zone 13N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32613"]]'
  lasheader = header_set_epsg(lasheader, 32613)
  
  lasheader[['Version Minor']] <- 4L
  lasheader[["Point Data Format ID"]] <- 6L
  
  lasheader = header_add_extrabytes(lasheader, pc$t, 't', 'waveform peak time')
  lasheader = header_add_extrabytes(lasheader, pc$sd, 'sd', 'stddev of waveform peak')
  lasheader = header_add_extrabytes(lasheader, pc$pise, 'pise', 'stderr of waveform intensity')
  lasheader = header_add_extrabytes(lasheader, pc$tse, 'tse', 'stderr of waveform peak time')
  lasheader = header_add_extrabytes(lasheader, pc$sdse, 'pise', 'stderr of waveform stddev')
  lasheader = header_add_extrabytes(lasheader, pc$pw, 'peakwidth', 'waveform peak width')
  lasheader = header_add_extrabytes(lasheader, pc$fs, 'frontslope', 'front slope')
  lasheader = header_add_extrabytes(lasheader, pc$ha, 'hat', 'leading half-amplitude time')
  
  lasheader = header_add_extrabytes(lasheader, pc$lowx, 'lowx', 'target x position - lead edge')
  lasheader = header_add_extrabytes(lasheader, pc$lowy, 'lowy', 'target y position - lead edge')
  lasheader = header_add_extrabytes(lasheader, pc$lowz, 'lowz', 'target z position - lead edge')
  lasheader = header_add_extrabytes(lasheader, pc$upx, 'upx', 'target x position - trail edge')
  lasheader = header_add_extrabytes(lasheader, pc$upy, 'upy', 'target y position - trail edge')
  lasheader = header_add_extrabytes(lasheader, pc$upz, 'upz', 'target z position - trail edge')
  
  lasheader = header_add_extrabytes(lasheader, pc$uncerUXpeak, 'uncerUXpeak', 'upper bound of CE95 of px')
  lasheader = header_add_extrabytes(lasheader, pc$uncerUYpeak, 'uncerUYpeak', 'upper bound of CE95 of py')
  lasheader = header_add_extrabytes(lasheader, pc$uncerUZpeak, 'uncerUZpeak', 'upper bound of CE95 of pz')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLXpeak, 'uncerLXpeak', 'lower bound of CE95 of px')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLYpeak, 'uncerLYpeak', 'lower bound of CE95 of py')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLZpeak, 'uncerLZpeak', 'lower bound of CE95 of pz')
  
  lasheader = header_add_extrabytes(lasheader, pc$uncerUXleading, 'uncerUXleading', 'upper bound of CE95 of lowx')
  lasheader = header_add_extrabytes(lasheader, pc$uncerUYleading, 'uncerUYleading', 'upper bound of CE95 of lowy')
  lasheader = header_add_extrabytes(lasheader, pc$uncerUZleading, 'uncerUZleading', 'upper bound of CE95 of lowz')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLXleading, 'uncerLXleading', 'lower bound of CE95 of lowx')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLYleading, 'uncerLYleading', 'lower bound of CE95 of lowy')
  lasheader = header_add_extrabytes(lasheader, pc$uncerLZleading, 'uncerLZleading', 'lower bound of CE95 of lowz')
  
  # lasheader = header_add_extrabytes(lasheader, pc$rn, 'rn', 'return number')
  # lasheader = header_add_extrabytes(lasheader, pc$nreturns, 'nreturns', 'number of returns on waveform')
  
  lasheader$`X scale factor` = 0.1
  lasheader$`Y scale factor` = 1
  lasheader$`Z scale factor` = .001
  lasheader$`X offset` = 0
  lasheader$`Y offset` = 0
  lasheader$`Z offset` = 0
  
  # If outpath specified, write to outpath
  if(!is.null(outpath)){
    lasfile = str_replace(basename(infile), '.csv', '.las')
    lasfile = file.path(outpath, lasfile)
  
  # Otherwise create a temp file for output
  } else {
    lasfile = file.path(tempdir(), "temp.las")
    }
  
  # Write las file out
  write.las(lasfile, lasheader, lasdata)
  rm(pc)
  rm(lasdata)
  rm(lasheader)
  gc()
}
