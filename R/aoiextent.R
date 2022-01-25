#' Find the Extent of an Area of Interest
#' 
#' This function allows you to clip waveform data ingested from the ENVI binary format to the extent of a local area of interest.
#' 
#' 
#' @param sfid A string corresponding to a shapefile filename; can be a shorthand plot ID, e.g. 'SG-NES2'
#' @param buff An integer indicating the buffer width around the AOI to include in the clip; default is 20 m
#' 
#' @return A list describing the extent of an area of interest in the form `[xmin, xmax, ymin, ymax]`

#'#' @examples
#' wf <- aoiextent(wfarrays, 'SG-NES2')
#' for w in wf:
#'     print(w[1:10])

#' @keywords waveform, geometry

#' @export
aoiextent <- function(sfid, sfdir, buff=20) {
  
  sfpath = list.files(sfdir, 
                        pattern = glob2rx(paste0(sfid,"*shp")),
                        full.names = T)
  sf = st_read(sfpath, quiet=T)
  plotbuff = st_buffer(sf,
                       buff,
                       endCapStyle = "SQUARE",
                       joinStyle = "MITRE",
                       mitreLimit = 0.05,)
  geoextent = as.list(extent(plotbuff))
  
  return(geoextent)
}
