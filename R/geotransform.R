#' geotransform
#'
#' The function allows you to convert parameters from decomposition to point cloud by using georeference data (generally it should come with waveform data).
#'   Detailed description of how to calculate can refer to Zhou, T., Popescu, S.C., Krause, K., Sheridan, R.D., Putman, E., 2017.
#'   Gold - a noveldeconvolution algorithm with optimization for waveform LiDAR processing.ISPRS Journal of Photogrammetry and Remote Sensing 129 (2017): 131-150.
#'   For the direct decomposition method, three kinds of position are calculated: leading edge, peak and trail edge. For the deconvolution and decomposition method,
#'   only the peak position is used to calculate the position. In additional, the uncertianty of the point cloud was provided based on detected peaks'
#'   95% confidence interval using deterministic method for waveform decomposition. This fucntion is suitable for NEON waveform lidar. For other kinds of dataset,
#'   you need to preprocess data to the same format or have the same required information or datsets.
#'
#' @param decomp the object from the decomposition or after deconvolution and decomposition.
#' @param geo the reference geolocation that is generally coming with waveform data and provided by the data provider.
#' @param source is determined by input data. If estimated parameters are from dcomposition, source = "decomposition". Otherwise will be deconvolution and decomposition.
#'   Default is decomposition.
#' @return A dataframe with columns. For the direct decompostion method, we will have 29 columns and for the deconvolution and decomposition method.
#'   17 columns were generated: index,pi,t,sd,pise,tse,sdse,px,py,pz,uncerUXpeak,uncerUYpeak,uncerUZpeak,uncerLXpeak,uncerLYpeak,uncerLZpeak,rn.
#' \item{\code{index}}{The index of waveform.}
#' \item{\code{pi}}{The estimated amplitude of an waveform componment.}
#' \item{\code{t}}{The estimated peak location of an waveform componment.}
#' \item{\code{sd}}{The estimated echo width of an waveform componment.}
#' \item{\code{pise}}{The standard error of the estimated amplitude.}
#' \item{\code{tse}}{The standard error of the estimated peak location.}
#' \item{\code{sdse}}{The standard error of the estimated echo width.}
#' \item{\code{px}}{Desired x position using peak locations.}
#' \item{\code{py}}{Desired y position using peak locations.}
#' \item{\code{pz}}{Desired y position using peak locations.}
#' \item{\code{lowx}}{Desired x position using leading edge locations.}
#' \item{\code{lowy}}{Desired y position using leading edge locations.}
#' \item{\code{lowz}}{Desired z position using leading edge locations.}
#' \item{\code{upx}}{Desired x position using trailing edge locations.}
#' \item{\code{upy}}{Desired y position using trailing edge locations.}
#' \item{\code{upz}}{Desired z position using trailing edge locations.}
#' \item{\code{uncerUXpeak}}{Upper bound of 95th confidence interval of px.}
#' \item{\code{uncerUYpeak}}{Upper bound of 95th confidence interval of py.}
#' \item{\code{uncerUZpeak}}{Upper bound of 95th confidence interval of pz.}
#' \item{\code{uncerLXpeak}}{Lower bound of 95th confidence interval of px.}
#' \item{\code{uncerLYpeak}}{Lower bound of 95th confidence interval of py.}
#' \item{\code{uncerLZpeak}}{Lower bound of 95thconfidence interval of pz.}
#' \item{\code{uncerUXleading}}{Upper bound of 95th confidence interval of lowx.}
#' \item{\code{uncerUYleading}}{Upper bound of 95th confidence interval of lowy.}
#' \item{\code{uncerUZleading}}{Upper bound of 95th confidence interval of lowz.}
#' \item{\code{uncerLXleading}}{Lower bound of 95th confidence interval of lowx.}
#' \item{\code{uncerLYleading}}{Lower bound of 95th confidence interval of lowy.}
#' \item{\code{uncerLZleading}}{Lower bound of 95th confidence interval of lowz.}
#' \item{\code{rn}}{The number of return for each waveform.}
#' By combining xyz, the users can get waveform-based point cloud using leading edge, peak and trail edge positions, and parameter uncertainty of the point cloud.
#' @import data.table
#' @import splitstackshape
#' @import sqldf
#' @importFrom stats ave
#' @export
#' @references
#'   Zhou, Tan*, Sorin C. Popescu, Keith Krause, Ryan D. Sheridan, and Eric Putman, 2017. Gold-A novel deconvolution algorithm with
#'   optimization for waveform LiDAR processing. ISPRS Journal of Photogrammetry and Remote Sensing 129 (2017):
#'   131-150. https://doi.org/10.1016/j.isprsjprs.2017.04.021
#' @examples
#'
#' data(geo)
#' data(decom_result)
#' data(decon_result)
#'##used part of data to show the results
#' decomp<-decom_result[1:80,]
#' geo<-geo[1:80]
#' deconp<-decon_result[1:80]
#'
#' ##the follwoing steps are reuired to conduct the geotransformation,
#' ##we need assign exactly same column names
#' geoindex=c(1:9,16)
#' colnames(geo)[geoindex]<-c("index","orix","oriy","oriz","dx","dy","dz","outref","refbin","outpeak")
#' ##for decomposition results
#'
#' geor<-geotransform(decomp,geo)
#'
#' ########for deconvolution and decomposition
#' decongeo<-geotransform(decomp=deconp,geo,source="deconvolution")
#'


geotransform <- function(decomp, geo, source="decomposition"){
  
  # Get georeference data
  decomp = data.table(decomp)

  # Check that indices in both dataframes are a set
  dc.id = unique(decomp$index)
  geo = geo[geo$index %in% dc.id,]
  
  ge.id = unique(geo$index)
  decomp = decomp[decomp$index %in% ge.id,]

  # Create a geolocation dataframe consistent with our decomposition data
  rindex<-table(decomp$index)
  geo0<-data.frame(geo, rindex)
  ngeo<-expandRows(geo0, "Freq")
  
  # Assign variables from geolocation columns
  orix<-ngeo$orix
  oriy<-ngeo$oriy
  oriz<-ngeo$oriz
  dx<-ngeo$dx
  dy<-ngeo$dy
  dz<-ngeo$dz
  refbin<-ngeo$refbin

  ## outref and outpeak become the reference points for decon+decom results
  outref<-ngeo$outref
  outpeak<-ngeo$outpeak

  # Find corresponding decomposition parameters file
  nr<-nrow(decomp)
  peak<-decomp$u
  sd<-decomp$sigma
  tse<-decomp$u_std
  low<-peak-sd*sqrt(2*log(2))
  up<-peak+sd*sqrt(2*log(2))

  ###4 begin to calculate

  #################to calculate the peak, leading and trail edge position
  if (source=="decomposition"){
    sgeo<-matrix(NA,nr,32)
    sgeo[,1:10]<-cbind(decomp$index,decomp$A,decomp$u,decomp$sigma,decomp$A_std,decomp$u_std,decomp$sigma_std,decomp$pw,decomp$fs,decomp$ha)
    sgeo[,11]<-(peak-refbin)*dx+orix
    sgeo[,12]<-(peak-refbin)*dy+oriy
    sgeo[,13]<-(peak-refbin)*dz+oriz
    sgeo[,14]<-(low-refbin)*dx+orix
    sgeo[,15]<-(low-refbin)*dy+oriy
    sgeo[,16]<-(low-refbin)*dz+oriz
    sgeo[,17]<-(up-refbin)*dx+orix
    sgeo[,18]<-(up-refbin)*dy+oriy
    sgeo[,19]<-(up-refbin)*dz+oriz

    #### to calculate the uncertainty of using peak,
    
    ####upper peak
    sgeo[,20]<-sgeo[,11] + 1.96*tse*dx;
    sgeo[,21]<-sgeo[,12] + 1.96*tse*dy;
    sgeo[,22]<-sgeo[,13] + 1.96*tse*dz;
    
    ####lower peak
    sgeo[,23]<-sgeo[,11] - 1.96*tse*dx;
    sgeo[,24]<-sgeo[,12] - 1.96*tse*dy;
    sgeo[,25]<-sgeo[,13] - 1.96*tse*dz;

    ####upper leading edge
    sgeo[,26]<-sgeo[,14] + 1.96*tse*dx;
    sgeo[,27]<-sgeo[,15] + 1.96*tse*dy;
    sgeo[,28]<-sgeo[,16] + 1.96*tse*dz;
    
    ####lower leading edge
    sgeo[,29]<-sgeo[,14] - 1.96*tse*dx;
    sgeo[,30]<-sgeo[,15] - 1.96*tse*dy;
    sgeo[,31]<-sgeo[,16] - 1.96*tse*dz;

    #
    ind<-decomp$index
    rn<-ave(ind, ind, FUN = seq_along)
    sgeo[,32]<-rn
    colnames(sgeo)<- c("index","pi","t","sd","pise","tse","sdse","pw","fs","ha","px","py","pz","lowx","lowy","lowz","upx","upy","upz",
                       "uncerUXpeak","uncerUYpeak","uncerUZpeak","uncerLXpeak","uncerLYpeak","uncerLZpeak",
                       "uncerUXleading","uncerUYleading","uncerUZleading","uncerLXleading","uncerLYleading","uncerLZleading","rn")
  } else {
    sgeo<-matrix(NA,nr,20)
    sgeo[,1:10]<-cbind(decomp$index,decomp$A,decomp$u,decomp$sigma,decomp$A_std,decomp$u_std,decomp$sigma_std,decomp$pw,decomp$fs,decomp$ha)
    sgeo[,11]<-(peak-refbin-outpeak+outref)*dx+orix
    sgeo[,12]<-(peak-refbin-outpeak+outref)*dy+oriy
    sgeo[,13]<-(peak-refbin-outpeak+outref)*dz+oriz

    #### to calculate the uncertainty of using peak,
    ####upper peak
    sgeo[,14]<-sgeo[,11] + 1.96*tse*dx;
    sgeo[,15]<-sgeo[,12] + 1.96*tse*dy;
    sgeo[,16]<-sgeo[,13] + 1.96*tse*dz;
    ####lower peak
    sgeo[,17]<-sgeo[,11] - 1.96*tse*dx;
    sgeo[,18]<-sgeo[,12] - 1.96*tse*dy;
    sgeo[,19]<-sgeo[,13] - 1.96*tse*dz;
    
    # Get return number
    ind<-decomp$index
    rn<-ave(ind, ind, FUN = seq_along)
    sgeo[,20]<-rn
    
    # Assign column names
    colnames(sgeo)<-c("index","pi","t","sd","pise","tse","sdse","pw","fs","ha","px","py","pz","uncerUXpeak","uncerUYpeak","uncerUZpeak","uncerLXpeak","uncerLYpeak","uncerLZpeak","rn")

  }

  return (sgeo)

}

