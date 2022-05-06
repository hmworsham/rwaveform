#'@export

# Function to normalize las points to interpolated ground surface
ground.normalize <- function(laspath, dh0=0.5, dhmax=3) {
  
  # Read las
  newlas = readLAS(laspath)

  # Classify ground points to create a normalization surface
  ws = seq(3, 9, 3)
  th = seq(0.1, 3, length.out = length(ws))
  zhangparams = util_makeZhangParam(dh0=0.5, dhmax=1.3)
  ws=zhangparams$ws
  th=zhangparams$th
  gclas = classify_ground(newlas, algorithm = pmf(ws, th), last_returns=T)
  
  # Normalize heights to surface
  nlas = normalize_height(gclas, tin())
  
  return(nlas)
}