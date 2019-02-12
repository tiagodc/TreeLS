#' @examples
#' file = system.file("extdata", "model_boles.laz", package="TreeLS")
#' tls = readTLS(file)
#' plot(tls)
#'
#' map = treeMap(tls)
#' plot(map, color='Radii')
#'
#' xymap = treePositions(map)
#' head(xymap)
