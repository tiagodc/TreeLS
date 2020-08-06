#' @examples
#' file = system.file("extdata", "pine_plot.laz", package="TreeLS")
#' tls = readTLS(file) %>%
#'   tlsNormalize %>%
#'   tlsSample
#'
#' x = plot(tls)
#'
#' map = treeMap(tls, map.hough(h_step = 1, max_h = 4))
#' add_treeMap(x, map, color='red')
#'
#' xymap = treeMap.positions(map)
