#' @examples
#' ### single tree
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' tls = stemPoints(tls)
#' df = stemSegmentation(tls)
#'
#' head(df)
#' tlsPlot(tls, df)
#'
#' ### forest plot
#' file = system.file("extdata", "pine_plot.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' # normalize the point cloud
#' tls = tlsNormalize(tls)
#'
#' # map the trees on a resampled point cloud so all trees have approximately the same point density
#' thin = tlsSample(tls, voxelize(0.02))
#' map = treeMap(thin, map.hough(min_density = 0.05))
#'
#' tls = stemPoints(tls, map)
#' df = stemSegmentation(tls, sgmt.ransac.circle(n=10))
#'
#' head(df)
#' tlsPlot(tls, df, map)
