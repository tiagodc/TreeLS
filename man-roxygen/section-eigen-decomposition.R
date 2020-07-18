#' @section Eigen Decomposition of Point Neighborhoods:
#'
#' Point filtering/classification methods that rely on eigen
#' decomposition rely on shape indices calculated for point 
#' neighborhoods (knn or voxel). To derive these shape indices, eigen 
#' decomposition is performed on the XYZ columns of a point cloud patch. 
#' Metrics related to object curvature are calculated upon ratios of the resulting 
#' eigen values, and metrics related to object orientation are caltulated from 
#' approximate normals obtained from the eigen vectors.
#'
#' For instance, a point neighborhood that belongs to a perfect flat
#' surface will have all of its variance explained by the first two eigen values, 
#' and none explained by the third eigen value. The 'normal' of such a surface,
#' i.e. the vector oriented in the direction orthogonal to the surface, 
#' is therefore represented by the third eigenvector.