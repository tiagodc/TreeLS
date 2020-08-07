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
#' and none explained by the third eigen value. The 'normal' of such surface,
#' i.e. the vector oriented in the direction orthogonal to the surface, 
#' is therefore represented by the third eigenvector.
#' 
#' Methods for both tree mapping and stem segmentation use those metrics, so in order 
#' to speed up the workflow one might apply \code{\link{fastPointMetrics}} to the point
#' cloud before other methods. The advantages of this approach are that users
#' can parameterize the point neighborhoods themselves when calculating their metrics. 
#' Those calculations won't be performed again internally in the tree mapping or stem 
#' denoising methods, reducing the overall processing time.
