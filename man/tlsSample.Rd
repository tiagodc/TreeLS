% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{tlsSample}
\alias{tlsSample}
\title{Samples a point cloud randomly or systematically}
\usage{
tlsSample(las, by = "voxel", val = 0.05)
}
\arguments{
\item{las}{\code{LAS} object}

\item{by}{sampling method: \emph{voxel} for systematic 3D sampling or \emph{random} for random sampling}

\item{val}{Sampling parameter value. For \code{by = 'voxel'}, \code{val} must be the voxel side length. For \code{bu = 'random'}, it must be the proportion of points to be kept - between 0 and 1.}
}
\value{
\code{LAS} object
}
\description{
Applies a random sample or voxel thinning algorithm te keep a fraction of the point cloud.
}