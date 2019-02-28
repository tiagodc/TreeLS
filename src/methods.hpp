//  ===============================================================================
//
//  Developers:
//
//  Tiago de Conto - ti@forlidar.com.br -  https://github.com/tiagodc/
//
//  COPYRIGHT: Tiago de Conto, 2019
//
//  This piece of software is open and free to use, redistribution and modifications
//  should be done in accordance to the GNU General Public License >= 3
//
//  Use this software as you wish, but no warranty is provided whatsoever. For any
//  comments or questions on TreeLS, please contact the developer (prefereably through my github account)
//
//  If publishing any work/study/research that used the tools in TreeLS,
//  please don't forget to cite the proper sources!
//
//  Enjoy!
//
//  ===============================================================================

#ifndef METHODS_HPP
#define METHODS_HPP

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Eigen/Dense>
#include <iostream>
#include "classes.hpp"

using namespace Rcpp;

template <typename AnyType>
void debugMsg(AnyType container);

vector<vector<double*> > rmatrix2cpp(NumericMatrix& cloud);

vector<double> getMinMax(vector<vector<double*> >& xyz);

vector<vector<double*> > cropCloud(vector<vector<double*> > cloud, double xCenter, double yCenter, double len = 1, bool circle = true, bool negative = false);

vector<bool> cropCloudFilter(vector<vector<double*> > cloud, double xCenter, double yCenter, double len = 1, bool circle = true, bool negative = false);

vector<bool> voxelFilter(vector<vector<double*> >& cloud, double voxel_spacing = 0.025);

vector<vector<vector<double*> > > getSlices(NumericMatrix& cloud, double zmin = 1, double zmax=3, double zstep = 0.5);

vector<vector<vector<double*> > > getSlices(vector<vector<double*> >& cloud, double zmin = 1, double zmax=3, double zstep = 0.5);

vector<vector<vector<double*> > > getSlices(vector<vector<double*> >& cloud, vector<unsigned int>& identifier);

Raster getCounts(vector<vector<double*> >& slice, double pixel_size);

vector<HoughCenters> getCenters(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

HoughCenters getSingleCenter(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

void assignTreeId(vector<HoughCenters>& disks, double distmax, double countDensity, unsigned int minLayers=1);

vector<HoughCenters> treeHough(vector<vector<double*> >& cppCloud, double h1 = 1, double h2 = 3, double hstep=0.5, double radius=0.25, double pixel=0.025, double density=0.1, unsigned int votes=3);

vector<double> ransacCircle(vector<vector<double*> >& cloud, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8);

#endif // METHODS_HPP
