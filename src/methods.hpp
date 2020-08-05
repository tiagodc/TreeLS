//  ===============================================================================
//
//  Developers:
//
//  Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
//
//  COPYRIGHT: Tiago de Conto, 2020
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

#include "algorithms.hpp"

void assignTreeId(vector<HoughCenters>& disks, double distmax, double countDensity, unsigned int minLayers=1);

vector<HoughCenters> getCenters(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

Raster getCounts(vector<vector<double> >& slice, double pixel_size);

HoughCenters getSingleCenter(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

vector<vector<vector<double> > > irlsPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints=100, double tolerance=0.05);

vector<vector<vector<double> > > irlsPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints=100, double tolerance=0.05);

vector<vector<double> > irlsStemCircle(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints=0, double tolerance=0.05);

vector<vector<double> > irlsStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints=100,  double tolerance=0.05);

vector<vector<vector<double> > > ransacPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05);

vector<vector<vector<double> > > ransacPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples=10, double pConfidence=0.95, double pInliers=0.9, double tolerance=0.05);

vector< vector<double> > ransacStemCircle(vector<vector<double> >& cloud, std::vector<unsigned int>& segments, std::vector<double>& radii, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05);

vector<vector<double> > ransacStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples=10, double pConfidence=0.95, double pInliers=0.9, double tolerance=0.05);

vector<HoughCenters> treeHough(vector<vector<double> >& cppCloud, double h1 = 1, double h2 = 3, double hstep=0.5, double radius=0.25, double pixel=0.025, double density=0.1, unsigned int votes=3);

vector<vector<double> > pointMetrics(vector<vector<double> >& cloud, vector<vector<unsigned int> >& idx, vector<bool> which_metrics);

vector<vector<double> > voxelMetrics(vector<vector<double> >& cloud, vector<vector<unsigned int> >& idx, vector<bool> which_metrics);

vector<unsigned long long int> voxelIndex(vector<vector<double> >& cloud, double voxel_spacing=0.05);

vector<vector<vector<double> > > treeEigenHough(vector<vector<double> >& cppEigenCloud, vector<unsigned int>& pointId, vector<unsigned int>& segId, double voxel_size, double max_rad, bool is2d = true, bool getSpace = false);

vector<vector<vector<double> > > plotEigenHough(vector<vector<double> >& cppEigenCloud, vector<unsigned int>& pointId, vector<unsigned int>& treeId, vector<unsigned int>& segId, double voxel_size, double max_rad, bool is2d = true, bool getSpace = false);

vector<unsigned int> treeIdsFromMap(vector<vector<double> >& xy, vector<vector<double> >& xymap, vector<unsigned int> ids, double length = 2.5, bool circle = true);

vector<vector<double> > bfStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples = 10, double pConfidence = 0.99, double pInliers = 0.8, double max_angle = 30, double tolerance = 0.05);

vector<vector<vector<double> > > bfPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples=10, double pConfidence = 0.95, double pInliers = 0.8, double max_angle = 30, double tolerance = 0.05);

#endif // METHODS_HPP
