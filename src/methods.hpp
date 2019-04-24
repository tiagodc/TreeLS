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

#include "algorithms.hpp"

Raster getCounts(vector<vector<double> >& slice, double pixel_size);

vector<HoughCenters> getCenters(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

HoughCenters getSingleCenter(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3);

void assignTreeId(vector<HoughCenters>& disks, double distmax, double countDensity, unsigned int minLayers=1);

vector<HoughCenters> treeHough(vector<vector<double> >& cppCloud, double h1 = 1, double h2 = 3, double hstep=0.5, double radius=0.25, double pixel=0.025, double density=0.1, unsigned int votes=3);

vector< vector<double> > ransacStemCircles(vector<vector<double> >& cloud, std::vector<unsigned int>& segments, std::vector<double>& radii, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05);

vector<vector<vector<double> > > ransacPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05);

#endif // METHODS_HPP
