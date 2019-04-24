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

#ifndef CLASSES_HPP
#define CLASSES_HPP

// [[Rcpp::depends(BH)]]

#include <boost/functional/hash.hpp>
#include <iostream>
#include <math.h>
#include <vector>
#include <RcppArmadillo.h>

using namespace std;

typedef unordered_set<array<unsigned int, 2>, boost::hash< array<unsigned int, 2> > > PixelSet;
typedef unordered_set<array<int, 3>, boost::hash< array<int, 3> > > VoxelSet;

class tempContainer{
  public:
    vector<bool> filter;
    vector<unsigned int> counts;
    vector<unsigned int> ids;
    vector<unsigned int> sections;
    vector<double> values;

    tempContainer(unsigned int n){
      setSizes(n);
    };

    void setSizes(unsigned int n){
      filter.resize(n, false);
      counts.resize(n, 0);
      values.resize(n, 0);
      ids.resize(n, 0);
      sections.resize(n, 0);
    };

    void clear(){
      filter.clear();
      filter.shrink_to_fit();
      counts.clear();
      counts.shrink_to_fit();
      ids.clear();
      ids.shrink_to_fit();
      sections.clear();
      sections.shrink_to_fit();
      values.clear();
      values.shrink_to_fit();
    }
};

class Raster{
  public:
    vector< vector<unsigned int> > matrix;
    unsigned int x_dim;
    unsigned int y_dim;
    unsigned int max_count;
    double pixel_size;
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;

    void setDims(){
      x_dim = abs( ceil( (max_x - min_x) / pixel_size ) );
      y_dim = abs( ceil( (max_y - min_y) / pixel_size ) );
    }

    void setMatrixSize(){
      matrix.resize(x_dim, vector<unsigned int>(y_dim));
    }

    vector<double> absCenter(unsigned int x, unsigned int y){
      double x_cen = min_x + (pixel_size/2) + ( x * pixel_size );
      double y_cen = min_y + (pixel_size/2) + ( y * pixel_size );
      vector<double> xy = {x_cen, y_cen};
      return xy;
    }

    vector<unsigned int> pixPosition(double x, double y){
      unsigned int x_pix = floor( (x - min_x) / pixel_size );
      unsigned int y_pix = floor( (y - min_y) / pixel_size );
      vector<unsigned int> xy = {x_pix, y_pix};
      return xy;
    }

    PixelSet rasterCircle(double radius, double cx, double cy){

      int n_points = ceil( (2 * PI * radius) / pixel_size );
      double angle_dist = 2 * PI / n_points;
      PixelSet pixels;

      for(double i = 0; i < 2*PI; i += angle_dist){
        double x = cos(i)*radius + cx;
        double y = sin(i)*radius + cy;

        vector<unsigned int> pxy = pixPosition(x, y);
        pixels.insert({pxy[0], pxy[1]});
      }

      return pixels;

    }

    void updateMatrix(double x, double y){
      if(x < min_x || x > max_x || y < min_y || y > max_y) return;

      vector<unsigned int> xy = pixPosition(x,y);

      if(xy[0] >= x_dim || xy[1] >= y_dim) return;

      if(++matrix[ xy[0] ][ xy[1] ] > max_count) max_count = matrix[ xy[0] ][ xy[1] ];
    }

    void cleanRadius(double x, double y, double radius){

      max_count = 0;

      if(x < min_x || x > max_x || y < min_y || y > max_y) return;

      vector<unsigned int> pix = pixPosition(x,y);
      vector<double> pixCenter = absCenter(pix[0], pix[1]);
      x = pixCenter[0];
      y = pixCenter[1];

      for(unsigned int i = 0; i < matrix.size(); ++i){
        for(unsigned int j = 0; j < matrix[i].size(); ++j){

          unsigned int* ct = &matrix[i][j];

          if(*ct == 0) continue;

          vector<double> pxy = absCenter(i, j);
          double dist = sqrt( pow( pxy[0] - x ,2) + pow( pxy[1] - y ,2) );

          if(dist > radius)
            *ct = 0;
          else if(*ct > max_count)
            max_count = *ct;

        }
      }
    }

};

class HoughCircle{
  public:
    double x_center;
    double y_center;
    double radius;
    unsigned int n_votes;
};

class HoughCenters{
  public:
    vector<HoughCircle> circles;
    HoughCircle main_circle;
    double avg_x;
    double avg_y;
    double aggregate_radius;
    double low_z;
    double up_z;
    unsigned int tree_id = 0;

    void getCenters(){

      if(circles.empty())
        return;

      double ax = 0;
      double ay = 0;
      main_circle = circles[0];

      for(auto& i : circles){
        ax += i.x_center;
        ay += i.y_center;

        if(i.n_votes > main_circle.n_votes){
          main_circle = i;
        }
      }
      avg_x = ax / circles.size();
      avg_y = ay / circles.size();
    }

};

#endif // CLASSES_HPP
