#ifndef CLASSES_HPP
#define CLASSES_HPP

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]

#include <boost/functional/hash.hpp>
#include <math.h>
#include <memory>
#include <Rcpp.h>
#include <vector>

using namespace std;

typedef unordered_set<array<unsigned int, 2>, boost::hash< array<unsigned int, 2> > > PixelSet;

class Raster{
  public:
    vector< vector<unsigned int> > matrix;
    unsigned int x_dim;
    unsigned int y_dim;
    unsigned int max_count;
    double pixel_size;
    double thickness;
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;

    void setMatrixSize(){
      matrix.resize(x_dim, vector<unsigned int>(y_dim));
    }

    vector<double> absCenter(unsigned int x, unsigned int y){
      double x_cen = ( min_x + (pixel_size/2) ) + ( x * pixel_size );
      double y_cen = ( min_y + (pixel_size/2) ) + ( y * pixel_size );
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

    void setDims(){
      x_dim = abs( ceil( (max_x - min_x) / pixel_size ) );
      y_dim = abs( ceil( (max_y - min_y) / pixel_size ) );
    }

};

class HoughCircle{
  public:
    double x_center;
    double y_center;
    double radius;
    int n_votes;
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
