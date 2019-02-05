#ifndef CLASSES_HPP
#define CLASSES_HPP

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]

#include <boost/functional/hash.hpp>
#include <math.h>
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
    float pixel_size;
    float thickness;
    float min_x;
    float max_x;
    float min_y;
    float max_y;
    float min_z;
    float max_z;

    void setMatrixSize(){
      matrix.resize(x_dim, vector<unsigned int>(y_dim));
    }

    vector<float> absCenter(unsigned int x, unsigned int y){
      float x_cen = ( min_x + (pixel_size/2) ) + ( x * pixel_size );
      float y_cen = ( min_y + (pixel_size/2) ) + ( y * pixel_size );
      vector<float> xy = {x_cen, y_cen};
      return xy;
    }

    vector<unsigned int> pixPosition(float x, float y){
      unsigned int x_pix = floor( (x - min_x) / pixel_size );
      unsigned int y_pix = floor( (y - min_y) / pixel_size );
      vector<unsigned int> xy = {x_pix, y_pix};
      return xy;
    }

    PixelSet rasterCircle(float radius, float cx, float cy){

      int n_points = ceil( (2 * PI * radius) / pixel_size );
      float angle_dist = 2 * PI / n_points;
      PixelSet pixels;

      for(float i = 0; i < 2*PI; i += angle_dist){
        float x = cos(i)*radius + cx;
        float y = sin(i)*radius + cy;

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
    float x_center;
    float y_center;
    float radius;
    int n_votes;
};

class HoughCenters{
  public:
    vector<HoughCircle> circles;
    HoughCircle* main_circle;
    float avg_x;
    float avg_y;
    float aggregate_radius;
    float low_z;
    float up_z;
    unsigned int tree_id = 0;

    void getCenters(){
      float ax = 0;
      float ay = 0;
      main_circle = &circles[0];

      for(auto& i : circles){
        ax += i.x_center;
        ay += i.y_center;

        if(i.n_votes > main_circle->n_votes){
          main_circle = &i;
        }
      }
      avg_x = ax / circles.size();
      avg_y = ay / circles.size();
    }

};

#endif // CLASSES_HPP
