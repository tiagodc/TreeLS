#ifndef CLASSES_HPP
#define CLASSES_HPP

#include <vector>

using namespace std;

class Raster{
  public:
    vector< vector<int> > matrix;
    int x_dim;
    int y_dim;
    int max_count;
    float pixel_size;
    double thickness;
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;

    void setMatrixSize(int nx, int ny){
      matrix.resize(nx, vector<int>(ny));
    };
};

class HoughCircle{
  public:
    float x_center;
    float y_center;
    double radius;
    int n_votes;
};

class HoughCenters{
  public:
    vector<HoughCircle> circles;
    int main_circle;
    double avg_x;
    double avg_y;
    float aggregate_radius;
    double low_z;
    double up_z;
    unsigned int tree_id = 0;

    void getCenters(){
      double ax = 0, ay = 0;
      main_circle = 0;
      for(unsigned i = 0; i < circles.size(); ++i){
        ax += circles[i].x_center;
        ay += circles[i].y_center;

        if(circles[i].n_votes > circles[main_circle].n_votes){
          main_circle = i;
        }
      }
      avg_x = ax / circles.size();
      avg_y = ay / circles.size();
    }

    HoughCircle getMainCircle(){
      return circles[main_circle];
    }
};

#endif // CLASSES_HPP
