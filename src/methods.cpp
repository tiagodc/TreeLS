#include <Rcpp.h>
#include <math.h>
#include <iostream>

#include "classes.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

vector<double> absCenter(int x, int y, double min_x, double min_y, float step){
  double x_cen = ( min_x + (step/2) ) + ( x * step );
  double y_cen = ( min_y + (step/2) ) + ( y * step );
  return vector<double> {x_cen, y_cen};
}

vector<int> pixPosition(double x, double y, double min_x, double min_y, float step){
  int x_pix = floor( (x - min_x) / step );
  int y_pix = floor( (y - min_y) / step );
  return vector<int> {x_pix, y_pix};
}

vector< vector<int> > rasterCircle(float radius, float pixel_size, double cx, double cy, double mx, double my){

  int n_points = ceil( (2 * PI * radius) / pixel_size );
  double angle_dist = 2 * PI / n_points;
  double x, y;
  vector<int> pxy(2);
  vector< vector<int> > pixels;

  for(double i = 0; i < 2*PI; i += angle_dist){
    x = cos(i)*radius + cx;
    y = sin(i)*radius + cy;

    pxy = pixPosition(x, y, mx, my, pixel_size);
    pixels.push_back(pxy);
  }

  return pixels;

}

Raster getCounts(NumericMatrix& slice , float pixel_size, float x_mid=0, float y_mid=0, float d_mid=-1){

  NumericMatrix::Column xcol = slice( _, 0);
  NumericMatrix::Column ycol = slice( _, 1);
  NumericMatrix::Column zcol = slice( _, 2);

  NumericMatrix::Column::iterator minx = std::min_element(xcol.begin(), xcol.end());
  NumericMatrix::Column::iterator miny = std::min_element(ycol.begin(), ycol.end());
  NumericMatrix::Column::iterator minz = std::min_element(zcol.begin(), zcol.end());
  NumericMatrix::Column::iterator maxx = std::max_element(xcol.begin(), xcol.end());
  NumericMatrix::Column::iterator maxy = std::max_element(ycol.begin(), ycol.end());
  NumericMatrix::Column::iterator maxz = std::max_element(zcol.begin(), zcol.end());

  Raster ras;
  ras.min_x = (d_mid > 0) ? (x_mid - d_mid) : *minx;
  ras.max_x = (d_mid > 0) ? (x_mid + d_mid) : *maxx;
  ras.min_y = (d_mid > 0) ? (y_mid - d_mid) : *miny;
  ras.max_y = (d_mid > 0) ? (y_mid + d_mid) : *maxy;
  ras.min_z = *minz;
  ras.max_z = *maxz;
  ras.pixel_size = pixel_size;
  ras.thickness = ras.max_z - ras.min_z;

  int xn = abs( ceil( (ras.max_x - ras.min_x) / pixel_size ) ) ;
  int yn = abs( ceil( (ras.max_y - ras.min_y) / pixel_size ) ) ;

  ras.x_dim = xn;
  ras.y_dim = yn;
  ras.setMatrixSize(xn, yn);

  int max_votes = 0;
  double x,y;
  int xCell;
  int yCell;

  for(int i = 0; i < slice.nrow(); i++){

    x = xcol[i];
    y = ycol[i];

    double dist = sqrt( pow(x - x_mid,2) + pow(y - y_mid,2) );
    if(d_mid > 0 && dist > d_mid) continue;

    xCell = abs( floor( (x - ras.min_x) / pixel_size ) );
    yCell = abs( floor( (y - ras.min_y) / pixel_size ) );

    if(xCell >= xn) xCell = xn-1;
    if(yCell >= yn) yCell = yn-1;

    ras.matrix[xCell][yCell] += 1;

    if(ras.matrix[xCell][yCell] > max_votes) max_votes = ras.matrix[xCell][yCell];
  }

  ras.max_count = max_votes;

  cout << "x: " << ras.matrix.size() << " y: " << ras.matrix[0].size() << endl;

  return ras;

}

vector<HoughCenters> getCenters(Raster* raster, float max_radius=0.25, float min_den=0.1, unsigned int min_votes=3){

  //count raster properties (&raster)
  unsigned int x_len = raster->x_dim;
  unsigned int y_len = raster->y_dim;
  int min_count = ceil( (raster->max_count)*min_den );

  //empty raster properties (votes)
  Raster empty_raster;
  empty_raster.min_x = raster->min_x - max_radius;
  empty_raster.min_y = raster->min_y - max_radius;
  empty_raster.max_x = raster->max_x + max_radius;
  empty_raster.max_y = raster->max_y + max_radius;
  empty_raster.pixel_size = raster->pixel_size;
  empty_raster.x_dim = abs( ceil( (empty_raster.max_x - empty_raster.min_x) / empty_raster.pixel_size ) );
  empty_raster.y_dim = abs( ceil( (empty_raster.max_y - empty_raster.min_y) / empty_raster.pixel_size ) );
  empty_raster.setMatrixSize(empty_raster.x_dim, empty_raster.y_dim);

  //get valid pixels
  vector< vector<unsigned int> > pixels;

  for(unsigned i = 0; i < x_len; ++i){
    for(unsigned j = 0; j < y_len; ++j){
      if(raster->matrix[i][j] >= min_count ){
        pixels.push_back( {i,j} );
      }
    }
  }

  //make circles centered in every valid pixel
  vector<HoughCenters> g_circles = {};
  HoughCenters g_cen;
  g_cen.circles = {};
  g_cen.avg_x = -1000;
  g_cen.avg_y = -1000;
  g_cen.low_z = raster->min_z;
  g_cen.up_z = raster->max_z;
  g_cen.aggregate_radius = max_radius * 3;
  g_circles.push_back(g_cen);

  vector< vector<int> > votes;
  HoughCircle hc;
  vector<double> center;
  vector< vector<int> > h_circle;
  set<unsigned long long int> pixel_set; //cod = 100.000*x + y
  vector<double> coor(2);
  unsigned int vx, vy;
  vector<int> pixel(2);
  vector<float> radii;

  for(float i = 0; i <= max_radius; i+=raster->pixel_size){
    radii.push_back(i);
  }

  for(auto& rad : radii){
    votes = empty_raster.matrix;
    hc.radius = rad;
    pixel_set = {};

    for(unsigned i = 0; i < pixels.size(); ++i){

      center = absCenter(pixels[i][0], pixels[i][1], raster->min_x, raster->min_y, raster->pixel_size);
      h_circle = rasterCircle(rad, empty_raster.pixel_size, center[0], center[1], empty_raster.min_x, empty_raster.min_y);

      for(unsigned j = 0; j < h_circle.size(); ++j){
        vx = h_circle[j][0];
        vy = h_circle[j][1];

        if( vx >= votes.size() || vy >= votes[0].size() ) continue;

        votes[ vx ][ vy ] += 1;

        if(votes[ vx ][ vy ] >= min_votes){
          pixel_set.insert( 100000*vx + vy );
        }
      }

    }

    for(auto& k : pixel_set){
      vx = floor(k / 100000);
      vy = k - 100000*vx;
      coor = absCenter(vx, vy, empty_raster.min_x, empty_raster.min_y, empty_raster.pixel_size);
      hc.x_center = coor[0];
      hc.y_center = coor[1];
      hc.n_votes = votes[vx][vy];
      //p_circles.push_back(hc);

      for(unsigned i = 0; i < g_circles.size(); ++i){
        double dist = sqrt( pow(hc.x_center - g_circles[i].avg_x, 2) + pow(hc.y_center - g_circles[i].avg_y, 2) );

        if(dist < g_cen.aggregate_radius){
          g_circles[i].circles.push_back(hc);
          break;
        }

        if(i == (g_circles.size()-1)){
          g_cen.circles = {};
          g_cen.circles.push_back(hc);
          g_cen.avg_x = hc.x_center;
          g_cen.avg_y = hc.y_center;
          g_circles.push_back(g_cen);
        }
      }
    }

  }

  return g_circles;

}

void getPreciseCenters(vector<HoughCenters>& circles){

  circles.erase(circles.begin());
  for(unsigned i = 0; i < circles.size(); ++i){
    circles[i].getCenters();
  }

}

void assignTreeId(vector<HoughCenters>& disks, float distmax, float countDensity, unsigned minLayers){

    vector<HoughCenters>::iterator dsk;
    dsk = disks.begin();

    unsigned maxCount = 0;
    for(unsigned i = 0; i < disks.size(); i++){
      unsigned tempSize = disks[i].circles.size();
      if(maxCount < tempSize){
        maxCount = tempSize;
      }
    }
    unsigned mincount = maxCount * countDensity;

  #ifdef _DEBUG
    cout << "discs over " << mincount << " points" << endl;
  #endif // _DEBUG

    unsigned id = 1;
    while( dsk != disks.end() ){
      if(dsk->tree_id > 0 || dsk->circles.size() < mincount ){
        dsk++;
        continue;
      }

      dsk->tree_id = id++;
      double x = dsk->getMainCircle().x_center;
      double y = dsk->getMainCircle().y_center;

      vector<HoughCenters>::iterator dsk2;
      dsk2 = disks.begin();
      while( dsk2 != disks.end()){
        if(dsk2->tree_id > 0 || dsk2->circles.size() < mincount ){
          dsk2++;
          continue;
        }

        double x2 = dsk2->getMainCircle().x_center;
        double y2 = dsk2->getMainCircle().y_center;

        double euclidean = sqrt( pow(x - x2, 2) + pow(y - y2, 2) );

        if(euclidean < distmax){
          dsk2->tree_id = dsk->tree_id;
        }

        dsk2++;
      }
      dsk++;
    }

  #ifdef _DEBUG
    cout << "n of trees: " << id << endl;
  #endif // _DEBUG

    unsigned treeStacks[id] = {0};
    for(unsigned i = 0; i < disks.size(); i++){
      unsigned index = disks[i].tree_id;
      if(index > 0){
        treeStacks[index]++;
      }
    }

    for(auto& dsk : disks){
      unsigned index = dsk.tree_id;
      unsigned indexCount = treeStacks[index];

      if(indexCount < minLayers){
        dsk.tree_id = 0;
      }
    }

};

List saveCloud(vector<HoughCenters>* coordinates){

  vector<double> xout;
  vector<double> yout;
  vector<double> zout;

  vector<HoughCenters>::iterator point;
  point = coordinates->begin();
  unsigned int marker = 1;
  int maxId = 0;

  unsigned int pt = 0;
  double z;
  NumericVector xyz;
  while (point != coordinates->end()){

    if(point->tree_id == 0){
      point++;
      continue;
    }

    vector<HoughCircle>::iterator c_point;
    c_point = point->circles.begin();

    z = (point->low_z + point->up_z)/2;
    xout.push_back(point->avg_x);
    yout.push_back(point->avg_y);
    zout.push_back(z);

    // laspoint.set_z( (point->low_z + point->up_z)/2 );
    // laspoint.set_x( point->avg_x );
    // laspoint.set_y( point->avg_y );
    //
    // laspoint.set_gps_time( point->circles.size() );
    // laspoint.set_user_data( point->tree_id );
    // laspoint.set_point_source_ID(marker);
    //
    // laspoint.set_synthetic_flag(1);

    if(point->tree_id > maxId){
      maxId = point->tree_id;
    }

    // laspoint.set_synthetic_flag(0);

    int counter = 0;
    while(c_point != point->circles.end()){

      xout.push_back(c_point->x_center);
      yout.push_back(c_point->y_center);
      zout.push_back(z);

      // laspoint.set_x( (*c_point).x_center );
      // laspoint.set_y( (*c_point).y_center );
      //
      // laspoint.set_intensity( c_point->n_votes );
      //
      // if( point->main_circle == counter++ ){
      //   laspoint.set_keypoint_flag(1);
      // }else{
      //   laspoint.set_keypoint_flag(0);
      // }

      c_point++;
    }
    point++;
    marker++;
  }

  for(unsigned i = 1; i <= maxId; i++){

    unsigned counter = 0;
    float sumX = 0;
    float sumY = 0;
    float sumAvgX = 0;
    float sumAvgY = 0;

    point = coordinates->begin();
    while (point != coordinates->end()){
      if(point->tree_id == i){
        counter++;

        sumX += point->getMainCircle().x_center;
        sumY += point->getMainCircle().y_center;

        sumAvgX += point->avg_x;
        sumAvgY += point->avg_y;
      }
      point++;
    }

    if(counter > 0){

      // laspoint.set_user_data(i);

      //assign MAIN circle average point
      float mainX = sumX / counter;
      float mainY = sumY / counter;

      xout.push_back(mainX);
      yout.push_back(mainY);
      zout.push_back(z);

      // laspoint.set_x(mainX);
      // laspoint.set_y(mainY);
      //
      // laspoint.set_keypoint_flag(1);
      // laspoint.set_synthetic_flag(0);

      //assign AVERAGE coordinate from circle means
      float meanX = sumAvgX / counter;
      float meanY = sumAvgY / counter;

      xout.push_back(meanX);
      yout.push_back(meanY);
      zout.push_back(z);

      // laspoint.set_x(meanX);
      // laspoint.set_y(meanY);
      //
      // laspoint.set_keypoint_flag(0);
      // laspoint.set_synthetic_flag(1);
    }
  }

  List out;
  out["X"] = xout;
  out["Y"] = yout;
  out["Z"] = zout;
  return out;
}

/////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List singleStack(NumericMatrix& las, float pixel=0.05, float rad_max=0.25, float min_den=0.1){

    vector<HoughCenters> treeMap;
    Raster ras;

    cout << "# rasterizing cloud's slice" << endl;
    ras = getCounts(las, pixel);

    cout << "# extracting center candidates" << endl;
    treeMap = getCenters(&ras);

    cout << "# extracting center estimates" << endl;
    getPreciseCenters(treeMap);

    cout << "# assigning tree IDs" << endl;
    assignTreeId(treeMap, rad_max, min_den, 1);

    cout << "# writing cloud of center candidates" << endl;

    return saveCloud(&treeMap);

}
