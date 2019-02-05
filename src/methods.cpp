#include <iostream>
#include "classes.hpp"

using namespace Rcpp;
using namespace std;

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

  ras.setDims();
  ras.setMatrixSize();

  unsigned int max_votes = 0;
  float x,y;
  unsigned int xCell;
  unsigned int yCell;

  for(int i = 0; i < slice.nrow(); i++){

    x = xcol[i];
    y = ycol[i];

    double dist = sqrt( pow(x - x_mid,2) + pow(y - y_mid,2) );
    if(d_mid > 0 && dist > d_mid) continue;

    xCell = abs( floor( (x - ras.min_x) / pixel_size ) );
    yCell = abs( floor( (y - ras.min_y) / pixel_size ) );

    if(xCell >= ras.x_dim) xCell = ras.x_dim-1;
    if(yCell >= ras.y_dim) yCell = ras.y_dim-1;

    ras.matrix[xCell][yCell]++;

    if(ras.matrix[xCell][yCell] > max_votes) max_votes = ras.matrix[xCell][yCell];
  }

  ras.max_count = max_votes;

  return ras;

}

vector<HoughCenters> getCenters(Raster* raster, float max_radius=0.25, float min_den=0.1, unsigned int min_votes=3){

  //count raster properties (&raster)
  unsigned int& x_len = raster->x_dim;
  unsigned int& y_len = raster->y_dim;
  unsigned int min_count = ceil( (raster->max_count)*min_den );

  //empty raster properties (votes)
  Raster empty_raster;
  empty_raster.min_x = raster->min_x - max_radius;
  empty_raster.min_y = raster->min_y - max_radius;
  empty_raster.max_x = raster->max_x + max_radius;
  empty_raster.max_y = raster->max_y + max_radius;
  empty_raster.pixel_size = raster->pixel_size;
  empty_raster.setDims();
  empty_raster.setMatrixSize();

  //get valid pixels
  vector< array<unsigned int,2> > pixels;

  for(unsigned i = 0; i < x_len; ++i){
    for(unsigned j = 0; j < y_len; ++j){
      if(raster->matrix[i][j] >= min_count )
        pixels.push_back( {i,j} );
    }
  }

  //make circles centered in every valid pixel
  vector<HoughCenters> g_circles;
  float inclusionRadius = max_radius * 3;

  for(float rad = raster->pixel_size; rad <= max_radius; rad+=raster->pixel_size){

    vector< vector<unsigned int> > votes = empty_raster.matrix;
    PixelSet pixel_set;

    for(auto& p : pixels){

      vector<float> center = raster->absCenter(p[0], p[1]);
      PixelSet h_circle = empty_raster.rasterCircle(rad, center[0], center[1]);

      for(auto& h : h_circle){
        unsigned int vx = h[0];
        unsigned int vy = h[1];

        if( vx >= votes.size() || vy >= votes[0].size() ) continue;

        votes[ vx ][ vy ]++;

        if(votes[ vx ][ vy ] >= min_votes){
          array<unsigned int, 2> vxy = {vx, vy};
          pixel_set.insert( vxy );
        }
      }
    }

    for(auto& k : pixel_set){
      unsigned int vx = k[0];
      unsigned int vy = k[1];
      vector<float> coor = empty_raster.absCenter(vx, vy);

      HoughCircle hc;
      hc.radius = rad;
      hc.x_center = coor[0];
      hc.y_center = coor[1];
      hc.n_votes = votes[vx][vy];

      // disposable when looking at single circle
      bool newCircle = true;
      for(auto& gc : g_circles){

        float dist = sqrt( pow(hc.x_center - gc.avg_x, 2) + pow(hc.y_center - gc.avg_y, 2) );

        if(dist < inclusionRadius){
          gc.circles.push_back(hc);
          newCircle = false;
          break;
        }
      }

      if(newCircle){
        HoughCenters g_cen;
        g_cen.circles = {hc};
        g_cen.avg_x = hc.x_center;
        g_cen.avg_y = hc.y_center;
        g_cen.low_z = raster->min_z;
        g_cen.up_z = raster->max_z;
        g_cen.aggregate_radius = inclusionRadius;
        g_circles.push_back(g_cen);
      }
    }
  }

  for(auto& i : g_circles) i.getCenters();

  return g_circles;

}

HoughCenters getSingleCenter(Raster* raster, float max_radius=0.25, float min_den=0.1, unsigned int min_votes=3){

  //count raster properties (&raster)
  unsigned int& x_len = raster->x_dim;
  unsigned int& y_len = raster->y_dim;
  unsigned int min_count = ceil( (raster->max_count)*min_den );

  //empty raster properties (votes)
  Raster empty_raster;
  empty_raster.min_x = raster->min_x - max_radius;
  empty_raster.min_y = raster->min_y - max_radius;
  empty_raster.max_x = raster->max_x + max_radius;
  empty_raster.max_y = raster->max_y + max_radius;
  empty_raster.pixel_size = raster->pixel_size;
  empty_raster.setDims();
  empty_raster.setMatrixSize();

  //get valid pixels
  vector< array<unsigned int,2> > pixels;

  for(unsigned i = 0; i < x_len; ++i){
    for(unsigned j = 0; j < y_len; ++j){
      if(raster->matrix[i][j] >= min_count )
        pixels.push_back( {i,j} );
    }
  }

  //make circles centered in every valid pixel
  HoughCenters circle_candidates;
  circle_candidates.low_z = raster->min_z;
  circle_candidates.up_z = raster->max_z;

  for(float rad = raster->pixel_size; rad <= max_radius; rad+=raster->pixel_size){

    vector< vector<unsigned int> > votes = empty_raster.matrix;
    PixelSet pixel_set;

    for(auto& p : pixels){

      vector<float> center = raster->absCenter(p[0], p[1]);
      PixelSet h_circle = empty_raster.rasterCircle(rad, center[0], center[1]);

      for(auto& h : h_circle){
        unsigned int vx = h[0];
        unsigned int vy = h[1];

        if( vx >= votes.size() || vy >= votes[0].size() ) continue;

        votes[ vx ][ vy ]++;

        if(votes[ vx ][ vy ] >= min_votes){
          array<unsigned int, 2> vxy = {vx, vy};
          pixel_set.insert( vxy );
        }
      }
    }

    for(auto& k : pixel_set){
      unsigned int vx = k[0];
      unsigned int vy = k[1];
      vector<float> coor = empty_raster.absCenter(vx, vy);

      HoughCircle hc;
      hc.radius = rad;
      hc.x_center = coor[0];
      hc.y_center = coor[1];
      hc.n_votes = votes[vx][vy];

      circle_candidates.circles.push_back(hc);
    }
  }

  circle_candidates.getCenters();

  return circle_candidates;

}

void assignTreeId(vector<HoughCenters>& disks, float distmax, float countDensity, unsigned minLayers){

    unsigned maxCount = 0;
    for(auto& i : disks){
      unsigned tempSize = i.circles.size();
      if(maxCount < tempSize){
        maxCount = tempSize;
      }
    }

    unsigned mincount = maxCount * countDensity;

    vector<HoughCenters>::iterator dsk;
    dsk = disks.begin();

    unsigned id = 1;
    while(dsk != disks.end()){

      if(dsk->tree_id > 0 || dsk->circles.size() < mincount ){
        dsk++;
        continue;
      }

      dsk->tree_id = id++;
      float x = dsk->main_circle->x_center;
      float y = dsk->main_circle->y_center;

      vector<HoughCenters>::iterator dsk2;
      dsk2 = disks.begin();
      while(dsk2 != disks.end()){

        if(dsk2->tree_id > 0 || dsk2->circles.size() < mincount ){
          dsk2++;
          continue;
        }

        float x2 = dsk2->main_circle->x_center;
        float y2 = dsk2->main_circle->y_center;

        float euclidean = sqrt( pow(x - x2, 2) + pow(y - y2, 2) );

        if(euclidean < distmax){
          dsk2->tree_id = dsk->tree_id;
        }

        dsk2++;
      }
      dsk++;
    }

    unsigned treeStacks[id] = {0};
    for(auto& i : disks){
      unsigned index = i.tree_id;
      if(index > 0) treeStacks[index]++;
    }

    for(auto& dsk : disks){
      unsigned index = dsk.tree_id;
      unsigned indexCount = treeStacks[index];
      if(indexCount < minLayers) dsk.tree_id = 0;
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

        sumX += point->main_circle->x_center;
        sumY += point->main_circle->y_center;

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
List singleStack(NumericMatrix& las, float pixel=0.05, float rad_max=0.25, float min_den=0.1, unsigned int min_votes = 3){

    vector<HoughCenters> treeMap;
    Raster ras;

    cout << "# rasterizing cloud's slice" << endl;
    ras = getCounts(las, pixel);

    cout << "# extracting center candidates" << endl;
    treeMap = getCenters(&ras, rad_max, min_den, min_votes);

    cout << "# assigning tree IDs" << endl;
    assignTreeId(treeMap, rad_max, min_den, 1);

    cout << "# writing cloud of center candidates" << endl;

    return saveCloud(&treeMap);

}

// [[Rcpp::export]]
List getCircle(NumericMatrix& las, float pixel=0.05, float rad_max=0.25, float min_den=0.1, unsigned int min_votes = 3){

  Raster ras = getCounts(las, pixel);
  HoughCenters circle = getSingleCenter(&ras, rad_max, min_den, min_votes);

  List out;
  out["x"] = circle.main_circle->x_center;
  out["y"] = circle.main_circle->y_center;
  out["rad"] = circle.main_circle->radius;
  out["votes"] = circle.main_circle->n_votes;

  return out;
}
