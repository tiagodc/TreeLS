#include <iostream>
#include "classes.hpp"

using namespace Rcpp;
using namespace std;

vector<vector<double*> > rmatrix2cpp(NumericMatrix& cloud){

  NumericMatrix::Column xcol = cloud( _, 0);
  NumericMatrix::Column ycol = cloud( _, 1);
  NumericMatrix::Column zcol = cloud( _, 2);

  vector<vector<double*> > xyz(3);

  for(int i = 0; i < cloud.nrow(); ++i){
    xyz[0].push_back(&xcol[i]);
    xyz[1].push_back(&ycol[i]);
    xyz[2].push_back(&zcol[i]);
  }

  return xyz;
}

vector<double> getMinMax(vector<vector<double*> >& xyz){

  vector<double> minmax(6);
  minmax[0] = minmax[1] = *xyz[0][0];
  minmax[2] = minmax[3] = *xyz[1][0];
  minmax[4] = minmax[5] = *xyz[2][0];

  for(unsigned int i = 1; i < xyz[0].size(); ++i){

    double x = *xyz[0][i];
    double y = *xyz[1][i];
    double z = *xyz[2][i];

    if(x < minmax[0]) minmax[0] = x; else if(x > minmax[1]) minmax[1] = x;
    if(y < minmax[2]) minmax[2] = y; else if(y > minmax[3]) minmax[3] = y;
    if(z < minmax[4]) minmax[4] = z; else if(z > minmax[5]) minmax[5] = z;
  }

  return minmax;

}

vector<vector<vector<double*> > > getSlices(NumericMatrix& cloud, double zmin = 1, double zmax=3, double zstep = 0.5){

  NumericMatrix::Column xcol = cloud( _, 0);
  NumericMatrix::Column ycol = cloud( _, 1);
  NumericMatrix::Column zcol = cloud( _, 2);

  unsigned int nlayers = ceil((zmax - zmin) / zstep);
  vector<vector<vector<double*> > > store(nlayers, vector<vector<double*> >(3));

  for(int i = 0; i < cloud.nrow(); ++i){

    double z = zcol[i];
    if(z < zmin || z >= zmax)
      continue;

    unsigned int n = floor((z - zmin) / zstep);

    store[n][0].push_back(&xcol[i]);
    store[n][1].push_back(&ycol[i]);
    store[n][2].push_back(&zcol[i]);
  }

  return store;

}

Raster getCounts(vector<vector<double*> >& slice , double pixel_size){

  vector<double*>& xcol = slice[0];
  vector<double*>& ycol = slice[1];
  vector<double> minmax = getMinMax(slice);

  Raster ras;
  ras.min_x = minmax[0];
  ras.max_x = minmax[1];
  ras.min_y = minmax[2];
  ras.max_y = minmax[3];
  ras.min_z = minmax[4];
  ras.max_z = minmax[5];

  ras.pixel_size = pixel_size;
  ras.thickness = ras.max_z - ras.min_z;

  ras.setDims();
  ras.setMatrixSize();
  ras.max_count = 0;

  double x,y;
  unsigned int xCell;
  unsigned int yCell;

  for(unsigned int i = 0; i < xcol.size(); ++i){

    x = *xcol[i];
    y = *ycol[i];

    xCell = abs( floor( (x - ras.min_x) / pixel_size ) );
    yCell = abs( floor( (y - ras.min_y) / pixel_size ) );

    if(xCell >= ras.x_dim) xCell = ras.x_dim-1;
    if(yCell >= ras.y_dim) yCell = ras.y_dim-1;

    ras.matrix[xCell][yCell]++;

    if(ras.matrix[xCell][yCell] > ras.max_count) ras.max_count = ras.matrix[xCell][yCell];
  }

  return ras;

}

vector<HoughCenters> getCenters(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3){

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
  double inclusionRadius = max_radius * 3;

  for(double rad = raster->pixel_size; rad <= max_radius; rad+=raster->pixel_size){

    vector< vector<unsigned int> > votes = empty_raster.matrix;
    PixelSet pixel_set;

    for(auto& p : pixels){

      vector<double> center = raster->absCenter(p[0], p[1]);
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
      vector<double> coor = empty_raster.absCenter(vx, vy);

      HoughCircle hc;
      hc.radius = rad;
      hc.x_center = coor[0];
      hc.y_center = coor[1];
      hc.n_votes = votes[vx][vy];

      // disposable when looking at single circle
      bool newCircle = true;
      for(auto& gc : g_circles){

        double dist = sqrt( pow(hc.x_center - gc.avg_x, 2) + pow(hc.y_center - gc.avg_y, 2) );

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

HoughCenters getSingleCenter(Raster* raster, double max_radius=0.25, double min_den=0.1, unsigned int min_votes=3){

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

  for(double rad = raster->pixel_size; rad <= max_radius; rad+=raster->pixel_size){

    vector< vector<unsigned int> > votes = empty_raster.matrix;
    PixelSet pixel_set;

    for(auto& p : pixels){

      vector<double> center = raster->absCenter(p[0], p[1]);
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
      vector<double> coor = empty_raster.absCenter(vx, vy);

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

void assignTreeId(vector<HoughCenters>& disks, double distmax, double countDensity, unsigned int minLayers=1){

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
      double x = dsk->main_circle.x_center;
      double y = dsk->main_circle.y_center;

      vector<HoughCenters>::iterator dsk2;
      dsk2 = disks.begin();

      while(dsk2 != disks.end()){

        if(dsk2->tree_id > 0 || dsk2->circles.size() < mincount ){
          dsk2++;
          continue;
        }

        double x2 = dsk2->main_circle.x_center;
        double y2 = dsk2->main_circle.y_center;

        double euclidean = sqrt( pow(x - x2, 2) + pow(y - y2, 2) );

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

List exportTreeMap(vector<HoughCenters>& coordinates){

  vector<double> xout;
  vector<double> yout;
  vector<double> zout;
  vector<double> radii;
  vector<bool> keyFlag;
  vector<bool> treeFlag;
  vector<unsigned short int> votes;
  vector<unsigned int> treeId;
  vector<unsigned int> discId;
  vector<unsigned int> nPoints;


  unsigned int diskCounter = 1;
  unsigned int maxId = 0;

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    if(point.tree_id > maxId)
      maxId = point.tree_id;

    double z = (point.low_z + point.up_z)/2;

    // main point
    xout.push_back(point.main_circle.x_center);
    yout.push_back(point.main_circle.y_center);
    zout.push_back(z);

    radii.push_back(point.main_circle.radius);
    votes.push_back(point.main_circle.n_votes);
    treeId.push_back(point.tree_id);
    nPoints.push_back(point.circles.size());
    discId.push_back(diskCounter);
    keyFlag.push_back(true);
    treeFlag.push_back(false);

    // other candidates
    for(auto& c_point : point.circles){

      xout.push_back(c_point.x_center);
      yout.push_back(c_point.y_center);
      zout.push_back(z);

      radii.push_back(c_point.radius);
      votes.push_back(c_point.n_votes);
      treeId.push_back(point.tree_id);
      nPoints.push_back(point.circles.size());
      discId.push_back(diskCounter);
      keyFlag.push_back(false);
      treeFlag.push_back(false);
    }
    diskCounter++;
  }

  double xSums[maxId] = {0};
  double ySums[maxId] = {0};
  unsigned int counters[maxId] = {0};

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    xSums[point.tree_id-1] += point.main_circle.x_center;
    ySums[point.tree_id-1] += point.main_circle.y_center;
    counters[point.tree_id-1]++;
  }

  for(unsigned int i = 0; i < maxId; ++i){

    if(counters[i] == 0)
      continue;

    double mainX = xSums[i] / counters[i];
    double mainY = ySums[i] / counters[i];

    xout.push_back(mainX);
    yout.push_back(mainY);
    zout.push_back(0);

    radii.push_back(0);
    votes.push_back(0);
    treeId.push_back(i+1);
    nPoints.push_back(0);
    discId.push_back(0);
    keyFlag.push_back(true);
    treeFlag.push_back(true);
  }

  List out;
  out["X"] = xout;
  out["Y"] = yout;
  out["Z"] = zout;
  out["Intensity"] = votes;
  out["PointSourceID"] = discId;
  out["Keypoint_flag"] = keyFlag;
  out["Radii"] = radii;
  out["TreeID"] = treeId;
  out["TreePosition"] = treeFlag;
  out["n"] = nPoints;

  return out;
}

/////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List getCircle(NumericMatrix& las, double pixel=0.05, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

  vector<vector<double*> > cloud = rmatrix2cpp(las);
  Raster ras = getCounts(cloud, pixel);
  HoughCenters circle = getSingleCenter(&ras, rad_max, min_den, min_votes);

  List out;
  out["x"] = circle.main_circle.x_center;
  out["y"] = circle.main_circle.y_center;
  out["rad"] = circle.main_circle.radius;
  out["votes"] = circle.main_circle.n_votes;

  return out;
}

// [[Rcpp::export]]
List singleStack(NumericMatrix& las, double pixel=0.05, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

    vector<HoughCenters> treeMap;
    Raster ras;

    vector<vector<double*> > stack = rmatrix2cpp(las);
    ras = getCounts(stack, pixel);
    treeMap = getCenters(&ras, rad_max, min_den, min_votes);
    assignTreeId(treeMap, rad_max, min_den);

    return exportTreeMap(treeMap);

}

// [[Rcpp::export]]
List stackMap(NumericMatrix& las, double hmin=1, double hmax=3, double hstep=0.5, double pixel=0.025, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

  vector<HoughCenters> treeMap;

  vector<vector<vector<double*> > > fullStack = getSlices(las, hmin, hmax, hstep);

  for(auto& stack : fullStack){

    if(stack[0].empty())
      continue;

    Raster ras = getCounts(stack, pixel);
    vector<HoughCenters> tempMap = getCenters(&ras, rad_max, min_den, min_votes);
    treeMap.insert(treeMap.end(), tempMap.begin(), tempMap.end());
  }

  unsigned int nlayers = 0.75 * (hmax - hmin) / hstep;
  assignTreeId(treeMap, rad_max, min_den, nlayers);

  return exportTreeMap(treeMap);

}
