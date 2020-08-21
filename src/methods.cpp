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

#include "methods.hpp"
#include <algorithm>

// convert point cloud slice to raster of point count
Raster getCounts(vector<vector<double> >& slice, double pixel_size){

  vector<double>& xcol = slice[0];
  vector<double>& ycol = slice[1];
  vector<double> minmax = getMinMax(slice);

  Raster ras;
  ras.min_x = minmax[0];
  ras.max_x = minmax[1];
  ras.min_y = minmax[2];
  ras.max_y = minmax[3];
  ras.min_z = minmax[4];
  ras.max_z = minmax[5];

  ras.pixel_size = pixel_size;

  ras.setDims();
  ras.setMatrixSize();
  ras.max_count = 0;

  for(unsigned int i = 0; i < xcol.size(); ++i){
    ras.updateMatrix(xcol[i], ycol[i]);
  }

  return ras;

}

// apply hough transform on Raster object
vector<HoughCenters> getCenters(Raster* raster, double max_radius, double min_den, unsigned int min_votes){

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

// hough transform for a single cylinder's point cloud
HoughCenters getSingleCenter(Raster* raster, double max_radius, double min_den, unsigned int min_votes){

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

// assign tree IDs to tree position point (from Raster) stacks
void assignTreeId(vector<HoughCenters>& disks, double distmax, double countDensity, unsigned int minLayers){

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

    vector<unsigned int> treeStacks(id, 0);
    for(auto& i : disks){
      unsigned index = i.tree_id;
      if(index > 0) treeStacks[index]++;
    }

    for(auto& dsk : disks){
      unsigned index = dsk.tree_id;
      unsigned indexCount = treeStacks[index];
      if(indexCount < minLayers) dsk.tree_id = 0;
    }

}

// single tree stem points detection
vector<HoughCenters> treeHough(vector<vector<double> >& cppCloud, double h1, double h2, double hstep, double radius, double pixel, double density, unsigned int votes){

  vector<double> bbox = getMinMax(cppCloud);
  vector<vector<double> > cloudSegment = getSlices(cppCloud, h1, h2, h2-h1)[0];

  if(cloudSegment[0].empty()){
    vector<HoughCenters> noPoints;
    return noPoints;
  }

  Raster raster = getCounts(cloudSegment, pixel);

  HoughCircle circle = getSingleCenter(&raster, radius, density, votes).main_circle;

  cloudSegment.clear();
  cloudSegment.shrink_to_fit();

  Raster startLayer;
  startLayer.min_x = bbox[0];
  startLayer.max_x = bbox[1];
  startLayer.min_y = bbox[2];
  startLayer.max_y = bbox[3];
  startLayer.pixel_size = pixel;
  startLayer.max_count = 0;
  startLayer.setDims();
  startLayer.setMatrixSize();

  unsigned int nLayers = ceil(bbox[5] / hstep);
  vector<Raster> treeRasters(nLayers, startLayer);

  for(unsigned int i = 0; i < cppCloud[0].size(); ++i){

    double x = cppCloud[0][i];
    double y = cppCloud[1][i];
    double z = cppCloud[2][i];

    if(z < 0) continue;

    unsigned int ptLayer = floor(z / hstep);
    Raster* alias = &treeRasters[ptLayer];

    if(alias->max_count == 0)
      alias->min_z = alias->max_z = z;
    else if(z < alias->min_x)
      alias->min_z = z;
    else if(z > alias->max_z)
      alias->max_z = z;

    alias->updateMatrix(x,y);
  }

  vector<HoughCenters> treeEstimates( treeRasters.size() );
  for(unsigned int i = 0; i < treeRasters.size(); ++i){

    if(circle.n_votes < votes) break;

    Raster* alias = &treeRasters[i];

    alias->cleanRadius(circle.x_center, circle.y_center, circle.radius + pixel*2);

    if(alias->max_count <= votes){
      if(i == 0)
        treeEstimates[i].main_circle = circle;
      else
        treeEstimates[i] = treeEstimates[i-1];
    }else
      treeEstimates[i] = getSingleCenter(alias, circle.radius + pixel*2, density , votes);

    if(alias->min_z > h2)
      circle = treeEstimates[i].main_circle;

  }

  return treeEstimates;

}

vector<vector<vector<double> > > treeEigenHough(vector<vector<double> >& cppEigenCloud, vector<unsigned int>& pointId, vector<unsigned int>& segId, double voxel_size, double max_rad, bool is2d, bool getSpace){

  IndexedCloudParts treeCloud(cppEigenCloud, pointId, segId);
  cppEigenCloud.clear();
  cppEigenCloud.shrink_to_fit();

  vector<vector<vector<double> > > results;

  for(auto& seg : treeCloud.parts){

    unsigned int id = seg.first;
    vector<vector<double> > vals = voxelCounter(seg.second.cloud, voxel_size, max_rad, is2d, getSpace);

    if(!getSpace){
      vals.push_back( treeCloud.ids2double(id) );
    };
    
    vals.push_back({ (double)id });
    results.push_back(vals);
  }

  return results;
}

// single stem methods
vector<vector<double> > ransacStemCircle(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

  vector<vector<vector<double> > > stemSlices = getChunks(cloud, segments);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<double> segRadii = idSortUnique(segments, radii);

  set<unsigned int> uniqueIds;
  for(auto& i : segments){
    uniqueIds.insert(i);
  }

  vector< vector<double> > estimates;

  for(unsigned int i = 0; i < stemSlices.size(); ++i){

    vector<vector<double> > slice = stemSlices[i];

    if(slice[0].size() <= nSamples) continue;

    double& hrad = segRadii[i];
    vector<double> temp = ransacCircle(slice, nSamples, pConfidence, pInliers, 20);

    double rdiff = abs(temp[2] - hrad);
    if(rdiff > tolerance){

      vector<double> bbox = getMinMax(slice);

      temp[0] = (bbox[1] + bbox[0])/2;
      temp[1] = (bbox[2] + bbox[3])/2;
      temp[2] = hrad;
      temp[3] = 0;
    }

    unsigned int id = *next(uniqueIds.begin(), i);
    temp.push_back(id);
    estimates.push_back(temp);
  }

  return estimates;

}

vector<vector<double> > ransacStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

  vector<vector<vector<double> > > stemSlices = getChunks(cloud, segments);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<double> segRadii = idSortUnique(segments, radii);

  set<unsigned int> uniqueIds;
  for(auto& i : segments){
    uniqueIds.insert(i);
  }

  vector< vector<double> > estimates;

  for(unsigned int i = 0; i < stemSlices.size(); ++i){

    vector<vector<double> > slice = stemSlices[i];

    // cout << "... seg " << i+1 << " of " << stemSlices.size() << ", n points: " << slice[0].size() << endl;

    if(slice[0].size() <= nSamples) continue;

    double& hrad = segRadii[i];
    vector<double> temp = ransacCylinder(slice, nSamples, pConfidence, pInliers);

    double rdiff = abs(temp[4] - hrad);
    if(rdiff > tolerance){
      temp[0] = 0;
      temp[1] = PI/2;
      temp[2] = 0;
      temp[3] = 0;
      temp[4] = hrad;
      temp[5] = 0;
    }

    unsigned int id = *next(uniqueIds.begin(), i);
    temp.push_back(id);
    estimates.push_back(temp);
  }

  return estimates;

}

vector<vector<double> > irlsStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints,  double tolerance){

  vector<vector<vector<double> > > stemSlices = getChunks(cloud, segments);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<double> segRadii = idSortUnique(segments, radii);

  set<unsigned int> uniqueIds;
  for(auto& i : segments){
    uniqueIds.insert(i);
  }

  vector< vector<double> > estimates;
  vector<double> initPars = {0, PI/2, 0, 0, 0};

  for(unsigned int i = 0; i < stemSlices.size(); ++i){

    vector<vector<double> > slice = stemSlices[i];

    // cout << "... seg " << i+1 << " of " << stemSlices.size() << ", n points: " << slice[0].size() << endl;

    if(slice[0].size() <= 5) continue;

    if(slice[0].size() > nPoints && nPoints > 0){
      double p = (double)nPoints / (double)slice[0].size();
      slice = randomPoints(slice, p);
    }

    double& hrad = segRadii[i];
    vector<double> temp = irlsCylinder(slice, initPars);

    double rdiff = abs(temp[4] - hrad);
    if(rdiff > tolerance){
      temp[0] = 0;
      temp[1] = PI/2;
      temp[2] = 0;
      temp[3] = 0;
      temp[4] = hrad;
      temp[5] = 0;
    }

    unsigned int id = *next(uniqueIds.begin(), i);
    temp.push_back(id);
    estimates.push_back(temp);
  }

  return estimates;

}

vector<vector<double> > irlsStemCircle(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints, double tolerance){

  vector<vector<vector<double> > > stemSlices = getChunks(cloud, segments);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<double> segRadii = idSortUnique(segments, radii);

  set<unsigned int> uniqueIds;
  for(auto& i : segments){
    uniqueIds.insert(i);
  }

  vector< vector<double> > estimates;

  for(unsigned int i = 0; i < stemSlices.size(); ++i){
    vector<vector<double> >& slice = stemSlices[i];

    if(slice[0].size() <= 3) continue;

    if(slice[0].size() > nPoints && nPoints > 0){
      double p = (double)nPoints / (double)slice[0].size();
      slice = randomPoints(slice, p);
    }

    double& hrad = segRadii[i];
    vector<double> initPars = eigenCircle(slice);
    vector<double> temp = irlsCircle(slice, initPars);

    double rdiff = abs(temp[2] - hrad);
    if(rdiff > tolerance){
      temp[0] = 0;
      temp[1] = 0;
      temp[2] = hrad;
      temp[3] = 0;
    }

    unsigned int id = *next(uniqueIds.begin(), i);
    temp.push_back(id);
    estimates.push_back(temp);
  }

  return estimates;

}

// plot-wise methods
vector<vector<vector<double> > > ransacPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;

  unsigned int progress_counter = 0;
  unsigned int n_ids = uniqueTotalCounter(treeId);
  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;
    progressPrinter("trees", progress_counter++, n_ids);

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = ransacStemCircle(tree, segs, segsRadii, nSamples, pConfidence, pInliers, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}

vector<vector<vector<double> > > ransacPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);
  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;
  
  unsigned int progress_counter = 0;
  unsigned int n_ids = uniqueTotalCounter(treeId);
  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;
    progressPrinter("trees", progress_counter++, n_ids);

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = ransacStemCylinder(tree, segs, segsRadii, nSamples, pConfidence, pInliers, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}

vector<vector<vector<double> > > irlsPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);
  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;

  unsigned int progress_counter = 0;
  unsigned int n_ids = uniqueTotalCounter(treeId);
  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;
    progressPrinter("trees", progress_counter++, n_ids);

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = irlsStemCylinder(tree, segs, segsRadii, nPoints, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}

vector<vector<vector<double> > > irlsPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nPoints, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);
  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;

  unsigned int progress_counter = 0;
  unsigned int n_ids = uniqueTotalCounter(treeId);
  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;
    progressPrinter("trees", progress_counter++, n_ids);

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = irlsStemCircle(tree, segs, segsRadii, nPoints, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}

vector<vector<double> > pointMetrics(vector<vector<double> >& cloud, vector<vector<unsigned int> >& idx, vector<bool> which_metrics){

  unsigned int ncol = idx.size();
  unsigned int nrow = idx[0].size();

  vector<vector<double> > out;
  for(unsigned int i = 0; i < nrow; ++i){

    vector<vector<double> > xyz(3);

    for(unsigned int j = 0; j < ncol; ++j){

      int cell = idx[j][i];
      if(cell-- == 0) break;

      xyz[0].push_back( cloud[0][cell] );
      xyz[1].push_back( cloud[1][cell] );
      xyz[2].push_back( cloud[2][cell] );
    }

    if(xyz[0].size() < 3){
      out.push_back({0});
    }else{
      vector<double> metrics = nnMetrics(xyz, which_metrics);
      out.push_back(metrics);
    }
  }

  return out;
}

vector<vector<double> > voxelMetrics(vector<vector<double> >& cloud, vector<vector<unsigned int> >& idx, vector<bool> which_metrics){

  unsigned int nvoxels = idx.size();
  vector<vector<double> > out;

  int i = 1;
  for(auto& vx : idx){

    if(vx.size() < 3){
      out.push_back({0});
      continue;
    }

    vector<vector<double> > xyz(3);

    for(auto& pt : vx){
      xyz[0].push_back( cloud[0][pt] );
      xyz[1].push_back( cloud[1][pt] );
      xyz[2].push_back( cloud[2][pt] );
    }

    vector<double> metrics = nnMetrics(xyz, which_metrics);
    out.push_back(metrics);
  }

  return out;
}

vector<unsigned long long int> voxelIndex(vector<vector<double> >& cloud, double voxel_spacing){

  typedef unsigned long long int llint;

  double xoffset = *min_element(cloud[0].begin(), cloud[0].end());
  double yoffset = *min_element(cloud[1].begin(), cloud[1].end());
  double zoffset = *min_element(cloud[2].begin(), cloud[2].end());

  VoxelGrid voxelRegistry(xoffset, yoffset, zoffset, voxel_spacing);
  vector<llint> indexer(cloud[0].size());

  for(unsigned int i = 0; i < indexer.size(); ++i){
    indexer[i] = voxelRegistry.getVoxelHash(cloud[0][i], cloud[1][i], cloud[2][i]);
  }

  llint mindex = *min_element(indexer.begin(), indexer.end());
  for(auto& i : indexer) i -= mindex;

  return indexer;
}

vector<vector<vector<double> > > plotEigenHough(vector<vector<double> >& cppEigenCloud, vector<unsigned int>& pointId, vector<unsigned int>& treeId, vector<unsigned int>& segId, double voxel_size, double max_rad, bool is2d, bool getSpace){

  IndexedCloudParts plotCloud(cppEigenCloud, pointId, treeId, segId);
  cppEigenCloud.clear();
  cppEigenCloud.shrink_to_fit();
  
  vector<vector<vector<double> > > plotResults;

  for(auto& tree : plotCloud.parts){

    unsigned int id = tree.first;

    vector<vector<vector<double> > > calcTree = treeEigenHough(tree.second.cloud, tree.second.uniqueIds, tree.second.indexer, voxel_size, max_rad, is2d, getSpace);

    for(auto& temp : calcTree){
      temp.push_back({ (double)id });
    }

    plotResults.insert(plotResults.end(), calcTree.begin(), calcTree.end());
  }

  return plotResults;
}

vector<unsigned int> treeIdsFromMap(vector<vector<double> >& xy, vector<vector<double> >& xymap, vector<unsigned int> ids, double length, bool circle){

  vector<unsigned int> treeIds( xy[0].size() );

  for(unsigned int i = 0; i < xy[0].size(); ++i){
    double x = xy[0][i];
    double y = xy[1][i];

    for(unsigned int j = 0; j < ids.size(); ++j){
      double xref = xymap[0][j];
      double yref = xymap[1][j];

      bool isInside = false;
      if(circle){
        double dist = sqrt( pow(x - xref, 2) + pow(y - yref, 2) );
        isInside = dist < length;
      }else{
        double xdist = abs(x - xref);
        double ydist = abs(y - yref);
        isInside = xdist < (length/2) && ydist < (length/2);
      }

      if(isInside){
        treeIds[i] = ids[j];
        break;
      }
    }
  }

  return treeIds;
}

vector<vector<double> > bfStemCylinder(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double max_angle, double tolerance){

  vector<vector<vector<double> > > stemSlices = getChunks(cloud, segments);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<double> segRadii = idSortUnique(segments, radii);

  set<unsigned int> uniqueIds;
  for(auto& i : segments){
    uniqueIds.insert(i);
  }

  vector< vector<double> > estimates;

  for(unsigned int i = 0; i < stemSlices.size(); ++i){

    vector<vector<double> > slice = stemSlices[i];

    // cout << "... seg " << i+1 << " of " << stemSlices.size() << ", n points: " << slice[0].size() << endl;

    if(slice[0].size() <= nSamples) continue;

    double& hrad = segRadii[i];
    vector<double> temp = bruteForceRansacCylinder(slice, nSamples, pConfidence, pInliers, 20, max_angle, true)[0];

    double rdiff = abs(temp[2] - hrad);
    if(rdiff > tolerance){
      temp[0] = 0;
      temp[1] = 0;
      temp[2] = hrad;
      temp[3] = 0;
      temp[4] = 0;
      temp[5] = 0;
    }

    unsigned int id = *next(uniqueIds.begin(), i);
    temp.push_back(id);
    estimates.push_back(temp);
  }

  return estimates;

}

vector<vector<vector<double> > > bfPlotCylinders(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double max_angle, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);
  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;

  unsigned int progress_counter = 0;
  unsigned int n_ids = uniqueTotalCounter(treeId);
  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;
    progressPrinter("trees", progress_counter++, n_ids);

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = bfStemCylinder(tree, segs, segsRadii, nSamples, pConfidence, pInliers, max_angle, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}
