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

#include "methods.hpp"

// convert r matrix to std::vector
vector<vector<double> > rmatrix2cpp(NumericMatrix& cloud){

  NumericMatrix::Column xcol = cloud( _, 0);
  NumericMatrix::Column ycol = cloud( _, 1);
  NumericMatrix::Column zcol = cloud( _, 2);

  vector<vector<double> > xyz(3);

  xyz[0].insert(xyz[0].begin(), xcol.begin(), xcol.end());
  xyz[1].insert(xyz[1].begin(), ycol.begin(), ycol.end());
  xyz[2].insert(xyz[2].begin(), zcol.begin(), zcol.end());

  return xyz;
}

// get min values from XYZ point cloud
vector<double> getMinMax(vector<vector<double> >& xyz){

  vector<double> minmax(6);
  minmax[0] = minmax[1] = xyz[0][0];
  minmax[2] = minmax[3] = xyz[1][0];
  minmax[4] = minmax[5] = xyz[2][0];

  for(unsigned int i = 1; i < xyz[0].size(); ++i){

    double x = xyz[0][i];
    double y = xyz[1][i];
    double z = xyz[2][i];

    if(x < minmax[0]) minmax[0] = x; else if(x > minmax[1]) minmax[1] = x;
    if(y < minmax[2]) minmax[2] = y; else if(y > minmax[3]) minmax[3] = y;
    if(z < minmax[4]) minmax[4] = z; else if(z > minmax[5]) minmax[5] = z;
  }

  return minmax;

}

// crop point cloud
vector<vector<double> > cropCloud(vector<vector<double> > cloud, double xCenter, double yCenter, double len, bool circle, bool negative){

  vector<vector<double> > keepCloud(3);

  for(unsigned int i = 0; i < cloud[0].size(); ++i){

    double& x = cloud[0][i];
    double& y = cloud[1][i];
    double& z = cloud[2][i];

    bool keep;
    if(circle){
      double dist = sqrt( pow(x - xCenter, 2) + pow(y - yCenter, 2) );
      keep = dist < len;
    }else{
      double dx = abs(x - xCenter);
      double dy = abs(y - yCenter);
      keep = dx < (len/2) && dy < (len/2);
    }

    if(negative) keep = !keep;

    if(keep){
      keepCloud[0].push_back(x);
      keepCloud[1].push_back(y);
      keepCloud[2].push_back(z);
    }

  }

  return keepCloud;

}

vector<bool> cropCloudFilter(vector<vector<double> > cloud, double xCenter, double yCenter, double len, bool circle, bool negative){

  vector<bool> keepCloud( cloud[0].size() );

  for(unsigned int i = 0; i < cloud[0].size(); ++i){

    double& x = cloud[0][i];
    double& y = cloud[1][i];

    bool keep;
    if(circle){
      double dist = sqrt( pow(x - xCenter, 2) + pow(y - yCenter, 2) );
      keep = dist < len;
    }else{
      double dx = abs(x - xCenter);
      double dy = abs(y - yCenter);
      keep = dx < (len/2) && dy < (len/2);
    }

   keepCloud[i] = negative ? !keep : keep;

  }

  cloud.clear();
  cloud.shrink_to_fit();

  return keepCloud;

}

// systematic voxel grid point filter (-thin_with_voxel)
vector<bool> voxelFilter(vector<vector<double> >& cloud, double voxel_spacing){

  double& xoffset = cloud[0][0];
  double& yoffset = cloud[1][0];
  double& zoffset = cloud[2][0];

  VoxelSet ledger;
  vector<bool> filter(cloud[0].size());

  for(unsigned int i = 0; i < cloud[0].size(); ++i){

   int nx = floor( (cloud[0][i] - xoffset) / voxel_spacing);
   int ny = floor( (cloud[1][i] - yoffset) / voxel_spacing);
   int nz = floor( (cloud[2][i] - zoffset) / voxel_spacing);

   array<int,3> voxel = {nx, ny, nz};
   filter[i] = ledger.insert(voxel).second;
  }

  cloud.clear();
  cloud.shrink_to_fit();

  return filter;

}

// split point cloud into horizontal slices
vector<vector<vector<double> > > getSlices(NumericMatrix& cloud, double zmin, double zmax, double zstep){

  NumericMatrix::Column xcol = cloud( _, 0);
  NumericMatrix::Column ycol = cloud( _, 1);
  NumericMatrix::Column zcol = cloud( _, 2);

  unsigned int nlayers = ceil((zmax - zmin) / zstep);
  vector<vector<vector<double> > > store(nlayers, vector<vector<double> >(3));

  for(int i = 0; i < cloud.nrow(); ++i){

    double z = zcol[i];
    if(z < zmin || z >= zmax)
      continue;

    unsigned int n = floor((z - zmin) / zstep);

    store[n][0].push_back(xcol[i]);
    store[n][1].push_back(ycol[i]);
    store[n][2].push_back(zcol[i]);
  }

  return store;

}

vector<vector<vector<double> > > getSlices(vector<vector<double> >& cloud, double zmin, double zmax, double zstep){

  vector<double>& xcol = cloud[0];
  vector<double>& ycol = cloud[1];
  vector<double>& zcol = cloud[2];

  unsigned int nlayers = ceil((zmax - zmin) / zstep);
  vector<vector<vector<double> > > store(nlayers, vector<vector<double> >(3));

  for(unsigned int i = 0; i < xcol.size(); ++i){

    double z = zcol[i];
    if(z < zmin || z >= zmax)
      continue;

    unsigned int n = floor((z - zmin) / zstep);

    store[n][0].push_back(xcol[i]);
    store[n][1].push_back(ycol[i]);
    store[n][2].push_back(zcol[i]);
  }

  return store;

}

// split point cloud according to some criterion
vector<vector<vector<double> > > getChunks(vector<vector<double> >& cloud, vector<unsigned int>& identifier){

  vector<double>& xcol = cloud[0];
  vector<double>& ycol = cloud[1];
  vector<double>& zcol = cloud[2];

  unsigned int minIndex = *min_element(identifier.begin(), identifier.end());
  unsigned int maxIndex = *max_element(identifier.begin(), identifier.end());
  unsigned int nlayers =  maxIndex - minIndex + 1;

  vector<vector<vector<double> > > store(nlayers, vector<vector<double> >(3));

  for(unsigned int i = 0; i < xcol.size(); ++i){

    unsigned int n = identifier[i] - minIndex;

    store[n][0].push_back(xcol[i]);
    store[n][1].push_back(ycol[i]);
    store[n][2].push_back(zcol[i]);
  }

  return store;

}

vector<vector<unsigned int> > partitionIndex(vector<unsigned int>& identifier, vector<unsigned int>& partitioner){

  unsigned int minIndex = *min_element(identifier.begin(), identifier.end());
  unsigned int maxIndex = *max_element(identifier.begin(), identifier.end());
  unsigned int nlayers =  maxIndex - minIndex + 1;

  vector< vector<unsigned int> > store(nlayers);

  for(unsigned int i = 0; i < identifier.size(); ++i){
    unsigned int n = identifier[i] - minIndex;
    store[n].push_back(partitioner[i]);
  }

  return store;

}

vector<vector<double> > partitionIndex(vector<unsigned int>& identifier, vector<double>& partitioner){

  unsigned int minIndex = *min_element(identifier.begin(), identifier.end());
  unsigned int maxIndex = *max_element(identifier.begin(), identifier.end());
  unsigned int nlayers =  maxIndex - minIndex + 1;

  vector< vector<double> > store(nlayers);

  for(unsigned int i = 0; i < identifier.size(); ++i){
    unsigned int n = identifier[i] - minIndex;
    store[n].push_back(partitioner[i]);
  }

  return store;

}

// get unique values according to an indexer
vector<unsigned int> idSortUnique(vector<unsigned int>& identifier, vector<unsigned int>& values){

  unsigned int minIndex = *min_element(identifier.begin(), identifier.end());
  unsigned int maxIndex = *max_element(identifier.begin(), identifier.end());
  unsigned int nlayers =  maxIndex - minIndex + 1;

  vector<unsigned int> store(nlayers);

  for(unsigned int i = 0; i < identifier.size(); ++i){
    unsigned int n = identifier[i] - minIndex;
    store[n] = values[i];
  }

  return store;

}

vector<double> idSortUnique(vector<unsigned int>& identifier, vector<double>& values){

  unsigned int minIndex = *min_element(identifier.begin(), identifier.end());
  unsigned int maxIndex = *max_element(identifier.begin(), identifier.end());
  unsigned int nlayers =  maxIndex - minIndex + 1;

  vector<double> store(nlayers);

  for(unsigned int i = 0; i < identifier.size(); ++i){
    unsigned int n = identifier[i] - minIndex;
    store[n] = values[i];
  }

  return store;

}

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

// fit ransac circle
vector<double> ransacCircle(vector<vector<double> >& cloud, unsigned int nSamples, double pConfidence, double pInliers){

  unsigned int kIterations = ceil(5 * log(1 - pConfidence) / log(1 - pow( pInliers, nSamples)));
  vector< vector<double> > allCircles( 4, vector<double>(kIterations) );
  unsigned int best = 0;

  Eigen::Matrix<double, Eigen::Dynamic, 3> tempMatrix;
  tempMatrix.resize(nSamples, 3);

  Eigen::Matrix<double, Eigen::Dynamic, 1> rhsVector;
  rhsVector.resize(nSamples, 3);

  for(unsigned int k = 0; k < kIterations; ++k){

    vector<unsigned int> random(nSamples);

    for(unsigned int i = 0; i < random.size(); ++i){
      unsigned int n;
      bool exists;

      do{
        n = floor( R::runif(0, cloud[0].size()) );
        exists = find(begin(random), end(random), n) != end(random);
      }while(exists);

      random[i] = n;

      tempMatrix(i, 0) = cloud[0][n];
      tempMatrix(i, 1) = cloud[1][n];
      tempMatrix(i, 2) = 1;
      rhsVector(i,0) = pow( tempMatrix(i,0), 2) + pow( tempMatrix(i,1), 2);
    }

    Eigen::Matrix<double, 3, 1> qrDecompose = tempMatrix.colPivHouseholderQr().solve(rhsVector);

    Eigen::Matrix<double, 3, 1> xyr;
    xyr(0,0) =  qrDecompose(0,0) / 2;
    xyr(1,0) =  qrDecompose(1,0) / 2;
    xyr(2,0) =  sqrt( ((pow( qrDecompose(0,0) ,2) + pow( qrDecompose(1,0) ,2)) / 4) + qrDecompose(2,0) );

    double sumOfSquares = 0;
    for(unsigned int i = 0; i < cloud[0].size(); ++i){
      double tempX = pow( cloud[0][i] - xyr(0,0) , 2);
      double tempY = pow( cloud[1][i] - xyr(1,0) , 2);
      sumOfSquares += pow( sqrt( tempX + tempY ) - xyr(2,0) , 2);
    }

    double circleError = sqrt(sumOfSquares / cloud[0].size());

    allCircles[0][k] = xyr(0,0);
    allCircles[1][k] = xyr(1,0);
    allCircles[2][k] = xyr(2,0);
    allCircles[3][k] = circleError;

    if(allCircles[3][k] < allCircles[3][best])
      best = k;

  }

  vector<double> bestFit = { allCircles[0][best] , allCircles[1][best] , allCircles[2][best] , allCircles[3][best] };

  return bestFit;

}

// fit ransac circles on stem cloud
vector<vector<double> > ransacStemCircles(vector<vector<double> >& cloud, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

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
    vector<double> temp = ransacCircle(slice, nSamples, pConfidence, pInliers);

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

// fit ransac circles over many stems
vector<vector<vector<double> > > ransacPlotCircles(vector<vector<double> >& cloud, vector<unsigned int>& treeId, vector<unsigned int>& segments, vector<double>& radii, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){

  vector<vector<vector<double> > > trees = getChunks(cloud, treeId);

  cloud.clear();
  cloud.shrink_to_fit();

  vector<unsigned int> uniqId = idSortUnique(treeId, treeId);
  vector<vector<unsigned int> > indices = partitionIndex(treeId, segments);
  vector<vector<double> > treeRadii = partitionIndex(treeId, radii);

  vector< vector< vector<double> > > treeEstimates;

  for(unsigned int i = 0; i < trees.size(); ++i){

    vector<unsigned int>& segs = indices[i];

    if(segs.empty()) continue;

    vector< vector<double> >& tree = trees[i];
    vector<double>& segsRadii = treeRadii[i];

    vector< vector<double> > temp = ransacStemCircles(tree, segs, segsRadii, nSamples, pConfidence, pInliers, tolerance);

    for(vector< vector<double> >::iterator t = temp.begin(); t != temp.end(); t++)
      t->push_back(uniqId[i]);

    treeEstimates.push_back(temp);
  }

  return treeEstimates;

}
