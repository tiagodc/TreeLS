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

#include "utils.hpp"

// statistics
void eigenDecomposition(vector<vector<double> >& cloud, vector<double>* eiVals, vector<vector<double> >* eiVecs){

  arma::mat armaCloud;

  for(unsigned int i = 0; i < cloud.size(); ++i){
    armaCloud.insert_cols(i, arma::vec(cloud[i]));
  }

  arma::mat coeff;
  arma::mat score;
  arma::vec latent;

  arma::princomp(coeff, score, latent, armaCloud);

  *eiVals = arma::conv_to<std::vector<double> >::from(latent);

  *eiVecs = {};
  for(unsigned int i = 0; i < cloud.size(); ++i){
    eiVecs->push_back( arma::conv_to<std::vector<double> >::from(coeff.col(i)) );
  }
}

vector<double> pointDistances(vector<vector<double> >& cloud){
  vector<double> dists;

  if(cloud[0].size() == 1){
    dists = {0};
    return dists;
  }

  for(unsigned int i = 0; i < cloud[0].size(); ++i){
    for(unsigned int j = i+1; j < cloud[0].size(); ++j){
      double sumsq = 0;
      for(unsigned int k = 0; k < cloud.size(); ++k) sumsq += pow(cloud[k][j] - cloud[k][i], 2);
      dists.push_back(sqrt(sumsq));
    }
  }

  return dists;
}

vector<double> nnMetrics(vector<vector<double> >& xyz, vector<bool> which){

      vector<double> eVal3d;
      vector<vector<double> > eVec3d;
      eigenDecomposition(xyz, &eVal3d, &eVec3d);
      vector<double> dists3d = pointDistances(xyz);

      vector<vector<double> > xy = {xyz[0], xyz[1]};
      vector<double> eVal2d;
      vector<vector<double> > eVec2d;
      eigenDecomposition(xy, &eVal2d, &eVec2d);
      vector<double> dists2d = pointDistances(xy);

      vector<double> z = {0,0,1};
      double zmean = accumulate(xyz[2].begin(), xyz[2].end(), 0.0) / xyz[2].size();
      double zmax = *max_element(xyz[2].begin(), xyz[2].end());
      double zmin = *min_element(xyz[2].begin(), xyz[2].end());
      double zsumsq = 0;
      for(auto& pt : xyz[2]) zsumsq += pow( pt - zmean ,2);
      double sumprod = 0;
      for(auto& ev : eVal3d) sumprod += ev * log(ev);
      
      double n_pts = xyz[0].size();
      double min_dist  = *min_element(dists3d.begin(), dists3d.end());
      double max_dist  = *max_element(dists3d.begin(), dists3d.end());
      double mean_dist = accumulate(dists3d.begin(), dists3d.end(), 0.0) / dists3d.size();
      double sd_dist = sqrt(variance(dists3d));
      double linearity = (eVal3d[0] - eVal3d[1]) / eVal3d[0];
      double planarity = (eVal3d[1] - eVal3d[2]) / eVal3d[0];;
      double scattering = eVal3d[2] / eVal3d[0];
      double omnivariance = pow((eVal3d[0] * eVal3d[1] * eVal3d[2]), (1.0/3.0));
      double anisotropy = (eVal3d[0] - eVal3d[2]) / eVal3d[0];
      double eigentropy = -sumprod;
      double eigen_sum_3d = eVal3d[0] + eVal3d[1] + eVal3d[2];
      double surface_variation = eVal3d[2] / eigen_sum_3d;
      double radius_knn_3d = *max_element(dists3d.begin(), dists3d.end()) / 2;
      double density_3d = n_pts / ((4.0/3.0) * M_PI * pow(radius_knn_3d, 3.0));
      double verticality = vecAngle(z, eVec3d[2]);
      double z_range = zmax - zmin;
      double z_sd = sqrt( zsumsq / xyz[2].size() );
      double radius_knn_2d = *max_element(dists2d.begin(), dists2d.end()) / 2;
      double density_2d = n_pts / (M_PI * pow(radius_knn_2d, 2.0));
      double eigen_sum_2d = eVal2d[0] + eVal2d[1];
      double eigen_ratio_2d = eVal2d[1] / eVal2d[0];

      vector<double> metrics = {};
      unsigned int metric_index = 0;

      if(which[metric_index++]) metrics.push_back(n_pts);      
      if(which[metric_index++]) metrics.push_back(min_dist);
      if(which[metric_index++]) metrics.push_back(max_dist);
      if(which[metric_index++]) metrics.push_back(mean_dist);
      if(which[metric_index++]) metrics.push_back(sd_dist);
      if(which[metric_index++]) metrics.push_back(linearity);
      if(which[metric_index++]) metrics.push_back(planarity);
      if(which[metric_index++]) metrics.push_back(scattering);
      if(which[metric_index++]) metrics.push_back(omnivariance);
      if(which[metric_index++]) metrics.push_back(anisotropy);
      if(which[metric_index++]) metrics.push_back(eigentropy);
      if(which[metric_index++]) metrics.push_back(eigen_sum_3d);
      if(which[metric_index++]) metrics.push_back(surface_variation);
      if(which[metric_index++]) metrics.push_back(radius_knn_3d);
      if(which[metric_index++]) metrics.push_back(density_3d);
      if(which[metric_index++]) metrics.push_back(verticality);
      if(which[metric_index++]) metrics.push_back(z_range);
      if(which[metric_index++]) metrics.push_back(z_sd);
      if(which[metric_index++]) metrics.push_back(radius_knn_2d);
      if(which[metric_index++]) metrics.push_back(density_2d);
      if(which[metric_index++]) metrics.push_back(eigen_sum_2d);
      if(which[metric_index++]) metrics.push_back(eigen_ratio_2d);

      for(auto& v : eVal3d){
        if(which[metric_index++]) metrics.push_back(v);
      }

      for(auto& a : eVec3d){
        for(auto& b : a){
          if(which[metric_index++]) metrics.push_back(b);
        }
      }

      return metrics;
}

double mad(vector<double> x, double c){

  double md = median(x);

  for(auto& i : x){
    i = abs(i - md);
  }

  return c * median(x);
}

double median(vector<double> x){
  sort(x.begin(), x.end());
  unsigned int i = round(x.size()/2);
  return x[i];
}

void tukeyBiSq(vector<double>& werrors, double b){  
  double s = mad(werrors);

  for(auto& i : werrors){    
    i /= s;
    i = abs(i) > b ? 0 : pow(1-(i/b)*(i/b),2);
  }
}

vector<double> xprod(vector<double>& a, vector<double>& b){

  vector<double> x = {
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  };

  return x;
}

double vecAngle(vector<double>& a, vector<double>& b){

  double dotprod = 0;
  double asqsum  = 0;
  double bsqsum  = 0;

  for(unsigned int i = 0; i < a.size(); ++i){
    dotprod += a[i]*b[i];
    asqsum  += a[i]*a[i];
    bsqsum  += b[i]*b[i];
  }

  double ang  = dotprod / ( sqrt(asqsum) * sqrt(bsqsum) );
  double cang = acos(ang) * 180/M_PI;
  return cang;
}

double variance(vector<double>& x){
  double mean = accumulate(x.begin(), x.end(), 0.0) / x.size();
  double var = 0.0;
  for(auto& i : x){
    var += pow(i - mean, 2);
  }
  var /= x.size();
  return var;
}

// conversions
vector<vector<double> > rmatrix2cpp(NumericMatrix& cloud){

  unsigned int ncol = cloud.ncol();
  vector<vector<double> > xyz(ncol);

  for(unsigned int i = 0; i < ncol; ++i){
    NumericMatrix::Column tempcol = cloud( _, i);
    xyz[i].insert(xyz[i].begin(), tempcol.begin(), tempcol.end());
  }

  return xyz;
}

vector<vector<unsigned int> > intmatrix2cpp(NumericMatrix& idxMatrix){

  unsigned int ncol = idxMatrix.ncol();
  vector<vector<unsigned int> > imat(ncol);

  for(unsigned int i = 0; i < ncol; ++i){
    NumericMatrix::Column tempcol = idxMatrix( _, i);
    imat[i].insert(imat[i].begin(), tempcol.begin(), tempcol.end());
  }

  return imat;
}

// point cloud manipulation
void bringOrigin(vector<vector<double> >& las){

  double x0 = median(las[0]);
  double y0 = median(las[1]);
  double z0 = median(las[2]);

  for(unsigned int i = 0; i < las[0].size(); ++i){
    las[0][i] -= x0;
    las[1][i] -= y0;
    las[2][i] -= z0;
  }

}

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

vector<vector<double> > randomPoints(vector<vector<double> >& cloud, double p){

  vector<vector<double> > reCloud(3);

  for(unsigned int i = 0; i < cloud[0].size(); ++i){

    double val = R::runif(0,1);
    if(val > p) continue;

    reCloud[0].push_back( cloud[0][i] );
    reCloud[1].push_back( cloud[1][i] );
    reCloud[2].push_back( cloud[2][i] );
  }

  return reCloud;

}

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

vector<vector<double> > voxelCounter(vector<vector<double> >& xyzNormals, double voxel, double max_rad, bool is2d, bool sendSpace){

  typedef unsigned long long int llint;

  double xmin = *min_element(xyzNormals[0].begin(), xyzNormals[0].end()) - max_rad;
  double ymin = *min_element(xyzNormals[1].begin(), xyzNormals[1].end()) - max_rad;
  double zmin = *min_element(xyzNormals[2].begin(), xyzNormals[2].end()) - max_rad;

  VoxelGrid voxelRegistry(xmin, ymin, zmin, voxel);

  vector<double> origCounts;
  vector<double> origDists;
  unsigned int nk = sendSpace ? 1 : 2;

  for(unsigned int k = 0; k < nk; ++k){
    for(unsigned int i = 0; i < xyzNormals[0].size(); ++i){

      double x  = xyzNormals[0][i];
      double y  = xyzNormals[1][i];
      double z  = xyzNormals[2][i];
      double e1 = xyzNormals[3][i];
      double e2 = xyzNormals[4][i];
      double e3 = xyzNormals[5][i];

      llint centerHash = is2d ? voxelRegistry.getPixelHash(x, y, z) : voxelRegistry.getVoxelHash(x, y, z);
      unordered_set<llint> uniqueIds;
      unsigned int maxCount = 0;
      double mainDist;

      for(double d = -max_rad; d < max_rad + voxel; d += voxel){
        double xtemp = x + d*e1;
        double ytemp = y + d*e2;
        double ztemp = z + d*e3;

        llint hash = is2d ? voxelRegistry.getPixelHash(xtemp, ytemp, ztemp) : voxelRegistry.getVoxelHash(xtemp, ytemp, ztemp);
        if(hash == centerHash) continue;

        if(k == 0){
          bool isFirst = uniqueIds.insert(hash).second;

          if(isFirst){
            if(is2d){
              voxelRegistry.updatePixelRegistry(xtemp, ytemp, ztemp);
            }else{
              voxelRegistry.updateVoxelRegistry(xtemp, ytemp, ztemp);
            }
          }
        }else{
          unsigned int tempCount = voxelRegistry.counter[hash];
          if(tempCount > maxCount){
            maxCount = tempCount;
            mainDist = abs(d);
          }
        }
      }

      if(k == 1){
        origCounts.push_back(maxCount);
        origDists.push_back(mainDist);
      }
    }
  }

  vector<vector<double> > nVoxel;

  if(!sendSpace){

    nVoxel = { origCounts, origDists };

  }else{
    for(auto& vox : voxelRegistry.counter){
      unsigned long long int hash = vox.first;
      array<unsigned int, 3> nxyz;
      array<unsigned int, 2> nxy;

      vector<double> row;
      row.push_back(vox.second);

      if(is2d){
        nxy = voxelRegistry.pixels[hash];
        row.push_back(nxy[0] * voxel + xmin);
        row.push_back(nxy[1] * voxel + ymin);
      }else{
        nxyz = voxelRegistry.voxels[hash];
        row.push_back(nxyz[0] * voxel + xmin);
        row.push_back(nxyz[1] * voxel + ymin);
        row.push_back(nxyz[2] * voxel + zmin);
      }

      nVoxel.push_back(row);
    }
  }

  return nVoxel;
}

//// split point cloud according to some criterion
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

vector<vector<vector<double> > > getFullChunks(vector<vector<double> >& cloud, vector<unsigned int>& identifier){

  set<unsigned int> uniqueIds;
  uniqueIds.insert(identifier.begin(), identifier.end());
  unsigned int nlayers = uniqueIds.size();
  unsigned int ncols = cloud.size();

  vector<vector<vector<double> > > store(nlayers, vector<vector<double> >(ncols));

  for(unsigned int i = 0; i < cloud[0].size(); ++i){
    unsigned int id = identifier[i];
    set<unsigned int>::iterator it = find(uniqueIds.begin(), uniqueIds.end(), id);
    unsigned int pos = distance(uniqueIds.begin(), it);

    for(unsigned int j = 0; j < ncols; ++j){
      store[pos][j].push_back( cloud[j][i] );
    }
  }

  return store;

}

vector<vector<unsigned int> > splitVector(vector<unsigned int>& to_split, vector<unsigned int>& split_by){

  set<unsigned int> usp(split_by.begin(), split_by.end());
  vector<vector<unsigned int> > parts(usp.size());

  for(unsigned int i = 0; i < to_split.size(); ++i){
    unsigned int val = to_split[i];
    unsigned int idx = split_by[i];

    set<unsigned int>::iterator it = find(usp.begin(), usp.end(), idx);
    unsigned int dst = distance(usp.begin(), it);

    parts[dst].push_back(val);
  }

  return parts;

}

//// split point cloud into horizontal slices
vector<vector<vector<double> > > getSlices(NumericMatrix& cloud, double zmin, double zmax, double zstep){
  vector<vector<double> > las = rmatrix2cpp(cloud);
  return getSlices(las, zmin, zmax, zstep);
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

// utilities
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

vector<vector<double> > fastApply(vector<vector<double> >& matrix, vector<string>& funcList){

  vector<vector<double> > calc( matrix[0].size(), vector<double>(funcList.size(), 0) );

  for(unsigned int i = 0; i < matrix[0].size(); ++i){

    vector<double> row;
    for(unsigned int k = 0; k < matrix.size(); ++k){
      double temp = matrix[k][i];
      if(temp < 0.00000001) break;
      row.push_back(temp);
    }

    if(row.empty()) continue;

    unsigned int j = 0;
    for(auto& f : funcList){

      if(f == "MedianDistance"){
        calc[i][j++] = median(row);
      }else if(f == "MinDistance"){
        calc[i][j++] = *min_element(row.begin(), row.end());
      }else if(f == "MaxDistance"){
        calc[i][j++] = *max_element(row.begin(), row.end());
      }else if(f == "MeanDistance"){
        calc[i][j++] = accumulate(row.begin(), row.end(), 0.0) / row.size();
      }else if(f == "VarDistance"){
        calc[i][j++] = variance(row);
      }else if(f == "SdDistance"){
        calc[i][j++] = sqrt(variance(row));
      }
    }
  }

  return calc;

}

vector<unsigned int> sortIndexes(vector<double>& values) {

  vector<unsigned int> idx(values.size());
  iota(idx.begin(), idx.end(), 0);

  sort(idx.begin(), idx.end(), [&values](unsigned int i1, unsigned int i2) {return values[i1] < values[i2];});

  return idx;
}

unsigned int uniqueTotalCounter(vector<unsigned int> input){
  unordered_set<unsigned int> uq(input.begin(), input.end());
  return uq.size();
}

void progressPrinter(string units, unsigned int n, unsigned int total){
  double percentage = 100.0 * (double)(++n) / (double)total;
  Rcpp::Rcout << fixed << setprecision(2) << "Progress (" << percentage << " %): " << n << " of " << total << " " << units << "\r" << flush;
}