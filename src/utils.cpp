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

#include "utils.hpp"

// statistics
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

void eigenDecomposition(vector<vector<double> >& cloud, vector<double>* eiVals, vector<vector<double> >* eiVecs) {

  arma::mat armaCloud;
  armaCloud.insert_cols(0, arma::vec(cloud[0]));
  armaCloud.insert_cols(1, arma::vec(cloud[1]));
  armaCloud.insert_cols(2, arma::vec(cloud[2]));

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

// conversions
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
