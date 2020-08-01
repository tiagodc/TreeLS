#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>

using namespace std;

// [[Rcpp::export]]
vector<double> pointDistances(){
  vector<double> dists;

  vector<vector<double> > cloud = {
      {1,2,3,4,5},
      {8,9,7,6,5},
      {12,13,34,45,5}
  };

  if(cloud[0].size() == 1){
    dists = {0};
    return dists;
  }

  for(unsigned int i = 0; i < cloud[0].size(); ++i){
    for(unsigned int j = i+1; j < cloud[0].size(); ++j){
      double sumsq = 0;
      for(unsigned int k=0; k < cloud.size(); ++k){ sumsq += pow(cloud[k][j] - cloud[k][i],2); }     
      dists.push_back(sqrt(sumsq));      
    }
  }

  return dists;
}
