#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <unordered_set>

using namespace std;

// [[Rcpp::export]]
void unique_vec(){
  vector<unsigned int> uq = {1,2,3,4,5,1,2,3,4,5,6,6,6,1,1};
  unordered_set<unsigned int> uqs(uq.begin(), uq.end());
  cout << uq.size() << " : " << uqs.size() << endl;
}

