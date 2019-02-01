#include <Rcpp.h>
#include <math.h>
#include <iostream>

#include "classes.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

Raster getCounts(NumericMatrix& slice , float pixel_size, double x_mid=0, double y_mid=0, double d_mid=-1){

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
