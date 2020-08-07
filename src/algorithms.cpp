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

// [[Rcpp::depends(RcppEigen)]]
#include <Eigen/Dense>

#ifndef USE_RCPP_ARMADILLO
  #define USE_RCPP_ARMADILLO
#endif

#include "optim.hpp"

#include "algorithms.hpp"

// eigen conversions
Eigen::Matrix<double, Eigen::Dynamic, 3> stl2eigenmat(vector<vector<double> >& xyz){
  Eigen::Matrix<double, Eigen::Dynamic, 3> eigenMat;
  eigenMat.resize(xyz[0].size(), 3);

  for(unsigned int i = 0; i < 3; ++i){
    eigenMat.col(i) = Eigen::VectorXd::Map(&xyz[i][0],xyz[i].size());
  }

  return eigenMat;
}

vector<vector<double> > eigenmat2stl(Eigen::Matrix<double, Eigen::Dynamic, 3>& mat){

  vector<vector<double> > xyz( 3, vector<double>(mat.rows()) );

  for(unsigned int i = 0; i < mat.cols(); ++i){
    Eigen::VectorXd::Map(&xyz[i][0], mat.rows()) = mat.col(i);
  }

  return xyz;
}

Eigen::Matrix<double, 3, 3> rotationMatrix(double ax, double ay, double az){

  Eigen::Matrix<double, 3, 3> Rx;
  Rx << 1, 0, 0,
        0, cos(ax), sin(ax),
        0, -sin(ax), cos(ax);

  Eigen::Matrix<double, 3, 3> Ry;
  Ry << cos(ay), 0, -sin(ay),
        0, 1, 0,
        sin(ay), 0, cos(ay);

  Eigen::Matrix<double, 3, 3> Rz;
  Rz << cos(az), sin(az), 0,
        -sin(az), cos(az), 0,
        0, 0, 1;

  Eigen::Matrix<double, 3, 3> rotMat = Rz * Ry * Rx;
  return rotMat;
}

vector<vector<double> > rotateCloud(vector<vector<double> >& xyz, double ax, double ay, double az){

  Eigen::Matrix<double, Eigen::Dynamic, 3> xyzMat = stl2eigenmat(xyz) * rotationMatrix(ax, ay, az);

  vector<vector<double> > rot = eigenmat2stl(xyzMat);

  return rot;
}

vector<vector<double> > bruteForceRansacCylinder(vector<vector<double> >& cloud, unsigned int nSamples, double pConfidence, double pInliers, unsigned int nBest, double maxAngle, bool bestOnly){

  vector<vector<double> > pars;
  if(bestOnly) pars.resize(1);
  double err = -1;

  for(double i = -maxAngle; i <= maxAngle; i += 1){
    for(double j = -maxAngle; j <= maxAngle; j += 1){
      vector<vector<double> > tempCloud = rotateCloud(cloud, i * PI/180, j * PI/180, 0.0);
      vector<double> ijPars = ransacCircle(tempCloud, nSamples, pConfidence, pInliers, nBest);

      ijPars.push_back(i);
      ijPars.push_back(j);

      if(bestOnly){
        if(err < 0 || ijPars[3] < err){
          err = ijPars[3];
          pars[0] = ijPars;
        }
      }else{
        pars.push_back(ijPars);
      }
    }
  }

  return pars;
}

// circle methods
vector<double> eigenCircle(vector<vector<double> >& cloud){
  unsigned int n = cloud[0].size();

  Eigen::Matrix<double, Eigen::Dynamic, 3> tempMatrix;
  tempMatrix.resize(n, 3);

  Eigen::Matrix<double, Eigen::Dynamic, 1> rhsVector;
  rhsVector.resize(n, 3);

  for(unsigned int i = 0; i < n; ++i){
    tempMatrix(i, 0) = cloud[0][i];
    tempMatrix(i, 1) = cloud[1][i];
    tempMatrix(i, 2) = 1;
    rhsVector(i,0) = pow( tempMatrix(i,0), 2) + pow( tempMatrix(i,1), 2);
  }

  Eigen::Matrix<double, 3, 1> qrDecompose = tempMatrix.colPivHouseholderQr().solve(rhsVector);

  Eigen::Matrix<double, 3, 1> xyr;
  xyr(0,0) =  qrDecompose(0,0) / 2;
  xyr(1,0) =  qrDecompose(1,0) / 2;
  xyr(2,0) =  sqrt( ((pow( qrDecompose(0,0) ,2) + pow( qrDecompose(1,0) ,2)) / 4) + qrDecompose(2,0) );

  double sumOfSquares = 0;
  for(unsigned int i = 0; i < n; ++i){
    double tempX = pow( cloud[0][i] - xyr(0,0) , 2);
    double tempY = pow( cloud[1][i] - xyr(1,0) , 2);
    sumOfSquares += pow( sqrt( tempX + tempY ) - xyr(2,0) , 2);
  }

  double circleError = sqrt(sumOfSquares / n);

  vector<double> bestFit = { xyr(0,1), xyr(1,0), xyr(2,0), circleError };
  return bestFit;
}

vector<double> ransacCircle(vector<vector<double> >& cloud, unsigned int nSamples, double pConfidence, double pInliers, unsigned int nBest){

  unsigned int kTimes = nBest + 5;
  unsigned int kIterations = kTimes * ceil(log(1 - pConfidence) / log(1 - pow( pInliers, nSamples)));
  vector< vector<double> > allCircles( 4, vector<double>(kIterations) );
  unsigned int best = 0;
  if(kIterations < nBest) nBest = 0;

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
    // for(unsigned int i = 0; i < cloud[0].size(); ++i){
    for(auto& i : random){
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

  vector<double> bestFit;

  if(nBest > 0){
    // double x = 0.0;
    // double y = 0.0;
    // double rad = 0.0;
    // double err = 0.0;
    unsigned int counter = 0;
    vector<double> x, y, rad, err;

    for(auto& k : sortIndexes(allCircles[3])){
      // x += allCircles[0][k];
      // y += allCircles[1][k];
      // rad += allCircles[2][k];
      // err += allCircles[3][k];

      x.push_back( allCircles[0][k] );
      y.push_back( allCircles[1][k] );
      rad.push_back( allCircles[2][k] );
      err.push_back( allCircles[3][k] );

      if(counter++ > nBest){
        // bestFit = { x / counter, y / counter, rad / counter, err / counter };
        bestFit = { median(x), median(y), median(rad), median(err) };
        break;
      }
    }

  }else{
    bestFit = { allCircles[0][best] , allCircles[1][best] , allCircles[2][best] , allCircles[3][best] };
  }

  return bestFit;

}

vector<double> irlsCircle(vector<vector<double> >& las, vector<double> initPars, double err_tol, unsigned int max_iter){

  arma::vec init(initPars);

  vector<double> weights(las[0].size(), 1);
  las.push_back(weights);

  double ssq = 0;
  unsigned int count = 0;
  bool converge = false;

  while(!converge){
    bool success = optim::nm(init,nmCircleDist,&las);
    if(!success) break;
    double err = nmCircleDist(init, nullptr, &las);
    vector<double> werr = circleDists(las, init);
    tukeyBiSq(werr);
    las[3] = werr;
    converge = abs(err - ssq) < err_tol || ++count == max_iter;
    ssq = err;
  }

  vector<double> pars = arma::conv_to<std::vector<double> >::from(init);
  pars.push_back(ssq);

  return pars;

}

vector<double> irlsCircleFit(NumericMatrix& cloud){
  vector<vector<double> > las = rmatrix2cpp(cloud);
  vector<double> init = eigenCircle(las);
  vector<double> res = irlsCircle(las, init);
  return res;
}

vector<double> nmCircleFit(vector<vector<double> >& las){

  arma::vec init(eigenCircle(las));

  bool success = optim::nm(init,nmCircleDist,&las);

  vector<double> pars = arma::conv_to<std::vector<double> >::from(init);
  double err =  success ? nmCircleDist(init, nullptr, &las) : 0;
  pars.push_back(err);

  return pars;

}

vector<double> circleDists(vector<vector<double> >& xyz, arma::vec& pars){

  double xInit = pars(0);
  double yInit = pars(1);
  double rad = pars(2);

  vector<double> sqDists(xyz[0].size());
  for(unsigned int i = 0; i < xyz[0].size(); ++i){
    double x = xyz[0][i];
    double y = xyz[1][i];

    double dst = sqrt( (x-xInit)*(x-xInit) + (y-yInit)*(y-yInit) ) - rad;
    sqDists[i] = (dst*dst);
  }

  return sqDists;
}

double nmCircleDist(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data){

  vector<vector<double> >* xyz = reinterpret_cast<vector<vector<double> >* >(opt_data);

  double xInit = vals_inp(0);
  double yInit = vals_inp(1);
  double rad = vals_inp(2);

  double distSum = 0;
  for(unsigned int i = 0; i < (*xyz)[0].size(); ++i){
    double x = (*xyz)[0][i];
    double y = (*xyz)[1][i];

    double dst = sqrt( (x-xInit)*(x-xInit) + (y-yInit)*(y-yInit) ) - rad;
    distSum += ((*xyz).size() == 4) ? (dst*dst)*(*xyz)[3][i] : (dst*dst);
  }

  return distSum;
}

// cylinder methods
vector<double> nmCylinderInit(vector<vector<double> >& las){

  // double x0 = median(las[0]);
  // double y0 = median(las[1]);
  // double z0 = median(las[2]);
  // double rho = sqrt(x0*x0 + y0*y0 + z0*z0);

  double rho=0;
  double theta = PI/2;
  double phi = 0;
  double alpha = 0;
  double r = 0;

  vector<double> pars = {rho, theta, phi, alpha, r};

  return pars;
}

vector<double> ransacCylinder(vector<vector<double> >& las, unsigned int nSamples, double pConfidence, double pInliers){

  bringOrigin(las);
  vector<vector<double> > tempCloud(3, vector<double>(nSamples));

  unsigned int kIterations = ceil(5 * log(1 - pConfidence) / log(1 - pow( pInliers, nSamples)));
  vector<double> bestfit;

  for(unsigned int k = 0; k < kIterations; ++k){

    vector<unsigned int> random(nSamples);

    for(unsigned int i = 0; i < random.size(); ++i){
      unsigned int n;
      bool exists;

      do{
        n = floor( R::runif(0, las[0].size()) );
        exists = find(begin(random), end(random), n) != end(random);
      }while(exists);

      random[i] = n;

      tempCloud[0][i] = las[0][n];
      tempCloud[1][i] = las[1][n];
      tempCloud[2][i] = las[2][n];
    }

    vector<double> fit = nmCylinderFit(tempCloud);

    if(k == 0) bestfit = fit;

    if(fit[5] < bestfit[5]){
      bestfit = fit;
    }

  }

  return bestfit;

}

vector<double> irlsCylinder(vector<vector<double> >& las, vector<double> initPars, double err_tol, unsigned int max_iter){

  // vector<vector<double> > las = rmatrix2cpp(cloud);
  bringOrigin(las);
  arma::vec init(initPars);

  vector<double> weights(las[0].size(), 1);
  las.push_back(weights);

  double ssq = 0;
  unsigned int count = 0;
  bool converge = false;

  while(!converge){
    bool success = optim::nm(init,nmCylinderDist,&las);
    double err = nmCylinderDist(init, nullptr, &las);

    vector<double> werr = cylDists(las, init);
    tukeyBiSq(werr);
    las[3] = werr;
    converge = abs(err - ssq) < err_tol || ++count == max_iter;
    ssq = err;
  }

  vector<double> pars = arma::conv_to<std::vector<double> >::from(init);
  pars.push_back(ssq);

  return pars;

}

vector<double> nmCylinderFit(vector<vector<double> >& las){

  bringOrigin(las);
  arma::vec init(nmCylinderInit(las));

  bool success = optim::nm(init,nmCylinderDist,&las);

  vector<double> pars = arma::conv_to<std::vector<double> >::from(init);
  double err =  success ? nmCylinderDist(init, nullptr, &las) : 0;
  pars.push_back(err);

  return pars;

}

vector<double> cylDists(vector<vector<double> >& xyz, arma::vec& pars){

  double rho = pars(0);
  double theta = pars(1);
  double phi = pars(2);
  double alpha = pars(3);
  double r = pars(4);

  vector<double> n = {cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta)};
  vector<double> ntheta = {cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta)};
  vector<double> nphi = {-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0};

  vector<double> nphibar = nphi;
  nphibar[0] /= sin(theta);
  nphibar[1] /= sin(theta);
  nphibar[2] /= sin(theta);

  vector<double> a = {
    ntheta[0] * cos(alpha) + nphibar[0] * sin(alpha),
    ntheta[1] * cos(alpha) + nphibar[1] * sin(alpha),
    ntheta[2] * cos(alpha) + nphibar[2] * sin(alpha)
  };

  vector<double> q = n;
  q[0] *= (r + rho);
  q[1] *= (r + rho);
  q[2] *= (r + rho);

  vector<double> sqDists(xyz[0].size());
  for(unsigned int i = 0; i < xyz[0].size(); ++i){
    double x = xyz[0][i];
    double y = xyz[1][i];
    double z = xyz[2][i];

    vector<double> iq = {x - q[0], y - q[1], z - q[2]};
    vector<double> xp = xprod(iq,a);
    double dst = sqrt( xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2] ) - r ;

    sqDists[i] = (dst*dst);

  }

  return sqDists;

}

double nmCylinderDist(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data){

  vector<vector<double> >* xyz = reinterpret_cast<vector<vector<double> >* >(opt_data);

  double rho = vals_inp(0);
  double theta = vals_inp(1);
  double phi = vals_inp(2);
  double alpha = vals_inp(3);
  double r = vals_inp(4);

  vector<double> n = {cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta)};
  vector<double> ntheta = {cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta)};
  vector<double> nphi = {-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0};

  vector<double> nphibar = nphi;
  nphibar[0] /= sin(theta);
  nphibar[1] /= sin(theta);
  nphibar[2] /= sin(theta);

  vector<double> a = {
    ntheta[0] * cos(alpha) + nphibar[0] * sin(alpha),
    ntheta[1] * cos(alpha) + nphibar[1] * sin(alpha),
    ntheta[2] * cos(alpha) + nphibar[2] * sin(alpha)
  };

  vector<double> q = n;
  q[0] *= (r + rho);
  q[1] *= (r + rho);
  q[2] *= (r + rho);

  double distSum = 0;
  for(unsigned int i = 0; i < (*xyz)[0].size(); ++i){
    double x = (*xyz)[0][i];
    double y = (*xyz)[1][i];
    double z = (*xyz)[2][i];

    vector<double> iq = {x - q[0], y - q[1], z - q[2]};
    vector<double> xp = xprod(iq,a);
    double dst = sqrt( xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2] ) - r ;

    distSum += ((*xyz).size() == 4) ? (dst*dst)*(*xyz)[3][i] : (dst*dst);
  }

  return distSum;
}
