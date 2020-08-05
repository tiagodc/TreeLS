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

#include <Eigen/Dense>

#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include "utils.hpp"

vector<double> circleDists(vector<vector<double> >& xyz, arma::vec& pars);

vector<double> cylDists(vector<vector<double> >& xyz, arma::vec& pars);

vector<double> eigenCircle(vector<vector<double> >& cloud);

vector<double> irlsCircle(vector<vector<double> >& las, vector<double> initPars, double err_tol = 1E-06, unsigned int max_iter = 100);

vector<double> irlsCircleFit(NumericMatrix& cloud);

vector<double> irlsCylinder(vector<vector<double> >& las, vector<double> initPars, double err_tol = 1E-06, unsigned int max_iter = 100);

double nmCircleDist(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);

vector<double> nmCircleFit(vector<vector<double> >& las);

double nmCylinderDist(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);

vector<double> nmCylinderFit(vector<vector<double> >& las);

vector<double> nmCylinderInit(vector<vector<double> >& las);

vector<double> ransacCircle(vector<vector<double> >& cloud, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, unsigned int nBest = 0);

vector<double> ransacCylinder(vector<vector<double> >& las, unsigned int nSamples=10, double pConfidence=0.99, double pInliers=0.8);

Eigen::Matrix<double, Eigen::Dynamic, 3> stl2eigenmat(vector<vector<double> >& xyz);

Eigen::Matrix<double, 3, 3> rotationMatrix(double ax, double ay, double az);

vector<vector<double> > eigenmat2stl(Eigen::Matrix<double, Eigen::Dynamic, 3>& mat);

vector<vector<double> > rotateCloud(vector<vector<double> >& xyz, double ax, double ay, double az);

vector<vector<double> > bruteForceRansacCylinder(vector<vector<double> >& cloud, unsigned int nSamples, double pConfidence, double pInliers, unsigned int nBest, double maxAngle = 45.0, bool bestOnly = false);

#endif // ALGORITHMS_HPP