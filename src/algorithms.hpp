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

#endif // ALGORITHMS_HPP