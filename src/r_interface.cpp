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

// export tree positions point stack
List exportTreeMap(vector<HoughCenters>& coordinates){

  vector<double> xout;
  vector<double> yout;
  vector<double> zout;
  vector<double> radii;
  vector<bool> keyFlag;
  vector<bool> treeFlag;
  vector<unsigned short int> votes;
  vector<unsigned int> treeId;
  vector<unsigned int> discId;
  // vector<unsigned int> nPoints;

  unsigned int diskCounter = 1;
  unsigned int maxId = 0;

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    if(point.tree_id > maxId)
      maxId = point.tree_id;

    double z = (point.low_z + point.up_z)/2;

    // main point
    xout.push_back(point.main_circle.x_center);
    yout.push_back(point.main_circle.y_center);
    zout.push_back(z);

    radii.push_back(point.main_circle.radius);
    votes.push_back(point.main_circle.n_votes);
    treeId.push_back(point.tree_id);
    // nPoints.push_back(point.circles.size());
    discId.push_back(diskCounter);
    keyFlag.push_back(true);
    treeFlag.push_back(false);

    // other candidates
    for(auto& c_point : point.circles){

      xout.push_back(c_point.x_center);
      yout.push_back(c_point.y_center);
      zout.push_back(z);

      radii.push_back(c_point.radius);
      votes.push_back(c_point.n_votes);
      treeId.push_back(point.tree_id);
      // nPoints.push_back(point.circles.size());
      discId.push_back(diskCounter);
      keyFlag.push_back(false);
      treeFlag.push_back(false);
    }
    diskCounter++;
  }

  vector<double> xSums(maxId, 0);
  vector<double> ySums(maxId, 0);
  vector<unsigned int> counters(maxId, 0);

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    xSums[point.tree_id-1] += point.main_circle.x_center;
    ySums[point.tree_id-1] += point.main_circle.y_center;
    counters[point.tree_id-1]++;
  }

  for(unsigned int i = 0; i < maxId; ++i){

    if(counters[i] == 0)
      continue;

    double mainX = xSums[i] / counters[i];
    double mainY = ySums[i] / counters[i];

    xout.push_back(mainX);
    yout.push_back(mainY);
    zout.push_back(0);

    radii.push_back(0);
    votes.push_back(0);
    treeId.push_back(i+1);
    // nPoints.push_back(0);
    discId.push_back(0);
    keyFlag.push_back(true);
    treeFlag.push_back(true);
  }

  List out;
  out["X"] = xout;
  xout.clear();
  xout.shrink_to_fit();

  out["Y"] = yout;
  yout.clear();
  yout.shrink_to_fit();

  out["Z"] = zout;
  zout.clear();
  zout.shrink_to_fit();

  out["Intensity"] = votes;
  votes.clear();
  votes.shrink_to_fit();

  out["PointSourceID"] = discId;
  discId.clear();
  discId.shrink_to_fit();

  out["Keypoint_flag"] = keyFlag;
  keyFlag.clear();
  keyFlag.shrink_to_fit();

  out["Radii"] = radii;
  radii.clear();
  radii.shrink_to_fit();

  out["TreeID"] = treeId;
  treeId.clear();
  treeId.shrink_to_fit();

  out["TreePosition"] = treeFlag;
  treeFlag.clear();
  treeFlag.shrink_to_fit();

  // out["n"] = nPoints;

  return out;
}

// [[Rcpp::export]]
LogicalVector thinCloud(NumericMatrix& las, double voxel = 0.025){
  vector<vector<double> > xyz = rmatrix2cpp(las);
  return wrap(voxelFilter(xyz, voxel));
}

// [[Rcpp::export]]
LogicalVector RCropCloud(NumericMatrix& las, double xCenter, double yCenter, double len, bool circle, bool negative){
  vector<vector<double> > xyz = rmatrix2cpp(las);
  return wrap( cropCloudFilter(xyz, xCenter, yCenter, len, circle, negative) );
}

// [[Rcpp::export]]
SEXP getHoughCircle(NumericMatrix& las, double pixel=0.05, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

  vector<vector<double> > cloud = rmatrix2cpp(las);
  Raster ras = getCounts(cloud, pixel);
  HoughCenters circle = getSingleCenter(&ras, rad_max, min_den, min_votes);

  vector<vector<double> > ledger;
  for(auto& i : circle.circles){
    vector<double> temp = {i.x_center, i.y_center, i.radius, double(i.n_votes)};
    ledger.push_back(temp);
  }

  // List out;
  // out["x"] = circle.main_circle.x_center;
  // out["y"] = circle.main_circle.y_center;
  // out["rad"] = circle.main_circle.radius;
  // out["votes"] = circle.main_circle.n_votes;

  return wrap( ledger );
}

// [[Rcpp::export]]
List singleStack(NumericMatrix& las, double pixel=0.05, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

  vector<HoughCenters> treeMap;
  Raster ras;

  vector<vector<double> > stack = rmatrix2cpp(las);
  ras = getCounts(stack, pixel);
  treeMap = getCenters(&ras, rad_max, min_den, min_votes);
  assignTreeId(treeMap, rad_max, min_den);

  return exportTreeMap(treeMap);

}

// [[Rcpp::export]]
List stackMap(NumericMatrix& las, double hmin=1, double hmax=3, double hstep=0.5, double pixel=0.025, double rad_max=0.25, double min_den=0.1, unsigned int min_votes = 3){

  vector<HoughCenters> treeMap;

  vector<vector<vector<double> > > fullStack = getSlices(las, hmin, hmax, hstep);

  for(auto& stack : fullStack){

    if(stack[0].empty())
      continue;

    Raster ras = getCounts(stack, pixel);
    vector<HoughCenters> tempMap = getCenters(&ras, rad_max, min_den, min_votes);
    treeMap.insert(treeMap.end(), tempMap.begin(), tempMap.end());
  }

  unsigned int nlayers = 0.75 * (hmax - hmin) / hstep;
  assignTreeId(treeMap, rad_max, min_den, nlayers);

  fullStack.clear();
  fullStack.shrink_to_fit();

  return exportTreeMap(treeMap);

}

// [[Rcpp::export]]
List houghStemPoints(NumericMatrix& las, double h1 = 1, double h2 = 3, double hstep=0.5, double radius=0.25, double pixel=0.025, double density=0.1, unsigned int votes=3){

  vector<vector<double> > cppCloud = rmatrix2cpp(las);
  vector<HoughCenters> treeEstimates = treeHough(cppCloud, h1, h2, hstep, radius, pixel, density, votes);

  if(treeEstimates.empty()){
    List noTree;
    return noTree;
  }

  tempContainer isStem(cppCloud[0].size());
  for(unsigned int i = 0; i < cppCloud[0].size(); ++i){

    double& x = cppCloud[0][i];
    double& y = cppCloud[1][i];
    double& z = cppCloud[2][i];

    if(z < 0) continue;

    unsigned int ptLayer = floor(z / hstep);

    HoughCircle* alias = &treeEstimates[ptLayer].main_circle;

    if(alias->n_votes < votes) continue;

    double dist = sqrt( pow(x - alias->x_center, 2) + pow(y - alias->y_center, 2) );

    if(dist < alias->radius + pixel*2 /* && dist > alias->radius - pixel*2 */){
      isStem.filter[i] = true;
      isStem.values[i] = alias->radius;
      isStem.counts[i] = alias->n_votes;
      isStem.sections[i] = ptLayer + 1;
    }
  }

  cppCloud.clear();
  cppCloud.shrink_to_fit();

  List output;
  output["Stem"] = isStem.filter;
  output["Segment"] = isStem.sections;
  output["Radius"] = isStem.values;
  output["Votes"] = isStem.counts;
  isStem.clear();

  return output;
}

// [[Rcpp::export]]
List houghStemPlot(NumericMatrix& las, NumericVector& ptIds, double h1 = 1, double h2 = 3, double hstep=0.5, double radius=0.25, double pixel=0.025, double density=0.1, unsigned int votes=3){

  // unordered_set<unsigned int> treeIds(pointIds.begin(), pointIds.end());
  // vector<unsigned int> treeIdsVec(treeIds.begin(), treeIds.end());

  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> pointIds = Rcpp::as< vector<unsigned int> >( ptIds );
  vector<vector<vector<double> > > treeList = getChunks(cloud, pointIds);
  // unordered_map<unsigned int, vector<HoughCenters> > denoisedTrees;
  vector<vector<HoughCenters>> denoisedTrees;

  for(auto& tree : treeList){

    // vector<vector<double> >& tree = treeList[i];
    if(tree[0].empty()) continue;

    vector<HoughCenters> denoised = treeHough(tree, h1, h2, hstep, radius, pixel, density, votes);

    if(denoised.empty()) continue;

    denoisedTrees.push_back(denoised);
  }

  tempContainer plotInfo( cloud[0].size() );
  for(unsigned int i = 0; i < cloud[0].size(); ++i){

    double& x = cloud[0][i];
    double& y = cloud[1][i];
    double& z = cloud[2][i];

    if(z < 0) continue;

    unsigned int ptLayer = floor(z / hstep);
    for(auto& tree : denoisedTrees){

      if(tree.size() <= ptLayer) continue;

      HoughCircle* tempCircle = &tree[ptLayer].main_circle;

      if(tempCircle->n_votes < votes) continue;

      double dist = sqrt( pow(x - tempCircle->x_center, 2) + pow(y - tempCircle->y_center, 2) );

      if(dist < tempCircle->radius + pixel*2 /* && dist > tempCircle->radius - pixel*2 */){
        plotInfo.filter[i] = true;
        plotInfo.values[i] = tempCircle->radius;
        plotInfo.counts[i] = tempCircle->n_votes;
        // plotInfo.ids[i] = tree.first;
        plotInfo.sections[i] = ptLayer + 1;
        break;
      }
    }
  }

  cloud.clear();
  cloud.shrink_to_fit();

  List output;
  output["Stem"]   = plotInfo.filter;
  // output["TreeID"] = plotInfo.ids;
  output["Segment"] = plotInfo.sections;
  output["Radius"] = plotInfo.values;
  output["Votes"]  = plotInfo.counts;
  plotInfo.clear();

  return output;

}

// [[Rcpp::export]]
NumericVector getCircleRansac(NumericMatrix& las, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  return wrap(ransacCircle(cloud, nSamples, pConfidence, pInliers));
}

// [[Rcpp::export]]
List ransacStemCircle(NumericMatrix& las, NumericVector& segs, NumericVector& rads, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<double> radii = Rcpp::as< vector<double> >( rads );
  return wrap(ransacStemCircle(cloud, segments, radii, nSamples, pConfidence, pInliers, tolerance));
}

// [[Rcpp::export]]
List irlsStemCylinder(NumericMatrix& las, NumericVector& segs, NumericVector& rads, unsigned int nPoints=500,  double tolerance=0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<double> radii = Rcpp::as< vector<double> >( rads );
  return wrap(irlsStemCylinder(cloud, segments, radii, nPoints, tolerance));
}

// [[Rcpp::export]]
List irlsStemCircle(NumericMatrix& las, NumericVector& segs, NumericVector& rads, unsigned int nSamples=500, double tolerance=0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(irlsStemCircle(cloud, segments, radii, nSamples, tolerance));
}

// [[Rcpp::export]]
List ransacStemCylinder(NumericMatrix& las, NumericVector& segs, NumericVector& rads, unsigned int nSamples=10, double pConfidence=0.95, double pInliers=0.8, double tolerance=0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(ransacStemCylinder(cloud, segments, radii, nSamples, pConfidence, pInliers, tolerance));
}

// [[Rcpp::export]]
List ransacPlotCircles(NumericMatrix& las, NumericVector& tId, NumericVector& segs, NumericVector& rads, unsigned int nSamples = 5, double pConfidence = 0.99, double pInliers = 0.8, double tolerance = 0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<unsigned int> treeId = Rcpp::as< vector<unsigned int> >( tId );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(ransacPlotCircles(cloud, treeId, segments, radii, nSamples, pConfidence, pInliers, tolerance));
}

// [[Rcpp::export]]
List ransacPlotCylinders(NumericMatrix& las, NumericVector& tId, NumericVector& segs, NumericVector& rads, unsigned int nSamples, double pConfidence, double pInliers, double tolerance){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<unsigned int> treeId = Rcpp::as< vector<unsigned int> >( tId );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(ransacPlotCylinders(cloud, treeId, segments, radii, nSamples, pConfidence, pInliers, tolerance));
}

// [[Rcpp::export]]
List irlsPlotCylinders(NumericMatrix& las, NumericVector& tId, NumericVector& segs, NumericVector& rads, unsigned int nPoints, double tolerance){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<unsigned int> treeId = Rcpp::as< vector<unsigned int> >( tId );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(irlsPlotCylinders(cloud, treeId, segments, radii, nPoints, tolerance));
}

// [[Rcpp::export]]
List irlsPlotCircles(NumericMatrix& las, NumericVector& tId, NumericVector& segs, NumericVector& rads, unsigned int nPoints, double tolerance){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<unsigned int> treeId = Rcpp::as< vector<unsigned int> >( tId );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(irlsPlotCircles(cloud, treeId, segments, radii, nPoints, tolerance));
}

// [[Rcpp::export]]
SEXP pointMetricsCpp(NumericMatrix& las, NumericMatrix& kIds, LogicalVector& whichMetrics){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<vector<unsigned int> > idx = intmatrix2cpp(kIds);
  vector<bool> wmt = Rcpp::as< vector<bool> >( whichMetrics );
  return wrap( pointMetrics(cloud, idx, wmt) );
}

// [[Rcpp::export]]
SEXP voxelIndex(NumericMatrix& las, double d){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  return wrap( voxelIndex(cloud, d) );
}

// [[Rcpp::export]]
List voxelMetrics(NumericMatrix& las, List& voxelIds, LogicalVector& whichMetrics){
  typedef vector< vector<unsigned int> > vvint;
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vvint idx = Rcpp::as< vvint >( voxelIds );
  vector<bool> wmt = Rcpp::as< vector<bool> >( whichMetrics );
  return wrap( voxelMetrics(cloud, idx, wmt) );
}

// [[Rcpp::export]]
SEXP treeEigenHough(NumericMatrix& las, NumericVector& ids, NumericVector& split_by, double voxel=0.05, double rad = 0.25, bool is2d = false, bool getSpace = false){
  vector<vector<double> > xyz = rmatrix2cpp(las);
  vector<unsigned int> stdIds(ids.begin(), ids.end());
  vector<unsigned int> splitter(split_by.begin(), split_by.end());

  return wrap( treeEigenHough(xyz, stdIds, splitter, voxel, rad, is2d, getSpace) );
}

// [[Rcpp::export]]
SEXP plotEigenHough(NumericMatrix& las, NumericVector& ids, NumericVector& split_by, NumericVector& resplit_by, double voxel=0.05, double rad = 0.25, bool is2d = false, bool getSpace = false){

  vector<vector<double> > xyz = rmatrix2cpp(las);
  vector<unsigned int> stdIds(ids.begin(), ids.end());
  vector<unsigned int> splitter1(split_by.begin(), split_by.end());
  vector<unsigned int> splitter2(resplit_by.begin(), resplit_by.end());

  return wrap( plotEigenHough(xyz, stdIds, splitter1, splitter2, voxel, rad, is2d, getSpace) );
}

// [[Rcpp::export]]
SEXP cppFastApply(NumericMatrix& matrix, StringVector& funcList){
  vector<vector<double> > cppmat = rmatrix2cpp(matrix);
  vector<string> funcs(funcList.begin(), funcList.end());
  return wrap( fastApply(cppmat, funcs) );
}

// [[Rcpp::export]]
SEXP cppCircleFit(NumericMatrix& las, std::string method = "qr", unsigned int n = 5, double p = 0.99, double inliers = 0.8, unsigned int nbest = 0){

  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<double> pars = {0};

  if(method == "irls"){
    pars = irlsCircleFit(las);
  }else if(method == "qr"){
    pars = eigenCircle(cloud);
  }else if(method == "nm"){
    pars = nmCircleFit(cloud);
  }else if(method == "ransac"){
    pars = ransacCircle(cloud, n, p, inliers, nbest);
  }

  return wrap( pars );
}

// [[Rcpp::export]]
SEXP cppCylinderFit(NumericMatrix& las, std::string method = "nm", unsigned int n = 10, double p = 0.95, double inliers = 0.9, double max_angle = 30, unsigned int n_best = 20){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<double> pars;

  double nmax = 100;
  if(method != "ransac" && method != "bf" && cloud[0].size() > nmax){
    double prop = nmax / (double)cloud[0].size();
    cloud = randomPoints(cloud, prop);
  }

  if(method == "irls"){
    vector<double> initPars = {0, PI/2, 0, 0, 0};
    pars = irlsCylinder(cloud, initPars);
  }else if(method == "nm"){
    pars = nmCylinderFit(cloud);
  }else if(method == "ransac"){
    pars = ransacCylinder(cloud, n, p, inliers);
  }else if(method == "bf"){
    pars = bruteForceRansacCylinder(cloud, n, p, inliers, n_best, max_angle, true)[0];
  }

  return wrap( pars );
}

// [[Rcpp::export]]
SEXP treeIdsFromMap(NumericMatrix& las, NumericMatrix& xycenters, NumericVector& uniqueIds, double length, bool circle){
  vector<vector<double> > xy = rmatrix2cpp(las);
  vector<vector<double> > xymap = rmatrix2cpp(xycenters);
  vector<unsigned int> ids = Rcpp::as< vector<unsigned int> >( uniqueIds );
  return wrap( treeIdsFromMap(xy, xymap, ids, length, circle) );
}

// [[Rcpp::export]]
SEXP bruteForceRansacCylinder(NumericMatrix& las, unsigned int nSamples, double pConfidence, double pInliers, unsigned int nBest, double maxAngle){
  vector<vector<double> > xyz = rmatrix2cpp(las);
  return wrap( bruteForceRansacCylinder(xyz, nSamples, pConfidence, pInliers, nBest, maxAngle) );
}

// [[Rcpp::export]]
List bfStemCylinder(NumericMatrix& las, NumericVector& segs, NumericVector& rads, unsigned int nSamples=10, double pConfidence=0.95, double pInliers=0.8, double max_angle = 30, double tolerance=0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(bfStemCylinder(cloud, segments, radii, nSamples, pConfidence, pInliers, max_angle, tolerance));
}

// [[Rcpp::export]]
List bfPlotCylinders(NumericMatrix& las, NumericVector& tId, NumericVector& segs, NumericVector& rads, unsigned int nSamples = 10, double pConfidence = 0.95, double pInliers = 0.8, double max_angle = 30, double tolerance = 0.05){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<unsigned int> segments = Rcpp::as< vector<unsigned int> >( segs );
  vector<unsigned int> treeId = Rcpp::as< vector<unsigned int> >( tId );
  vector<double> radii = Rcpp::as< vector<double> >( rads );  
  return wrap(bfPlotCylinders(cloud, treeId, segments, radii, nSamples, pConfidence, pInliers, max_angle,tolerance));
}