# require(TreeLS)
# require(magrittr)
# require(lidR)
# require(rgl)
# # open artificial sample file
# file = 'test_data/Parcela.las'
# tls = readTLS(file, filter='-thin_with_voxel 0.025')
#
# # normalize the point cloud
# tls = tlsNormalize(tls, keepGround = T)
#
# # extract the tree map from a systematically sampled point cloud
# thin = tlsSample(tls, 'voxel', 0.025)
# map = treeMap(thin)
#
# # visualize the tree map in 2D and 3D
# xymap = treePositions(map, TRUE)
# plot(map, color='Radii')
#
# # classify the stem points
# tls = stemPoints_plot(tls, map)
#
# a = lasfilter(tls, Stem)
#
# plot(tls, color='Stem', size=.5)
# plot(a, color='TreeID', clear_artifacts=F)
