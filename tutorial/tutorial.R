require(rgl)
source('https://raw.githubusercontent.com/tiagodc/TreeLS/master/stem_extract/base_functions.R')
source('https://raw.githubusercontent.com/tiagodc/TreeLS/master/stem_extract/pre_sets.R')

path = 'https://raw.githubusercontent.com/tiagodc/TreeLS/master/tutorial/'

trn = 1 # 1 == spruce , 2 == pine

tree = read.table(url(paste(path,'SW_tree_',trn,'.xyz',sep='')), header = T)
tree = Vfilter(tree, l.int = .1, thr = 500)

proof = read.table(paste(path, 'tutorial_field.txt', sep=''), header = T)


#M1

M1.trunk = pref_HT(tree = tree, l.int = .5, cell.size = .025, min.den = .3)

M1.stem = fit_RANSAC_circle(trunk = M1.trunk, l.int = .5, cut.rad = .01)

bg3d('black')
rgl.points(M1.stem, col='blue')
rgl.points(M1.trunk, col='red')
rgl.points(tree, col='grey', size=1)

M1.conic = taper.mod(stem = M1.stem, tree = tree, l.int = .5, method = 'circle', rad.max = .3)
M1.eval = real.vs.est(tree = tree, stem = M1.stem, 
                      d = proof$d[proof$tree == trn], H = proof$h[proof$tree == trn], 
                      model = M1.conic, rad.max = .3, method = 'circle', cyl.len = .5, plot.dd = T)


#M2

M2.trunk = pref_SD(tree = tree, k = 20, flat.min = .9, ang.tol = 10, l.int = .5, freq.ratio = .2)

M2.stem = fit_IRTLS(trunk = M2.trunk, c.len = .5, max.rad = .3, s.height = 1, speed.up = T)

bg3d('black')
rgl.points(M2.stem, col='blue')
rgl.points(M2.trunk, col='red')
rgl.points(tree, col='grey', size=1)

M2.conic = taper.mod(stem = M2.stem, tree = tree, l.int = .5, method = 'cylinder', rad.max = .3)
M2.eval = real.vs.est(tree = tree, stem = M2.stem,
                      d = proof$d[proof$tree == trn], H = proof$h[proof$tree == trn],
                      model = M2.conic, rad.max = .3, method = 'cylinder', cyl.len = .5, plot.dd = T)


#M3

M3.trunk = pref_VN(tree = tree, noise1.rad = .05, noise2.rad = .1, flat.min = .9, ang.tol = 10, neighborhood = 4, largest.cov = NULL, axis.dist = .5)

M3.stem = fit_RANSAC_cylinder(trunk = M3.trunk, c.len = .5, h.init = 1, max.rad = .3, timesN = 2)

bg3d('black')
rgl.points(M3.stem, col='blue')
rgl.points(M3.trunk, col='red')
rgl.points(tree, col='grey', size=1)

M3.conic = taper.mod(stem = M3.stem, tree = tree, l.int = .5, method = 'cylinder', rad.max = .3)
M3.eval = real.vs.est(tree = tree, stem = M3.stem,
                      d = proof$d[proof$tree == trn], H = proof$h[proof$tree == trn],
                      model = M3.conic, rad.max = .3, method = 'cylinder', cyl.len = .5, plot.dd = T)
