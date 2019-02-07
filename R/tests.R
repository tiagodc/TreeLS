require(magrittr)
require(lidR)
require(rgl)

las = 'd:/Projects/TLStools/test_clouds/gerdau.laz'

las %<>% readLAS(select='XYZ', filter='-thin_with_voxel 0.025')

# plot(las)

m = treeMap(las@data %>% as.matrix, pixel = 0.025, min_den = 0.05, min_votes = 3, hstep = 0.5, hmin = 0, hmax = 3) %>% do.call(what=cbind) %>% as.data.frame
m$Intensity %<>% as.integer
m$Keypoint_flag %<>% as.logical
m$PointSourceID %<>% as.integer
m %<>% LAS

plot(m, color='Radii', clear_artifacts=F, size=3)

clear3d()
rgl.points(las@data, size=.5)
rgl.points(m@data, col='red')

dim(m@data)
m %<>% lasfilterduplicates()
dim(m@data)

require(TreeLS)
a = cyl(100000, 1, 0.4, 0.01)
clear3d()
rgl.points(a)

px = 0.025
getCircle(a, px) %>% as.double
b = a %>% makeRaster(cell.size = px) %>% hough(pixel_size = px)
b$centers[ which.max(b$centers[,4]), ] %>% as.double

