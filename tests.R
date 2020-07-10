Rcpp::sourceCpp('src/r_interface.cpp', rebuild = T)
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')
source('R/plot_extras.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('./test_data/Parcela.las', filter='-keep_random_fraction 0.2')

benchmarkme::get_ram()
as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo",intern=TRUE))

sz = ls() %>% sapply(function(x) get(x) %>% object.size) %>% sum
sz / 1000000
a = gc(full = T)
rm(las)
gc()
