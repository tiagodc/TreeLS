#### v2.0.0
- naming convention changed for all methods used inside other functions - all start with a 3 letter prefix (`ptm`, `map`, `ptm`, `trp`...).
- some methods were renamed: `tlsTransform`, `treeMap.positions`.
- new methods added por tree mapping, stem points and stem segmentation, contemplating several new possible use cases and forest scenarios (complex forests, high accuracy LiDAR sensors, irregular stem shapes and so on).
- `treeMap.merge` function introduced and applied internally in `treeMap`, allowing detection and merging of tree duplicates and forked stems into the same tree. 
- extra step added to the inventory workflow and new method added: `treePoints`. This way whole trees are assigned TreeIDs in a plot, not just their stem points.
- new method added for quick forest inventory metrics: `tlsInventory`.
- new functionalities introduced: `fastPointMetrics`, `shapeFit`, `nnFilter`, `writeTLS`, `shapeFit.forks`
- all RANSAC methods have been tuned internally for extra robustness and stability when estimating diameters of noisy point clouds
- `tlsPlot` completely refactored to dynamically interpret any unnamed inputs and plot them accordingly. 
- `add_*` functions added, making plotting of 3D products much more flexible and compatible with `lidR` standards.
- actual 3D cylinder plotting added for all stem segmentation methods
- dependency issues fixed updated
- `lidR` internal calls updated with new naming standards
- bug fixes and performance enhancements

#### v1.0.1
- bug fix #3