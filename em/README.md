# Segmentation Details

Current segmentation is performed using watershed segmentation. Initial segmentation result usually has too many regions, one for each local maxima in density map. To get rid of this, a scale-space grouping method is performed in order to reach the desired result, which is a region for each molecule domain. 

## Parameters

* Smoothing steps: Specifies the number of smoothing steps performed to original density map.
* Step size: Size in voxels for


## EMD-1010
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1010%20.png" width="256">


## EMD-1364
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1364%20.png" width="256">

## EMD-5017
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/5017.png" width="256">

## EMD-2596
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/2596.png" width="256">
