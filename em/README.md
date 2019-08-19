# Segmentation Details

Current segmentation is performed using watershed segmentation. Initial segmentation result usually has too many regions, one for each local maxima in density map. To get rid of this, a scale-space grouping method is performed in order to reach the desired result, which is a region for each molecule domain. 

## Parameters

* Smoothing steps: Specifies the number of Gaussian smoothing steps performed to original density map.
* Step size: Standard deviation in voxels to perform smoothing.


# Segger Comparison
## EMD-1010
### Step size = 1, Steps = 3
Our Method             |  Segger
:-------------------------:|:-------------------------:
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1010.png" width="512">  |  <img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1010-chimera.png" width="512">
4 regions | 4 regions


## EMD-1364
### Step size = 1, Steps = 3
Our Method             |  Segger
:-------------------------:|:-------------------------:
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1364.png" width="512">  |  <img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1364-chimera.png" width="512">
4 regions | 4 regions

## EMD-5017
### Step size = 1, Steps = 3
Our Method             |  Segger
:-------------------------:|:-------------------------:
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/5017.png" width="512">  |  <img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/5017-chimera.png" width="512">
3 regions | 2 regions

## EMD-1048
### Step size = 1, Steps = 5
Our Method             |  Segger
:-------------------------:|:-------------------------:
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1048.png" width="512">  |  <img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/1048-chimera.png" width="512">
111 regions | 186 regions

## EMD-2596
### Step size = 1, Steps = 4
Our Method             |  Segger
:-------------------------:|:-------------------------:
<img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/2596.png" width="512">  |  <img src="https://github.com/tecdatalab/biostructure/raw/em-visualizer/em/data/2596-chimera.png" width="512">
55 regions | 94 regions

