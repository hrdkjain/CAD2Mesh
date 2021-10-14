# CAD2Mesh
This project contains matlab script for conversion of CAD model to manifold surface mesh. Parts of the code have been taken from the following repositories:
- [Learning Geometry Images](https://github.com/sinhayan/learning_geometry_images)
- [Skeleton3D](https://de.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d)

## Usage
#### Download
```
git clone https://github.com/hrdkjain/Matlab-Functions.git
git clone https://github.com/hrdkjain/CAD2Mesh.git
```
#### Environment Requirement 
The project has been tested on Ubuntu 20.04 with following requirements:
- Matlab (tested on version R2021a)
- meshlabserver (tested on version [2020.09](https://github.com/cnr-isti-vclab/meshlab/releases/download/Meshlab-2020.09/MeshLabServer2020.09-linux.AppImage))
#### Project Dependencies
This project depends on [Matlab-Functions](https://github.com/hrdkjain/m_libaries).
#### Execution
- Copy `include.m.example` to `include.m` and modify the paths. 
- The `main.m` script reads CAD models stored in `srcDir`, converts them to CAD to mesh and saves them to `dstDir`. 
- After performing the processing, there might be few non-manifold edges left in the surface mesh, which are removed by `meshlabserver`. 
- Often to run meshlabserver from Matlab, `LD_LIBRARY_PATH` environment variable needs to be set to qt library path. Use this command to find the qt tool path `qtchooser -print-env | grep QTTOOLDIR`. 

##  Projects
If you like this project you also might be interested in other projects which use CAD2Mesh as preprocessing stage for [3D surface mesh reconstruction](https://github.com/hrdkjain/LearningSymmetricShapes) or [mesh generation](https://github.com/hrdkjain/GenIcoNet).