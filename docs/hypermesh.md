# Hypermesh: Meshing The High-Dimensional Spaces for Data Analysis and Visualization

![](./images/regular_simplex_subdivision.svg)

Hypermesh is a header-only C++17 library that meshes high-dimensional spaces for data analysis and visualization, such as feature extraction and tracking, flow map computation, and topology analysis. 

The library is still under development.  Stay tuned. 

## Basic Components

* High-dimensional mesh and mesh element data structures
* High-dimensional array

## Supported Mesh Types

* Regular grid mesh that consists of *n*-cubes in arbitrary dimensions
* Tessellated regular grid mesh that consists of *n*-simplices in arbitrary dimensions
* (TODO) Unstructured *n*-simplex mesh in arbitrary dimensions
* (TODO) *n*-simplex-*m*-prism mesh in arbitrary dimensions, e.g.
  * 2-simplex-1-prism mesh that generalizes given triangular mesh into 3D
  * 3-simplex-1-prism mesh that generalize given tetrahedral mesh into 4D


