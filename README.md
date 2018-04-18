# Scale Space Construction on 3D Surfaces

We diffuse a signal (mean curvature) on explicit and implicit surfaces. We use the [mesh Laplacian](https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP) for explicit surface diffusion and the [Implicit Closest Point Representation](https://people.maths.ox.ac.uk/macdonald/icpm.pdf) for implicit surface diffusion.

We propose a version of Gaussian for explicit surfaces called the mesh Gaussian. We also extend the method of Gaussian diffusion on implicit surfaces first theorized by [Merriman and Ruuth](ftp://ftp.math.ucla.edu/pub/camreport/cam06-32.pdf).


## Getting Started

ext/ directory contains external libraries I've included for the your convenience

src/ directoy contains support functions used by the main programs

main/ directoy contains the main programs

models/ is empty by design, but will populate as you run programs. It is easier to load saved matrices (_e.g._ mesh Laplacian), than always constructing them from scratch. models/ directory will build for specific surfaces
* Plane
* Circle
* Sphere
* General Models


### Prerequisites

You will need a working copy of [MATLAB](https://www.mathworks.com/products/matlab.html). and [gcc 4.7](https://gcc.gnu.org/gcc-4.7/)

### Installing

Only the mesh Laplacian from ext/MeshLP/ needs to be mexed in Matlab. While in Matlab, follow these commands

```
cd ext/MeshLP/

mexcommands.m
```

Both the mesh and the cotangent Laplacian will compile. Now your ready to run some example code.


## Authors
* [Andrew Rhodes](https://github.com/AndrewRhodes) - Initial Implementation - West Virginia University Ph.D. Student

## External Libraries

* [Colin Macdonald](https://github.com/cbm755) - [Implicit Closest Point Representation](https://people.maths.ox.ac.uk/macdonald/icpm.pdf) at [Github](https://github.com/cbm755/cp_matrices)

* [Jian Sun](https://github.com/areslp) - [mesh Laplacian](http://www.cs.jhu.edu/~misha/Fall09/Belkin08.pdf) at [Github](https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP)

* Wil O.C. Ward - [icosphere.m](https://www.mathworks.com/matlabcentral/fileexchange/50105-icosphere)

* [Gabriel Peyr](https://github.com/gpeyre) - [mesh file readers](https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_graph)



## Future Work

Currently all code is written in MATLAB, though I would like to someday convert it other more open source languages.
