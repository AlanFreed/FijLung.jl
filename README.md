# FijLung

This module provides the deformation gradient *F<sub>ij</sub>* at three locations in the right lung of an individual who has suffered trauma  produced by waves traversing and reflecting across their thorax cavity caused by a high-energy impact to the rib cage from a blunt object. The 1 direction points from the chest towards the left arm. The 2 direction points from the torso towards the head along the spine. And the 3 direction points from the spine towards the breast bone.

These data are used in a book currently being written by the author of this software, i.e.,

Freed, A.D., and Clayton, J.D., *Application of Laplace Stretch in Alveolar Mechanics: A Case Study of Blast and Blunt Trauma*.

To use this module, you will need to add the following private repositories to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/FijLung.jl")
```
## Locations

These locations lie on the same horizontal plane of a right lung, roughly placed midway between the top and bottom surfaces of a lung, i.e., they lie within a 1-3 plane.

1. Location 1:

Just under the visceral pleura, adjacent to the location of impact.

2. Location 2:

Deep within the right lung.

3. Location 3:

Next to a bronchiole tube.

Location 1 is closest to the point of impact, while Location 3 is the furthest from this point.

## Functions

The following functions supply the raw data that are fit with splines. These data, which reside on [github](https::/github.com/AlanFreed/FijLung.jl), are downloaded at runtime.

The following functions return an array of physical scalars (`PhysicalFields.ArrayOfPhysicalScalars`) and an array of physical tensors (`PhysicalFields.ArrayOfPhysicalTensors`) that hold the time *t* (in seconds) and the nine components of a deformation gradient *F<sub>ij</sub>*, respectively, at the 1<sup>st</sup> location sequenced over an interval of 34 milliseconds.

```
function t_loc1()::ArrayOfPhysicalScalars
function F_loc1()::ArrayOfPhysicalTensors
```

The following functions returnsan array of physical scalars (`PhysicalFields.ArrayOfPhysicalScalars`) and an array of physical tensors (`PhysicalFields.ArrayOfPhysicalTensors`) that hold the time *t* (in seconds) and the nine components of a deformation gradient *F<sub>ij</sub>*, respectively, at the 2<sup>nd</sup> location sequenced over an interval of 34 milliseconds.

```
function t_loc2()::ArrayOfPhysicalScalars
function F_loc2()::ArrayOfPhysicalTensors
```

The following functions return an array of physical scalars (`PhysicalFields.ArrayOfPhysicalScalars`) and an array of physical tensors (`PhysicalFields.ArrayOfPhysicalTensors`) that hold the time *t* (in seconds) and the nine components of a deformation gradient *F<sub>ij</sub>*, respectively, at the 3<sup>rd</sup> location sequenced over an interval of 34 milliseconds.

```
function t_loc3()::ArrayOfPhysicalScalars
function F_loc3()::ArrayOfPhysicalTensors
```

## Type

This data structure holds data that have been fit with splines describing the nine individual components of a deformation gradient *F<sub>ij</sub>*, plus its first d*F<sub>ij</sub>*/d*t* and second d<sup>2</sup>*F<sub>ij</sub>*/d*t*<sup>2</sup> derivatives in time.

```
struct SplineF
    loc::Integer                  # loc       Lung location, i.e., 1, 2 or 3.
    N::  Integer                  # N         Number of nodes spanning time.
    t::  ArrayOfPhysicalScalars   # t         The time at these nodes.
    F::  ArrayOfPhysicalTensors   # F         Their deformation gradients.
    F′:: ArrayOfPhysicalTensors   # dF/dt     Their first derivatives in time.
    F′′::ArrayOfPhysicalTensors   # d²F/dt²   Their second derivatives in time.
end
```

### Constructor

The constructor used to create an instance of `SplineF` takes on the form
```
splineF = SplineF(location::Integer, nodes::Integer)
```
which creates a new data structure that holds data for time (`splineF.t`), the deformation gradient (`splineF.F`), its first derivative in time (`splineF.F′`) and its second derivative in time (`splineF.F′′`) at a specified `location`, viz., `splineF.loc` = 1, 2 or 3. Thses splines interpolate the raw data with a density of `splineF.N = nodes,` uniformly spaced in time. There are 340 knots, or raw data points, that the B-spline fits. Typically, the number of `nodes` is a value much greater than this, e.g., 1000. The deformation gradient `SplineF.F` comes from a cubic B-spline of these raw data. Its first derivative `SplineF.F′` is therefore described by a quadratic B-spline, while its second derivative `SplineF.F′′` is described by a linear B-spline.

### Test

An example of how these data can be used, e.g., plotted, can be found in the *test* directory in file *testFijLung.jl*.