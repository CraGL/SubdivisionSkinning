Posing tool for "Skinning Cubic Bezier Splines and Catmull-Clark Subdivision
Surfaces"

# Compiling
This code depends on:

 - libigl
   - eigen
   - stdlib
   - cgal
   - tetgen
   - mosek (optional)
 - libsubdivision_skinning
  - OpenSubDiv

To compile libigl and its dependencies see the [libigl
instructions](https://github.com/libigl/libigl). Remember that libigl is a
_header only_ library. No compilation is necessary. CGAL and tetgen must be
compiled. Academics can install Mosek for free, but this is optional. Just
defined the flag `-DIGL_NO_MOSEK` to disable Mosek support when building this
project.

# Running
Prepare a coarse subdivision surface cage as a .obj file `cage.obj` and a skeleton using
the [.tgf file
format](http://igl.ethz.ch/projects/libigl/file-formats/tgf.html) of libigl
`skeleton.tgf`.
You could use the `libigl/examples/skeleton-builder/` example to build your
skeleton with a GUI.

Run this program for the first time to compute bounded biharmonic weights:

    ./subdivgui cage.obj skeleton.tgf weights.dmat

This will attempt to _clean_ the model by meshing self-intersections, and
filling holes. Then it will compute a tetrahedral mesh of the surface's solid
volume. Finally it will compute bounded biharmonic weights for each skeleton
bone. When finished the weights will be saved to the provided path
`weights.dmat` in the [.dmat file
format](http://igl.ethz.ch/projects/libigl/file-formats/dmat.html) of libigl.

The GUI will now start and you may follow the standard output instructions to
interact and change visualization settings.

Subsequent runs of 

    ./subdivgui cage.obj skeleton.tgf weights.dmat

will load weights from `weights.dmat` rather than recompute them.
