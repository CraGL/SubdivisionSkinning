# Skinning Catmull-Clark Subdivision Surfaces

This code implements the subdivision surface portion of our SIGGRAPH Asia 2014 paper
[Skinning Cubic Bézier Splines and Catmull-Clark Subdivision Surfaces](http://cs.gmu.edu/~ygingold/splineskin/)
by [Songrun Liu](http://cs.gmu.edu/~sliu11/), [Alec Jacobson](http://www.cs.columbia.edu/~jacobson/), and [Yotam Gingold](http://cs.gmu.edu/~ygingold/).
The implementation consists of a library and a GUI for posing subdivision surfaces with skeletons

## Compiling

This code depends on:

- [OpenSubDiv](http://graphics.pixar.com/opensubdiv) 2.x
- [libigl](https://github.com/libigl/libigl)
    - included sub-dependencies:
        - AntTweakBar
        - embree
        - tetgen
        - yimg
        - tetgen
    - not included sub-dependencies:
        - CGAL (e.g. `brew install cgal`)
    - mosek (optional)
- [eigen](http://eigen.tuxfamily.org/) (e.g. `brew install eigen`)

Academics can install Mosek for free, but this is optional.
The CMakeFiles define the flag `-DIGL_NO_MOSEK` to disable Mosek support when building this
project.

## Download and compile the OpenSubdiv 2.x dependency
    git clone https://github.com/PixarAnimationStudios/OpenSubdiv.git
    (
        cd OpenSubdiv
        &&
        git reset --hard 25cee425f32758a7e0e8812da628007a8eeecce6
        &&
        mkdir build
        &&
        cd build
        &&
        cmake -DCMAKE_BUILD_TYPE=Release ..
        &&
        make
        )


## Download libigl and compile the third-party dependencies
    git clone https://github.com/libigl/libigl.git
    ( cd libigl/external/AntTweakBar/src && make )
    ( cd libigl/external/embree && make build && cd build && cmake .. && make )
    ( cd libigl/external/embree && make build && cd build && cmake .. && make )
    ( cd libigl/external/tetgen && make )
    ( cd libigl/external/yimg && make )


## Compile this project (libsubdivision_skinning and subdivgui)
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
    cd ..

# Using the library

The library is compiled into
    build/lib/libsubdivision_skinning.{a,dylib}

The library interface is exposed through a pure-C interface. The header file describing the functions is in
    lib/subdivision_skinning_wrapper.h

The included GUI is a usage example: `subdivgui/subdivgui.cpp`

# Running the included GUI
Prepare a coarse subdivision surface cage as a .obj file `cage.obj` and a skeleton using
the [.tgf file
format](http://igl.ethz.ch/projects/libigl/file-formats/tgf.html) of libigl
`skeleton.tgf`.
You could use the `libigl/examples/skeleton-builder/` example to build your
skeleton with a GUI.

Run this program for the first time to compute bounded biharmonic weights:

    ./build/subdivgui/subdivgui cage.obj skeleton.tgf weights.dmat

You can run the included torus example:
    
    ./build/subdivgui/subdivgui subdivgui/torus{.obj,.tgf,-weights.dmat}

This will attempt to _clean_ the model by meshing self-intersections, and
filling holes. Then it will compute a tetrahedral mesh of the surface's solid
volume. Finally it will compute bounded biharmonic weights for each skeleton
bone. When finished the weights will be saved to the provided path
`weights.dmat` in the [.dmat file
format](http://igl.ethz.ch/projects/libigl/file-formats/dmat.html) of libigl.

The GUI will now start and you may follow the standard output instructions to
interact and change visualization settings.

Subsequent runs of 

    ./build/subdivgui/subdivgui cage.obj skeleton.tgf weights.dmat

will load weights from `weights.dmat` rather than recompute them.
