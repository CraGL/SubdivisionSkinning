.PHONY: all

GG=clang++
#CFLAGS += -DIGL_STATIC_LIBRARY
CFLAGS += -Wall -Wno-deprecated-declarations -Wno-deprecated-register
#CFLAGS += -O0
CFLAGS += -O3 -DNDEBUG
CFLAGS+=-std=c++11 -g -DVERBOSE
MOSEKPLATFORM=osx64x86
MOSEKVERSION=7
# msse4.2 is necessary for me to get embree to compile correctly
AFLAGS=-m64 -msse4.2


all: obj subdivgui

.PHONY: subdivgui

$(info Hello, $(LOGNAME))
ifeq ($(LOGNAME),ajx)
	LIB+=/opt/local/lib/gcc47/libstdc++.dylib
	LIBIGL=/usr/local/igl/libigl/
	DEFAULT_PREFIX=/opt/local
	ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar_clang -framework AppKit
else
	IGL_NO_MOSEK=1
	LIBIGL=/Users/yotam/Work/ext/libigl
	DEFAULT_PREFIX=/usr/local
## I don't have Alec's patched GLUT with GLUT_ACTIVE_COMMAND
	CFLAGS += -DGLUT_ACTIVE_COMMAND=8
	ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
endif
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
LIBIGL_INC=-I$(LIBIGL)/include
# LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglembree -liglcgal -ligltetgen -liglbbw -liglmosek
LIBIGL_LIB=-L$(LIBIGL)/lib
EMBREE_INC=-I$(LIBIGL)/external/embree/embree/ -I$(LIBIGL)/external/embree/
EMBREE_LIB=-L$(LIBIGL)/external/embree/build -lembree -lsys

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include

TETGEN=$(LIBIGL)/external/tetgen
TETGEN_LIB=-L$(TETGEN) -ltet 
TETGEN_INC=-I$(TETGEN)

# YIMAGE Library
YIMG=$(LIBIGL)/external/yimg/
YIMG_LIB=-L$(YIMG) -lyimg -lz -L/usr/X11/lib -lpng -bind_at_load
YIMG_INC=-I/usr/X11/include -I$(YIMG) 

CARBON_LIB=-framework Carbon

# Use free glut for mouse scrolling
#FREE_GLUT=/opt/local/
#FREE_GLUT_INC=-I$(FREE_GLUT)/include
#FREE_GLUT_LIB=-L$(FREE_GLUT)/lib -lglut
GLUT_LIB=-framework GLUT -framework OpenGL
GLUT_INC=

ifdef IGL_NO_MOSEK
CFLAGS+=-DIGL_NO_MOSEK
else
# Adjust your mosek paths etc. accordingly
ifndef MOSEKPLATFORM
  MOSEKPLATFORM=osx64x86
endif
ifndef MOSEKVERSION
  MOSEKVERSION=7
endif
IGLMOSEK=../mosek/
IGLMOSEK_INC=-I$(IGLMOSEK)/
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64
endif

CGAL=$(DEFAULT_PREFIX)
CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
CGAL_FLAGS=-frounding-math -fsignaling-nans 
CFLAGS+=$(CGAL_FLAGS)

SUBDIV=.
SUBDIV_INC=-I$(SUBDIV)
SUBDIV_LIB=-L$(SUBDIV) -lsubdivision_skinning -losdCPU

INC+=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC) $(GLUT_INC) ${CGAL_INC} ${TETGEN_INC} $(MOSEK_INC) $(EMBREE_INC) $(SUBDIV_IN) $(YIMG_INC)
LIB+=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(MATLAB_LIB) ${CGAL_LIB} $(CARBON_LIB) $(TETGEN_LIB) $(MOSEK_LIB) $(EMBREE_LIB) $(SUBDIV_LIB) $(YIMG_LIB)

CPP_FILES=$(wildcard ./*.cpp)
#OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 
OBJ_FILES=obj/subdivgui.o obj/clean.o obj/robust_bbw.o obj/subdiv_weights.o obj/Mesh.o


subdivgui: obj $(OBJ_FILES)
	$(GG) $(OPENMP) $(AFLAGS) $(CFLAGS) -o subdivgui $(LIB) $(OBJ_FILES) 

obj:
	mkdir -p obj

obj/%.o: %.cpp
	$(GG) $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	$(GG) $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f subdivgui
