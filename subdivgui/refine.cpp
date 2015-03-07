#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
export PS4=""
set -o xtrace
TEMP="$0.cpp"
#clang++ -O3 -std=c++11 -o .main $TEMP \
#g++ -O3 -std=c++11 -o .main $TEMP \
#  -fopenmp \
printf "//" | cat - $0 >$TEMP
clang++ -g -O3 -DNDEBUG -std=c++11 -o .main $TEMP \
  -I. libosdCPU.a libsubdivision_skinning.a \
  -I/opt/local/include/eigen3/ -I/usr/local/igl/libigl/include \
  -I/usr/local/igl/libigl/external/AntTweakBar/include \
  -I/Applications/MATLAB_R2014a.app/extern/include/ \
  -I/usr/local/igl/libigl/external/AntTweakBar/src \
  -I/usr/local/igl/libigl/external/glfw/include \
  -L/usr/local/igl/libigl/external/glfw/lib -lglfw3 -framework QuartzCore -framework IOKit \
  -I/usr/local/igl/libigl/external/Singular_Value_Decomposition/ \
  -I/opt/local/include/ -I/usr/include/ \
  -L/opt/local/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt \
  -L/Applications/MATLAB_R2014a.app/extern/lib/ \
  -framework OpenGL \
  -framework GLUT \
  -framework AppKit \
  -L/opt/local/lib -lboost_thread-mt -lboost_system-mt \
  -L/Applications/MATLAB_R2014a.app/bin/maci64/ -lmx -lmat -lmex -leng \
  -L/usr/local/igl/libigl/external/AntTweakBar/lib -lAntTweakBar \
  -L/opt/local/lib -lboost_program_options-mt && ./.main "$@"
#-DIGL_STATIC_LIBRARY -L/usr/local/igl/libigl/lib -liglviewer -ligl -liglmatlab -liglsvd3x3 -msse4.2 -fopenmp \
#rm -f .main
rm -f $TEMP
# END BASH SCRIPT
exit
*/

// This is a tiny program to read in a coarse mesh and write out a refined
// mesh.

#include <igl/readOBJ.h>
#include <igl/REDRUM.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/pathinfo.h>
#include <igl/file_exists.h>
#include <Eigen/Core>
#include <iostream>


extern "C"
{
#include "subdivision_skinning_wrapper.h"
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  MatrixXd CV,V;
  MatrixXi CQ,Q;

  int level = 3;

  string filename = "torus.obj", output_filename = "";
  switch(argc)
  {
    case 4:
      output_filename = argv[3];
    case 3:
      level = atoi(argv[2]);
    case 2:
      // Read and prepare mesh
      filename = argv[1];
      break;
    default:
      cerr<<"Usage:"<<endl<<"    ./example model.obj [level] [ouput.obj]"<<endl;
      cout<<endl<<"Opening default mesh..."<<endl;
  }

  // Load subdivision cage
  string dir,_1,_2,name;
  pathinfo(filename,dir,_1,_2,name);
  {
    vector<vector<int > > vCQ;
    vector<vector<double> > vCV;
    if(!readOBJ(filename,vCV,vCQ))
    {
      return EXIT_FAILURE;
    }
    list_to_matrix(vCV,CV);
    if(!list_to_matrix(vCQ,4,-1,CQ))
    {
      cerr<<"Error: "<<filename<<" contains high valence facet."<<endl;
      return EXIT_FAILURE;
    }
  }

  if(output_filename.size() == 0)
  {
    output_filename = dir+"/"+name+"-output.obj";
    if(file_exists(output_filename.c_str()))
    {
      cout<<YELLOWGIN("Output set to overwrite "<<output_filename)<<endl;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // Tell subdiv lib about cage
  /////////////////////////////////////////////////////////////////////////
  void* eval;
  const Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> CVT = 
    CV.cast<subdivision_evaluator_real_t>().transpose();
  const subdivision_evaluator_real_t * vs = CVT.data();
  const int num_vs = CV.rows();
  const MatrixXi CQT = CQ.cast<int>().transpose();
  const int * faces = CQT.data();
  const int num_faces = CQ.rows();
  eval = new_subdivision_evaluator( 
    num_vs, (subdivision_evaluator_real_t*)vs, num_faces, (int*)faces,level);
    
  /////////////////////////////////////////////////////////////////////////
  // Ask subdiv lib for refined mesh
  /////////////////////////////////////////////////////////////////////////
  const int num_refined_faces = num_refined_quad_faces_of_subdivision_evaluator( eval );
  int* refined_faces = new int[ num_refined_faces*4 ];
  get_refined_quad_faces_of_subdivision_evaluator( eval, refined_faces );
  const auto & QT = Map<MatrixXi>(refined_faces,4,num_refined_faces);
  Q = QT.transpose();
  delete[] refined_faces;
    
  const int num_refined_vs = num_refined_vertices_of_subdivision_evaluator( eval );
  subdivision_evaluator_real_t* refined_vs = new subdivision_evaluator_real_t[ num_refined_vs*3 ];
  get_refined_vertices_of_subdivision_evaluator( eval, refined_vs );
  const auto & VT = Map<Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> >
    (refined_vs,3,num_refined_vs);
  V = VT.cast<double>().transpose();
  delete[] refined_vs;
    
  if(!writeOBJ(output_filename,V,Q))
  {
    cerr<<"Failed to write to "<<output_filename<<endl;
    return EXIT_FAILURE;
  }

  // Clean up
  delete_subdivision_evaluator( eval );
  return EXIT_SUCCESS;
}
