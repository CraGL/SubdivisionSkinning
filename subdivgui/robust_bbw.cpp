#include "robust_bbw.h"
#include "clean.h"

#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>
#include <igl/tetgen/mesh_with_skeleton.h>
#include <igl/boundary_conditions.h>
#include <igl/bbw/bbw.h>
#include <igl/REDRUM.h>

bool robust_bbw(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  Eigen::MatrixXd & W)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  // clean mesh
  MatrixXd CV;
  MatrixXi CF;
  VectorXi IM;
  if(!clean(V,F,CV,CF,IM))
  {
    return false;
  }
  MatrixXd TV;
  MatrixXi TT;
  // compute tet-mesh
  {
    MatrixXi _1;
#ifdef VERBOSE
    cerr<<"mesh_with_skeleton"<<endl;
    //writeOBJ("mesh_with_skeleton.obj",CV,CF);
#endif
    if(!mesh_with_skeleton(CV,CF,C,{},BE,{},10,"pq1.5Y",TV,TT,_1))
    {
      cout<<REDRUM("tetgen failed.")<<endl;
      return false;
    }
    //writeMESH("mesh_with_skeleton.mesh",TV,TT,MatrixXi());
  }
  // Finally, tetgen may have still included some insanely small tets.
  // Just ignore these during weight computation (and hope they don't isolate
  // any vertices).
  {
    const MatrixXi oldTT = TT;
    VectorXd vol;
    volume(TV,TT,vol);
    const double thresh = 1e-17;
    const int count = (vol.array()>thresh).cast<int>().sum();
    TT.resize(count,oldTT.cols());
    int c = 0;
    for(int t = 0;t<oldTT.rows();t++)
    {
      if(vol(t)>thresh)
      {
        TT.row(c++) = oldTT.row(t);
      }
    }
  }

  // compute weights
  VectorXi b;
  MatrixXd bc;
  if(!boundary_conditions(TV,TT,C,{},BE,{},b,bc))
  {
    cout<<REDRUM("boundary_conditions failed.")<<endl;
    return false;
  }
  // compute BBW
  // Default bbw data and flags
  BBWData bbw_data;
  bbw_data.verbosity = 1;
#ifdef IGL_NO_MOSEK
  bbw_data.qp_solver = QP_SOLVER_IGL_ACTIVE_SET;
  bbw_data.active_set_params.max_iter = 10;
#else
  bbw_data.mosek_data.douparam[MSK_DPAR_INTPNT_TOL_REL_GAP]=1e-10;
  bbw_data.qp_solver = QP_SOLVER_MOSEK;
#endif
#ifdef VERBOSE
  cerr<<"bbw"<<endl;
#endif
  // Weights matrix
  if(!bbw(TV,TT,b,bc,bbw_data,W))
  {
    return false;
  }
  //writeDMAT("W.dmat",W);
  //writeDMAT("IM.dmat",(IM.array()+1).eval());
  // Normalize weights to sum to one
  normalize_row_sums(W,W);
  MatrixXd oldW = W;
  for(int i = 0;i<IM.size();i++)
  {
    if(IM(i)<V.rows())
    {
      W.row(IM(i)) = oldW.row(i);
    }
  }
  // keep first #V weights
  W.conservativeResize(V.rows(),W.cols());
  return true;
}

