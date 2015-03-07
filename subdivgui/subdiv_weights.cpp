#include "subdiv_weights.h"
#include "robust_bbw.h"
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/embree/bone_heat.h>

extern "C"
{
#include "subdivision_skinning_wrapper.h"
}

bool subdiv_weights(
  const Eigen::MatrixXd & CV,
  const Eigen::MatrixXi & CQ,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  const int computation_level,
  const int evaluation_level,
  const WeightsType type,
  Eigen::MatrixXd & W)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;

  typedef subdivision_evaluator_real_t Scalar;
  typedef Eigen::Matrix<Scalar,Dynamic,Dynamic> MatrixXS;

  void* computation_eval;
  void* evaluation_eval;

  // 1. Create evaluator from control mesh with level L1 appropriate for
  // computing weights.

  {
    const MatrixXS CVT = CV.cast<Scalar>().transpose();
    const Scalar * vs = CVT.data();
    const int num_vs = CV.rows();
    const MatrixXi CQT = CQ.cast<int>().transpose();
    const int * faces = CQT.data();
    const int num_faces = CQ.rows();
    computation_eval = new_subdivision_evaluator( 
      num_vs, (Scalar*)vs, num_faces, (int*)faces, computation_level);
  }

  // 2. Compute refined mesh R1 (V,Q)
  MatrixXd V;
  MatrixXi Q,F;
  int * refined_faces;
  int num_refined_faces;
  Scalar * refined_vs;
  int num_refined_vs;
  {
    num_refined_faces = 
      num_refined_quad_faces_of_subdivision_evaluator( computation_eval );
    refined_faces = new int[ num_refined_faces*4 ];
    get_refined_quad_faces_of_subdivision_evaluator( 
      computation_eval, refined_faces );
    const auto & QT = Map<MatrixXi>(refined_faces,4,num_refined_faces);
    Q = QT.transpose();
    polygon_mesh_to_triangle_mesh(Q,F);
    num_refined_vs = 
      num_refined_vertices_of_subdivision_evaluator( computation_eval );
    refined_vs = new Scalar[ num_refined_vs*3 ];
    get_refined_vertices_of_subdivision_evaluator(computation_eval,refined_vs);
    const auto & VT = Map<MatrixXS >
      (refined_vs,3,num_refined_vs);
    V = VT.cast<double>().transpose();
  }

  // 3. Compute weights W1 on R1 (V,F)
  MatrixXd W1;
  switch(type)
  {
    default:
      assert(false && "Unknown weighting type.");
    case WEIGHTS_TYPE_BBW:
      if(!robust_bbw(V,F,C,BE,W1))
      {
        return false;
      }
      break;
    case WEIGHTS_TYPE_BONE_HEAT:
      if(!bone_heat(V,F,C,{},BE,{},W1))
      {
        return false;
      }
      break;
  }

  // Let L2 be the level appropriate for our skinning engine.
  if(computation_level == evaluation_level)
  {
    W = W1;
    return true;
  }
  // 4. Create evaluator from R1 with L2-L1.
  {
    assert(num_refined_vs == W1.rows());
    evaluation_eval = new_subdivision_evaluator( 
      num_refined_vs, (Scalar*)refined_vs, 
      num_refined_faces, (int*)refined_faces, 
      evaluation_level-computation_level);
  }

  // 5. Evaluate R1 using the weights as the control mesh positions to get weights W2 for L2.
  {
    const int num_weights = W1.cols();
    const int num_refined_vs = 
      num_refined_vertices_of_subdivision_evaluator( evaluation_eval );
    Scalar* refined_W = new Scalar[ num_refined_vs*num_weights ];
    const MatrixXS W1T = W1.cast<Scalar>().transpose();
    const Scalar * coarse_W = W1T.data();
    get_refined_vertices_of_subdivision_evaluator_with_control_vertices( 
      evaluation_eval,
      num_weights,
      coarse_W,
      refined_W);
    const auto & WT = Map<MatrixXS >(refined_W,num_weights,num_refined_vs);
    W = WT.cast<double>().transpose();
    delete[] refined_W;
  }


  delete[] refined_faces;
  delete[] refined_vs;
  delete_subdivision_evaluator( computation_eval );
  delete_subdivision_evaluator( evaluation_eval );
  return true;
}
