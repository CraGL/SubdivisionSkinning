#include <Eigen/Core>
enum WeightsType
{
  WEIGHTS_TYPE_BBW = 0,
  WEIGHTS_TYPE_BONE_HEAT = 1,
  NUM_WEIGHTS_TYPES = 2
};
bool subdiv_weights(
  const Eigen::MatrixXd & CV,
  const Eigen::MatrixXi & CQ,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  const int computation_level,
  const int evaluation_level,
  const WeightsType type,
  Eigen::MatrixXd & W);
  
