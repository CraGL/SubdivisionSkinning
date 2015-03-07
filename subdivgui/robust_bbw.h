#ifndef ROBUST_BBW_H
#define ROBUST_BBW_H
#include <Eigen/Core>
bool robust_bbw(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  Eigen::MatrixXd & W);
#endif
