// Given a surface mesh (V,F), compute the (approximated) biharmonic distance
// between all pairs of vertices.
//
// Inputs:
//   V  #V by 3 list of vertex positions
//   F  #F by 3 list of triangle indices into the rows of V
//   k  number of eigenvectors to use for computing the approximated distance
// Outputs:
//   D  symmetric #V by #V matrix of distances where the entry at (i,j) is the
//      approximated biharmonic distance between vertex i and j

#include <Eigen/Core>

void biharmonic_distance_approx(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int k,
  Eigen::MatrixXd &D);