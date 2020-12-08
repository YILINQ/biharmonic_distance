#include "../include/biharmonic_distance.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/massmatrix.h"
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

void biharmonic_distance(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd &D)
{
  int N = V.rows();

  SparseMatrix<double> A, A_inv, Lc;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, A);
  igl::invert_diag(A, A_inv);
  igl::cotmatrix(V, F, Lc);

  SparseMatrix<double> Ld;
  Ld = A_inv * Lc;

  // compute Lc * A.inverse * Lc
  MatrixXd LcA_Lc = Lc * Ld;
  LcA_Lc.row(0).setZero();
  LcA_Lc.col(0).setZero();
  LcA_Lc(0, 0) = 1;

  // construct J matrix (from paper)
  MatrixXd J;
  J.resize(N, N);
  Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N);
  J = Eigen::MatrixXd::Identity(N, N) - (1.0 / N) * (Ones * Ones.transpose());
  J.row(0).setZero();

  // now LcA_Lc is invertible
  MatrixXd Gd;
  Gd = LcA_Lc.llt().solve(J);

  VectorXd off = (1.0 / N) * Gd.colwise().sum();
  MatrixXd offset = (off * Ones.transpose()).transpose();
  Gd -= offset;

  // dist D(i, j)^2 = Gd(i, i) + Gd(j, j) - 2 * Gd(i, j)
  VectorXd diag = Gd.diagonal();
  MatrixXd dd = diag * Ones.transpose();
  D.resize(N, N);
  D = sqrt((dd + dd.transpose() - 2 * Gd).array());
}
