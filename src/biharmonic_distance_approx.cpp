#include "../include/biharmonic_distance_approx.h"
#include "../spectra-0.8.1/include/Spectra/MatOp/SparseCholesky.h"
#include "../spectra-0.8.1/include/Spectra/MatOp/SparseSymMatProd.h"
#include "../spectra-0.8.1/include/Spectra/SymGEigsSolver.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/massmatrix.h"
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>

using namespace Eigen;
using namespace Spectra;
using namespace std;

void biharmonic_distance_approx(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int k,
  Eigen::MatrixXd &D)
{
  int N = V.rows();
  D.resize(N, N);

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

  // approximate approach
  int K = k > N ? N : k;

  // in practice, setting convergence speed to be between 2K and 5K leads to the best result
  int conv_speed = K * 5;
  SparseSymMatProd<double> op(Lc);
  SparseCholesky<double> Bop(A);
  SymGEigsSolver<double, LARGEST_ALGE, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, K, conv_speed);

  geigs.init();
  int nconv = geigs.compute(10000);
  VectorXd evalues;
  MatrixXd evectors;

  if (geigs.info() == SUCCESSFUL) {
    evalues = geigs.eigenvalues();
    evectors = geigs.eigenvectors();
  } else {
    cout<<"Failed to solve the generalized eigenvalue problem.\n";
    return;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double d = 0.0;
      for (int a = 0; a < K; a++) {
        d += 1.0 / (evalues(a) * evalues(a)) * (evectors(i, a) - evectors(j, a)) * (evectors(i, a) - evectors(j, a));
      }
      D(i, j) = sqrt(d);
    }
  }
  return;
}
