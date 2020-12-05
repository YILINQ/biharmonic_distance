#include "../include/biharmonic_distance.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/pinv.h"

using namespace Eigen;
using namespace std;
using namespace igl;



void biharmonic_distance(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd &D)
{
	int N = V.rows();

	D.resize(N, N);
	SparseMatrix<double> A, A_inv, Lc;
	igl::massmatrix(V, F, MASSMATRIX_TYPE_VORONOI, A);
	igl::invert_diag(A, A_inv);
	igl::cotmatrix(V, F, Lc);

	SparseMatrix<double> Ld;
	Ld = A_inv * Lc;

	// compute Lc * A.inverse * Lc.inverse
	MatrixXd LcA_Lc_ = Lc * Ld;
	LcA_Lc_.row(0) *= 0;
	LcA_Lc_.col(0) *= 0;
	// set 1st row and col to be 0
	// LcA_Lc_.prune([](int i, int j, double){
	// 	return i!=0 && j!=0;
	// });
	// set intersection to be 1
	LcA_Lc_(0, 0) = 1;

	MatrixXd LcA_Lc_inverse = LcA_Lc_.inverse();

	// make J matrix
	MatrixXd J;
	J.resize(N, N);
	Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N);
	J = Eigen::MatrixXd::Identity(N, N) - (1.0 / N) * (Ones * Ones.transpose());
	J.row(0) *= 0;
	// SparseMatrix<double> JJ = J.sparseview();

	// // set 1st row of J to be 0
	// JJ.prune([](int i, int j, double){
	// 	return i!=0;
	// })

	// now LcA_Lc_ is invertible
	MatrixXd Gd;
	Gd.resize(N, N);

	for(int j = 0; j < N; j++){
		Gd.col(j) = LcA_Lc_inverse * J.col(j);
		Gd.col(j) -= (Ones.transpose() * Gd.col(j)) * Ones / (Ones.transpose()* Ones);
		// Gd.col(j) -= Gd.col(j).sum() * Ones;
	}

	// dist D(i, j)^2 = Gd(i, i) + Gd(j, j) - 2*Gd(i, j)
	// D = -2 * Gd;
	
	// set D;
	VectorXd diag = Gd.diagonal();
	MatrixXd dd = diag * Ones.transpose();

	D = sqrt((dd + dd.transpose() - 2*Gd).array());

	// for(int i = 0; i < N; i++){
	// 	for(int j = 0; j < N; j++){
	// 		D(i, j) = sqrt(Gd(i, i) + Gd(j, j) - 2*Gd(i, j));
	// 	}
	// }



}	