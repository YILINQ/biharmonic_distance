#include "../include/biharmonic_distance.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/pinv.h"
#include "igl/eigs.h"

#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace Spectra;
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
	igl::massmatrix(V, F, MASSMATRIX_TYPE_BARYCENTRIC, A);
	igl::invert_diag(A, A_inv);
	igl::cotmatrix(V, F, Lc);


	SparseMatrix<double> Ld;
	Ld = A_inv * (Lc);
	// cout << "A_inv \n" << A_inv <<sep;
	// cout << "Lc \n" << Lc <<sep;

	// compute Lc * A.inverse * Lc.inverse
	MatrixXd LcA_Lc_ = (Lc) * Ld;
	
	LcA_Lc_.row(0).setZero();
	LcA_Lc_.col(0).setZero();

	LcA_Lc_(0, 0) = 1;

	// TODO approximate approach
	int K = 200;
	if(N < K){
		K = N;
	}

	
	// approximate approach


	// make J matrix
	MatrixXd J;
	J.resize(N, N);
	Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N);
	J = Eigen::MatrixXd::Identity(N, N) - (1.0 / N) * (Ones * Ones.transpose());
	J.row(0).setZero();


	// now LcA_Lc_ is invertible
	MatrixXd Gd;
	// cout << LcA_Lc_ << endl;
	Gd = LcA_Lc_.llt().solve(J);
	// cout << Gd << sep;

	VectorXd off = (1.0 / N) * Gd.colwise().sum();
	MatrixXd offset = (off * Ones.transpose()).transpose();
	Gd -= offset;
	// for(int j = 0; j < N; j++){
	// 	Gd.col(j) -= (1.0 / N) * (Gd.col(j).sum() * Ones);
	// }

	// dist D(i, j)^2 = Gd(i, i) + Gd(j, j) - 2*Gd(i, j)
	
	// set D;
	VectorXd diag = Gd.diagonal();
	MatrixXd dd = diag * Ones.transpose();
	D = sqrt((dd + dd.transpose() - 2*Gd).array());

}	