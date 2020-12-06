#include "../include/biharmonic_distance.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/pinv.h"
#include "igl/eigs.h"

using namespace Eigen;
using namespace std;
using namespace igl;



void biharmonic_distance(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd &D)
{
	std::string sep = "\n----------------------------------------\n";
	int N = V.rows();

	D.resize(N, N);
	SparseMatrix<double> A, A_inv, Lc;
	igl::massmatrix(V, F, MASSMATRIX_TYPE_VORONOI, A);
	igl::invert_diag(A, A_inv);
	igl::cotmatrix(V, F, Lc);


	SparseMatrix<double> Ld;
	Ld = A_inv * (-Lc);

	// compute Lc * A.inverse * Lc.inverse
	MatrixXd LcA_Lc_ = (-Lc) * Ld;
	
	LcA_Lc_.row(0).setZero();
	LcA_Lc_.col(0).setZero();
	// set 1st row and col to be 0
	// LcA_Lc_.prune([](int i, int j, double){
	// 	return i!=0 && j!=0;
	// });
	// set intersection to be 1
	LcA_Lc_(0, 0) = 1;

	// approximate approach
	// int K = 150; 

	// MatrixXd evecs;
	// VectorXd evals;
	// igl::eigs(Lc, A, K, EIGS_TYPE_SM, evecs, evals);
	// cout << evecs.rows() << " " << evecs.cols() << endl;	



	MatrixXd LcA_Lc_inverse = LcA_Lc_.inverse();
	// igl::pinv(LcA_Lc_, -1, LcA_Lc_inverse);
	// cout << LcA_Lc_ << sep;
	// make J matrix
	MatrixXd J;
	J.resize(N, N);
	Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N);
	J = Eigen::MatrixXd::Identity(N, N) - (1.0 / N) * (Ones * Ones.transpose());
	J.row(0).setZero();
	// SparseMatrix<double> JJ = J.sparseview();

	// // set 1st row of J to be 0
	// JJ.prune([](int i, int j, double){
	// 	return i!=0;
	// })

	// now LcA_Lc_ is invertible
	MatrixXd Gd;
	// Gd = LcA_Lc_.llt().solve(J);
	Gd = LcA_Lc_inverse * J;
	cout << Gd << sep;
	for(int j = 0; j < N; j++){
		
		Gd.col(j) -= (Ones.transpose() * Gd.col(j)) * Ones / (Ones.transpose()* Ones);
		
	}

	// dist D(i, j)^2 = Gd(i, i) + Gd(j, j) - 2*Gd(i, j)
	
	// set D;
	VectorXd diag = Gd.diagonal();
	MatrixXd dd = diag * Ones.transpose();
	D = sqrt((dd + dd.transpose() - 2*Gd).array());


}	