#include "../include/biharmonic_distance.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/pinv.h"
#include "igl/eigs.h"

#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include <../spectra-0.8.1/include/Spectra/GenEigsSolver.h>
#include <../spectra-0.8.1/include/Spectra/MatOp/SparseGenMatProd.h>


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
	Eigen::SparseMatrix<double> A, A_inv, Lc;
	igl::massmatrix(V, F, MASSMATRIX_TYPE_BARYCENTRIC, A);
	igl::invert_diag(A, A_inv);
	igl::cotmatrix(V, F, Lc);


	Eigen::SparseMatrix<double> Ld;
	Ld = A_inv * (Lc);

	// compute Lc * A.inverse * Lc
	Eigen::MatrixXd LcA_Lc_ = (Lc) * Ld;
	
	LcA_Lc_.row(0).setZero();
	LcA_Lc_.col(0).setZero();

	LcA_Lc_(0, 0) = 1;

	// approximate approach
	int K = 200;
	if (K > N){
		K = N;
	}
	cout << "1\n";
	SparseGenMatProd<double> op(Ld);
	GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, 3, 6);
	cout << "11\n";

	eigs.init();
	cout << "111\n";

    int nconv = eigs.compute();
    Eigen::VectorXcd evalues;

    if(eigs.info() == SUCCESSFUL)
        evalues = eigs.eigenvalues();

   	std::cout << "Eigenvalues found:\n" << evalues << std::endl;
   	return;

	// make J matrix
	Eigen::MatrixXd J;
	J.resize(N, N);
	Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N);
	J = Eigen::MatrixXd::Identity(N, N) - (1.0 / N) * (Ones * Ones.transpose());
	J.row(0).setZero();


	// now LcA_Lc_ is invertible
	Eigen::MatrixXd Gd;
	Gd = LcA_Lc_.llt().solve(J);

	Eigen::VectorXd off = (1.0 / N) * Gd.colwise().sum();
	Eigen::MatrixXd offset = (off * Ones.transpose()).transpose();
	Gd -= offset;

	// dist D(i, j)^2 = Gd(i, i) + Gd(j, j) - 2*Gd(i, j)
	// set D;
	Eigen::VectorXd diag = Gd.diagonal();
	Eigen::MatrixXd dd = diag * Ones.transpose();
	D = sqrt((dd + dd.transpose() - 2*Gd).array());

}	