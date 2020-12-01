#include "../include/signing_the_unsigned.h"
#include "igl/knn.h"
#include "igl/octree.h"



using namespace Eigen;
using namespace std;

void signing_the_unsigned(
    const Eigen::MatrixXd & P,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
	// code here

	// compute K-nn
	// K ranges from 12 - 30
	int K = 19;
	vector<vector<int> > O_PI;
	MatrixXi O_CH;
	MatrixXd O_CN;
	VectorXd O_W;
	MatrixXi KNN;
	igl::octree(P, O_PI, O_CH, O_CN, O_W);
	igl::knn(P, K, O_PI, O_CH, O_CN, O_W, KNN);

	// compute unsigned distance
	VectorXd dU;
	dU.resize(P.rows());
	for(int i = 0; i < P.rows(); i++){
		double u = 0;
		for(int j = 0; j < KNN.cols(); j++){
			u += (P.row(i) - P.row(KNN(i, j))).norm();
		}
		dU(i) = 1.0 / K * sqrt(u);

	}
}