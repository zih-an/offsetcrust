#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//Eigen::VectorXd solve_ldlt(const Eigen::MatrixXd& CoeffMat,
//	const Eigen::VectorXd& right,
//	double pinvtoler);
Eigen::VectorXd solve_ldlt(const Eigen::SparseMatrix<double>& CoeffMat,
	const Eigen::VectorXd& right,
	double pinvtoler);
Eigen::VectorXd solve_ldlt(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& ldlt,
	const Eigen::VectorXd& right,
	double pinvtoler);