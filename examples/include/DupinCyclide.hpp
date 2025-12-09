#pragma once

#include <offset_crust_neq.hpp>
#include <igl/sparse_voxel_grid.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <cmath>
#include <random>

using namespace Offset3D;





class DupinCyclide {

private:
	Eigen::MatrixXd V, VN, RGrad;
	Eigen::VectorXd Radius;
	Eigen::MatrixXi F;

	bool isConstant = false;

	double 
		a = 1, 
		b = 0.98, 
		c = 0.199, 
		d = 0.3; 
	double R_max = d + c;
	double R_min = d - c;

	double implicit_func(Eigen::Vector3d& pt);
	Eigen::Vector3d parametric_func(double u, double v);


	double offset_radius(double u, double v);
	double pt2surfDis(Eigen::Vector3d& pt);
	Eigen::Vector3d grad(Eigen::Vector3d& pt);
	Eigen::Vector3d surface_gradient(double u, double v);
	std::pair<double, double> estimate_initial_uv(const Eigen::Vector3d& p);

	void sample(
		Eigen::MatrixXd& outV,
		Eigen::MatrixXd& outVN,
		Eigen::VectorXd& outR,
		Eigen::MatrixXd& outRGrad,
		const int n_desired = 100
	);

	double get_F(
		Eigen::VectorXd& pt, 
		int fid, 
		Eigen::Vector3d& grad
	);

	
	void random_sample_para(int n);


public:
	DupinCyclide() {};

	void mc(int res = 300);
	void para(int u_steps = 13, int v_steps=11);

	void smooth_main();
	void mesh_main();
	void mat(int u_steps = 10, int v_steps = 10);

};

