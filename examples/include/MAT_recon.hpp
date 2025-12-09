#pragma once

#include <offset_crust_neq.hpp>
#include <igl/sparse_voxel_grid.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <cmath>
#include <random>


namespace Offset3D {
	struct PInfo {
		Eigen::VectorXd pt;
		double pr;
		std::vector<Eigen::VectorXd> qt;
		std::vector<Eigen::VectorXd> qn;
		double q_offset = 1e-6;

		PInfo() {}
		PInfo(const Eigen::VectorXd& p, double r)
			:pt(p), pr(r) {}
		PInfo(const Eigen::VectorXd& p, double r, const Eigen::VectorXd& q, const Eigen::VectorXd& n)
			:pt(p), pr(r)
		{
			qt.push_back(q);
			qn.push_back(n);
		}
		void insert(const Eigen::VectorXd& q, const Eigen::VectorXd& n) {
			qt.push_back(q);
			qn.push_back(n);
		}
		void refresh_qt() {
			for (int i = 0; i < qt.size(); i++) {
				qt[i] = pt + q_offset * qn[i];
			}
		}

	};


	class MATRecon : public OffsetCrust {

		std::vector<Eigen::Vector3d> MA_v;
		std::vector<double> MA_r;
		std::vector<std::vector<int>> MA_f, MA_e;
		Eigen::MatrixXd maV;
		Eigen::MatrixXi maF, maE;


		int blue_noise_desired = 50000;
		int seg_angle = 10;
		double blue_noise_radius;
		double eps = 1e-6;

		bool VERBOSE_ON = true;
		MetaData VERBOSE;

		void read_mat(std::string matfile);

		void preprocessing_MA(
			Eigen::MatrixXd& maV,
			Eigen::MatrixXi& maF,
			Eigen::MatrixXi& maE,
			Eigen::VectorXd VDIS,
			Eigen::MatrixXd& outP,
			std::vector<std::vector<Eigen::VectorXd>>& outPNormals,
			Eigen::VectorXd& outRadius
		);

		bool sampleTriangleBand(
			double blue_noise_radius,
			int fid,
			std::vector<PInfo>& output,
			double q_orient,
			std::vector<int>& blueids,
			Eigen::MatrixXd& P_blue,
			Eigen::VectorXi& blue_validFlags
		);



	public:
		MATRecon(std::string filename);


	};

}


