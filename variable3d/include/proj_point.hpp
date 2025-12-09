#pragma once

#include <cgal_types_3d.hpp>
#include <Eigen/Core>
#include <PQP.h>


namespace Offset3D {
	class ProjPoint {
	private:
		PQP_Model m;
		Eigen::MatrixXd FaceNormals, __V;
		Eigen::MatrixXi __F;
		bool __IS_INSIDE;

	public:
		ProjPoint(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& FN, bool flag) {
			FaceNormals = FN;
			__F = F;
			__V = V;

			m.BeginModel();
			for (int f_id = 0; f_id < F.rows(); f_id++) {
				int v0 = F(f_id, 0),
					v1 = F(f_id, 1),
					v2 = F(f_id, 2);

				PQP_REAL tri[3][3];
				tri[0][0] = V(v0, 0), tri[0][1] = V(v0, 1), tri[0][2] = V(v0, 2);
				tri[1][0] = V(v1, 0), tri[1][1] = V(v1, 1), tri[1][2] = V(v1, 2);
				tri[2][0] = V(v2, 0), tri[2][1] = V(v2, 1), tri[2][2] = V(v2, 2);
				m.AddTri(tri[0], tri[1], tri[2], f_id);
			}
			m.EndModel();
			__IS_INSIDE = flag;
		}

		ProjPoint(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
			__F = F;
			__V = V;

			m.BeginModel();
			for (int f_id = 0; f_id < F.rows(); f_id++) {
				int v0 = F(f_id, 0),
					v1 = F(f_id, 1),
					v2 = F(f_id, 2);

				PQP_REAL tri[3][3];
				tri[0][0] = V(v0, 0), tri[0][1] = V(v0, 1), tri[0][2] = V(v0, 2);
				tri[1][0] = V(v1, 0), tri[1][1] = V(v1, 1), tri[1][2] = V(v1, 2);
				tri[2][0] = V(v2, 0), tri[2][1] = V(v2, 1), tri[2][2] = V(v2, 2);
				m.AddTri(tri[0], tri[1], tri[2], f_id);
			}
			m.EndModel();
		}

		//query point
		double querypt(
			Eigen::VectorXd& query_ptV,
			Eigen::VectorXd& pt,
			Eigen::VectorXd& normal,
			bool& onface
		) {
			PQP_DistanceResult res;
			double point[3] = {
				query_ptV(0), query_ptV(1), query_ptV(2)
			};

			PQP_Distance(&res, &m, point, 0.0, 0.0);

			if (res.pos_flag == 0) onface = true;
			else onface = false;

			int fid = res.last_tri->id;

			normal = FaceNormals.row(fid);
			Eigen::VectorXd face_v0 = __V.row(__F(fid, 0));
			

			Eigen::VectorXd proj_p1(3), proj_p2(3);
			proj_p1 << res.P1()[0], res.P1()[1], res.P1()[2];
			pt = proj_p1;
			
			if (__IS_INSIDE) normal = -normal;
			return res.Distance();

		}

		double getdis(Eigen::VectorXd& query_ptV) {
			PQP_DistanceResult res;
			double point[3] = {
				query_ptV(0), query_ptV(1), query_ptV(2)
			};

			PQP_Distance(&res, &m, point, 0.0, 0.0);
			return res.distance;
		}

		double getdis(
			Eigen::VectorXd& query_ptV,
			Eigen::VectorXd& pt
		) {
			PQP_DistanceResult res;
			double point[3] = {
				query_ptV(0), query_ptV(1), query_ptV(2)
			};

			PQP_Distance(&res, &m, point, 0.0, 0.0);

			Eigen::VectorXd proj_p1(3);
			proj_p1 << res.P1()[0], res.P1()[1], res.P1()[2];
			pt = proj_p1;

			return res.distance;
		}


		bool  querypt_interpolate(
			Eigen::VectorXd& query_ptV,
			Eigen::VectorXd& base_pt,
			Eigen::VectorXd& normal,  //output
			double epsilon
		) {
			PQP_DistanceResult res;
			double point[3] = {
				query_ptV(0), query_ptV(1), query_ptV(2)
			};

			PQP_Distance(&res, &m, point, 0.0, 0.0);

			int fid = res.last_tri->id;
			normal = FaceNormals.row(fid);
			if (__IS_INSIDE) normal = -normal;

			if (res.pos_flag >= 4 && res.pos_flag <= 6) {
				int vid = res.pos_flag - 4;
				Eigen::VectorXd face_v = __V.row(__F(fid, vid));
				if ((base_pt - face_v).norm() <= 1e-7) {
					return true;
				}
			}
			

			return false;
		}
		
	};


}
