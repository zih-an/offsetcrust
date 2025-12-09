#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <cgal_types_3d.hpp>
#include <face.hpp>
#include <map>
#include <vector>
#include <stack>
#include <queue>
#include <Eigen/Core>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/random_points_on_mesh.h>
#include <utils.hpp>
#include <OutputOBJ.hpp>
#include <igl/winding_number.h>
#include <igl/writeOFF.h>
#include <limits>
#include <igl/blue_noise.h>
#include <proj_point.hpp>
#include <edgePool.hpp>
#include <random>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/resolve_duplicated_faces.h>
#include <CGAL/assertions.h>
#include <igl/blue_noise.h>
#include <igl/barycentric_interpolation.h>
#include <igl/doublearea.h>



namespace Offset3D {
	class OffsetCrust {
		bool VERBOSE_ON = true;
		MetaData VERBOSE;
		double diagonal_length;

	private:
		Tree tree;
		std::list<Triangle> triangles;


	private:
		Eigen::MatrixXd SPHV, SPHVN;
		Eigen::MatrixXi SPHF;
		const double RADIUS = 1;

		ProjPoint* PQP;


	private:
		Eigen::MatrixXd OrgVertices, OrgFaceNormals;  //cols=3
		Eigen::MatrixXi OrgFaces;
		double q_offset = 1e-6,
			OFFSET;
		bool __IS_INSIDE, __outputPoly, __useCentroid;

		int	num_corner_points = 5;


		Eigen::MatrixXd PrePVerts, PreQVerts;
		std::map<int, int> qi2pi;
		std::map<int, std::vector<int>> pi2qi;
		int PSize;


		/*parameters*/
		double
			diheralBar,
			lambda,
			ratio;
		int blue_noise_desired, seg_angle;


		CGAL::Bbox_3 bbox = CGAL::Bbox_3(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5);



	private:
		void compute_power_diagram(
			Eigen::MatrixXd& AllVertices,
			Eigen::MatrixXd& QNormals,
			Eigen::VectorXi& blue_validFlags,
			double w_diff,
			std::string vorFilename
		);

		void extract_voronoi(
			Regular_triangulation& rt,
			Eigen::MatrixXd& QNormals,
			std::vector<VorFaceTri*>& vor_diagram //output (voronoi)
		);


		void optim_vor_diagram(
			std::vector<VorFaceTri*>& vor_diagram,
			Eigen::MatrixXd& QNormals,
			Eigen::MatrixXd& AllVertices,
			std::string vorFilename
		);



	private:
		void preprocessing_with_sphere(
			Eigen::MatrixXd& Vertices,
			Eigen::MatrixXi& Faces,
			Eigen::MatrixXd& OutAllVertices,
			Eigen::MatrixXd& OutQNormals,
			Eigen::VectorXi& blue_validFlags,
			double q_orient = 1
		);

		void preprocessing_sharp(
			Eigen::MatrixXd& Vertices,
			Eigen::MatrixXi& Faces,
			Eigen::MatrixXd& OutAllVertices,
			Eigen::MatrixXd& OutQNormals,
			Eigen::VectorXi& blue_validFlags,
			double q_orient = 1
		);


		void generateWeights(
			Eigen::MatrixXd& AllVertices,
			Eigen::VectorXd& AllWeights, //output
			double d_diff
		);

		void QEM_K(
			Eigen::VectorXd& normal,
			double d,
			Eigen::MatrixXd& K
		);

		void postproc_remove_dup(
			Eigen::MatrixXd& polyV,
			std::vector<std::vector<int>>& polyF,
			Eigen::MatrixXd& newV,
			Eigen::MatrixXi& newF
		);



	public:
		OffsetCrust(
			std::string filename,
			int _num_angle = 10,
			int blue_points = 70000,
			double corner_ratio = 0.05,
			double _lambda = 0.01
		);
		~OffsetCrust() {
			delete PQP;
		}
		void core_main(
			double d,
			bool is_inside,
			std::string outOBJFile,
			std::string outMetaFile,
			bool relative_dis,
			bool outPoly,
			bool useCentroid,
			bool do_sharp,
			double diheral
		);

	};


}

