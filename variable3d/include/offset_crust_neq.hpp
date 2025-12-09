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
#include <igl/cotmatrix.h>
#include <CGAL/assertions.h>
#include <igl/grad.h>
#include <igl/blue_noise.h>
#include <igl/barycentric_interpolation.h>
#include <igl/doublearea.h>
#include <igl/heat_geodesics.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <igl/dual_contouring.h>




namespace Offset3D {
	struct Source {
		int id;
		double dis;

		Source(){}
		Source(int _id, double _dis) :id(_id), dis(_dis) { }
	};

	class OffsetCrust {
		bool VERBOSE_ON = true;
		MetaData VERBOSE;
		double diagonal_length;
		ProjPoint* proj_model;

	private:
		Tree tree;
		std::list<Triangle> triangles;


	protected:
		//sphere center=0,0,0; r=1
		Eigen::MatrixXd SPHV, SPHVN;
		Eigen::MatrixXi SPHF;
		const double RADIUS = 1;


	protected:
		Eigen::MatrixXd OrgVertices, OrgFaceNormals;  //cols=3
		Eigen::MatrixXi OrgFaces;
		double q_offset = 1e-6,
			OFFSET;
		bool __IS_INSIDE, __outputPoly, __useCentroid, __relativedis;

		Eigen::MatrixXd PrePVerts, PreQVerts;
		std::map<int, int> qi2pi;
		std::map<int, std::vector<int>> pi2qi;
		std::map<int, double> pi2dis;
		int PSize;
		Eigen::AlignedBox<double, 3> bounding_box;


		/*parameters*/
		int	num_corner_points = 5;
		double
			diheralBar,
			lambda,
			ratio;
		int blue_noise_desired, seg_angle;


		CGAL::Bbox_3 bbox = CGAL::Bbox_3(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5);



	private:
		void init_tree();

		void compute_power_diagram(
			Eigen::MatrixXd& AllVertices,
			Eigen::MatrixXd& QNormals,
			Eigen::VectorXi& blue_validFlags,
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




	protected:
		std::string __filename;

		Eigen::VectorXd __VDIS;
		Eigen::MatrixXd gradient;

		std::tuple<double, double, double> tri_interpolation(
			Eigen::Vector3d& v0,
			Eigen::Vector3d& v1,
			Eigen::Vector3d& v2,
			int fid,
			Eigen::Vector3d& p
		);
		double calc_tri_area(
			Eigen::Vector3d& v1,
			Eigen::Vector3d& v2,
			Eigen::Vector3d& v3
		);

		double get_F(Eigen::VectorXd& pt, int fid);
		bool get_VDIS(
			std::string outOrgFile, 
			std::string colorbar,
			std::vector<Source>& sources);
		bool get_gallery_VDIS(
			std::string outOrgFile, 
			std::string colorbar
		);
		bool get_VDIS_geodesic(
			std::string outOrgFile,
			std::string colorbar,
			std::vector<Source>& sources
		);


		void get_normals(
			Eigen::MatrixXd& FaceNormals
		);


	public:
		OffsetCrust();
		OffsetCrust(
			std::string filename,
			int _num_angle = 10,
			int blue_points = 70000,
			double corner_ratio = 0.05,
			double _lambda = 0.01
		);
		~OffsetCrust() {
			delete proj_model;
		}

		void core_main(
			std::vector<Source>& sources,
			bool inside,
			std::string outOBJFile,
			std::string outMetaFile,
			bool relative_dis,
			bool outPoly,
			bool useCentroid,
			bool do_sharp,
			bool runGalleryDis,
			double diheral
		);

		void core_main(
			Eigen::MatrixXd& V,
			Eigen::MatrixXd& NV,
			Eigen::VectorXd& Radius,
			Eigen::MatrixXd& Grad,
			std::string filename,
			bool isConstant = false
		);

		void core_main(
			Eigen::MatrixXd& P,
			std::vector<std::vector<Eigen::VectorXd>>& PNormals,
			Eigen::VectorXd& Radius,
			std::string filename
		);

		static void rotateN(
			std::vector<Eigen::VectorXd>& normals,
			Eigen::Vector3d& grad
		);

	};


}

