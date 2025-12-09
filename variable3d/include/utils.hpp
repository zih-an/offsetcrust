#pragma once

#include <Eigen/Core>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <nanoflann.hpp>
#include <nanoflann_util.h>
#include <igl/random_points_on_mesh.h>
#include <igl/per_face_normals.h>



class UTILS {
public:
	static std::string data_dir;
	static std::string input_dir, out_dir, sphere_file, colorbar_file;

	static void normalization(Eigen::MatrixXd& M) {
		M = (M - M.minCoeff() * Eigen::MatrixXd::Ones(M.rows(), M.cols())) / (M.maxCoeff() - M.minCoeff());
	}
	static void normalization_2(Eigen::MatrixXd& M) {
		M = (M - M.minCoeff() * Eigen::MatrixXd::Ones(M.rows(), M.cols())) / (M.maxCoeff() - M.minCoeff());
		M = M * 2 - Eigen::MatrixXd::Ones(M.rows(), M.cols());  // (-1,1)
	}

	static void normalize_centroid(Eigen::MatrixXd& vertices) {
		Eigen::RowVector3d centroid = vertices.colwise().mean();

		vertices.rowwise() -= centroid;
		double maxDistance = vertices.rowwise().norm().maxCoeff();
		if (maxDistance > 0) {
			vertices /= maxDistance;
		}
	}

	static bool isPointInTriangle(
		const Eigen::VectorXd& P,
		const Eigen::VectorXd& A,
		const Eigen::VectorXd& B,
		const Eigen::VectorXd& C)
	{
		Eigen::Vector3d v0 = C - A;
		Eigen::Vector3d v1 = B - A;
		Eigen::Vector3d v2 = P - A;

		double dot00 = v0.dot(v0);
		double dot01 = v0.dot(v1);
		double dot02 = v0.dot(v2);
		double dot11 = v1.dot(v1);
		double dot12 = v1.dot(v2);

		double inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

		double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
		if (u < 0 || u > 1) // if u out of range, return directly
		{
			return false;
		}

		double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
		if (v < 0 || v > 1) // if v out of range, return directly
		{
			return false;
		}

		return u + v <= 1;
	}

	static void outputXYZ(
		Eigen::MatrixXd& Vertices,
		Eigen::MatrixXd& Normals,
		std::string filename
	) {
		assert(Vertices.rows() == Normals.rows());
		if (Vertices.rows() != Normals.rows())
		{
			std::cout << "size error" << std::endl;
			return;
		}

		std::cout << "output xyz..." << std::endl;
		std::ofstream file(UTILS::out_dir + filename + "_normal.xyz");
		for (int i = 0; i < Vertices.rows(); i++) {
			file <<
				std::to_string(Vertices(i, 0)) + " "
				+ std::to_string(Vertices(i, 1)) + " "
				+ std::to_string(Vertices(i, 2)) + " "
				+ std::to_string(Normals(i, 0)) + " "
				+ std::to_string(Normals(i, 1)) + " "
				+ std::to_string(Normals(i, 2)) + "\n";
		}
		file.close();
		std::cout << UTILS::out_dir + filename + "_normal.xyz " << "output done." << std::endl;
	}

	static void outputOBJ_with_normal(
		Eigen::MatrixXd& Vertices,
		Eigen::MatrixXd& Normals,
		std::string filename,
		double len = 0.05
	) {
		assert(Vertices.rows() == Normals.rows());
		std::cout << "output xyz..." << std::endl;
		std::ofstream file(UTILS::data_dir + filename + "_line.obj");
		int cnt = 0;
		for (int i = 0; i < Vertices.rows(); i++) {
			Eigen::VectorXd
				v1 = Vertices.row(i),
				v2 = Vertices.row(i),
				n = Normals.row(i);
			v2 += len * n;
			file <<
				"v " + std::to_string(v1(0)) + " " + std::to_string(v1(1)) + " " + std::to_string(v1(2)) + "\n";
			file <<
				"v " + std::to_string(v2(0)) + " " + std::to_string(v2(1)) + " " + std::to_string(v2(2)) + "\n";
			file <<
				"l " + std::to_string(cnt + 1) + " " + std::to_string(cnt + 2) + "\n";
			cnt += 2;
		}
		file.close();


		std::ofstream filev(UTILS::data_dir + filename + "_verts.obj");
		for (int i = 0; i < Vertices.rows(); i++) {
			filev <<
				"v "
				+ std::to_string(Vertices(i, 0)) + " "
				+ std::to_string(Vertices(i, 1)) + " "
				+ std::to_string(Vertices(i, 2)) + "\n";
		}
		filev.close();

	}

	static void outputOBJ_with_lines(
		Eigen::MatrixXd& Vertices,
		Eigen::MatrixXd& Lines,
		std::string filename
	) {
		std::cout << "output xyz..." << std::endl;
		std::ofstream file(UTILS::data_dir + filename + "_line.obj");
		int cnt = 0;
		for (int i = 0; i < Lines.rows() - 1; i++) {
			Eigen::VectorXd
				v1 = Lines.row(i),
				v2 = Lines.row(i + 1);
			file <<
				"v " + std::to_string(v1(0)) + " " + std::to_string(v1(1)) + " " + std::to_string(v1(2)) + "\n";
			file <<
				"v " + std::to_string(v2(0)) + " " + std::to_string(v2(1)) + " " + std::to_string(v2(2)) + "\n";
			file <<
				"l " + std::to_string(cnt + 1) + " " + std::to_string(cnt + 2) + "\n";
			cnt += 2;
		}
		file.close();


		std::ofstream filev(UTILS::data_dir + filename + "_verts.obj");
		for (int i = 0; i < Vertices.rows(); i++) {
			filev <<
				"v "
				+ std::to_string(Vertices(i, 0)) + " "
				+ std::to_string(Vertices(i, 1)) + " "
				+ std::to_string(Vertices(i, 2)) + "\n";
		}
		filev.close();

	}

	static void outputOBJ_with_lines(
		Eigen::MatrixXd& Lines,
		std::string filename
	) {
		std::cout << "output xyz..." << std::endl;
		std::ofstream file(UTILS::data_dir + filename + "_line.obj");
		int cnt = 0;
		for (int i = 0; i < Lines.rows() - 1; i += 2) {
			Eigen::VectorXd
				v1 = Lines.row(i),
				v2 = Lines.row(i + 1);
			file <<
				"v " + std::to_string(v1(0)) + " " + std::to_string(v1(1)) + " " + std::to_string(v1(2)) + "\n";
			file <<
				"v " + std::to_string(v2(0)) + " " + std::to_string(v2(1)) + " " + std::to_string(v2(2)) + "\n";
			file <<
				"l " + std::to_string(cnt + 1) + " " + std::to_string(cnt + 2) + "\n";
			cnt += 2;
		}
		file.close();
	}

	static void outputXYZ(
		Eigen::MatrixXd& Vertices,
		std::string filename
	) {
		std::cout << "output xyz... " << UTILS::out_dir + filename + "_verts.xyz" << std::endl;
		std::ofstream file(UTILS::out_dir + filename + "_verts.xyz");
		for (int i = 0; i < Vertices.rows(); i++) {
			file <<
				std::to_string(Vertices(i, 0)) + " "
				+ std::to_string(Vertices(i, 1)) + " "
				+ std::to_string(Vertices(i, 2)) + "\n";
		}
		file.close();
	}

	static void random_sample(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::MatrixXd& outV,
		Eigen::MatrixXd& outVN,
		const int n_desired
	) {
		Eigen::MatrixXd FN;
		igl::per_face_normals(V, F, FN);

		// Heuristic to  determine radius from desired number 
		Eigen::MatrixXd B;
		Eigen::VectorXi FI;
		Eigen::MatrixXd P_white;
		igl::random_points_on_mesh(n_desired, V, F, B, FI, outV);
		outVN.conservativeResize(outV.rows(), outV.cols());
		for (int i = 0; i < FI.size(); i++) {
			int id = FI(i);
			outVN.row(i) = FN.row(id);
		}
	}

	static double calc_tri_area(
		Eigen::Vector3d& v1,
		Eigen::Vector3d& v2,
		Eigen::Vector3d& v3
	) {
		Eigen::Vector3d AB = v2 - v1;
		Eigen::Vector3d AC = v3 - v1;
		return 0.5 * AB.cross(AC).norm();
	}


	static std::tuple<double, double, double> tri_interpolation(
		Eigen::Vector3d& v0,
		Eigen::Vector3d& v1,
		Eigen::Vector3d& v2,
		int fid,
		Eigen::Vector3d& p
	) {
		double area_ABC = UTILS::calc_tri_area(v0, v1, v2);
		double area_PBC = UTILS::calc_tri_area(p, v1, v2),
			area_PCA = UTILS::calc_tri_area(p, v2, v0),
			area_PAB = UTILS::calc_tri_area(p, v0, v1);

		// bc
		double alpha = area_PBC / area_ABC;
		double beta = area_PCA / area_ABC;
		double gamma = area_PAB / area_ABC;

		return std::make_tuple(alpha, beta, gamma);
	}


	static Eigen::Vector3d tri_grad(
		Eigen::Vector3d& vv0,
		Eigen::Vector3d& vv1,
		Eigen::Vector3d& vv2,
		double val0, double val1, double val2
	) {
		Eigen::Vector3d
			e_ik = vv0 - vv2,
			e_ji = vv1 - vv0,
			n = e_ik.cross(e_ji);
		double S2 = n.norm(); //OFacesArea[fid] * 2
		n.normalize();
		Eigen::Vector3d
			vec_ik = n.cross(e_ik),
			vec_ji = n.cross(e_ji);
		double
			d_ji = val1 - val0,
			d_ki = val2 - val0;

		Eigen::Vector3d tmp = (d_ji / S2) * vec_ik + (d_ki / S2) * vec_ji;

		return tmp;
	}


	static void rand_double(double r1, double r2, Eigen::MatrixXd& samples, int n) {
		assert(samples.rows() == n);
		assert(samples.cols() == 1);
		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<double> dis(r1, r2);

		for (int i = 0; i < n; ++i) {
			samples(i) = dis(gen);
		}
	}




	template <class MatT>
	static Eigen::Matrix<typename MatT::Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime>
		pseudoinverse(const MatT& mat, typename MatT::Scalar tolerance = typename MatT::Scalar{ 1e-4 }) // choose appropriately
	{
		typedef typename MatT::Scalar Scalar;
		auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
		const auto& singularValues = svd.singularValues();
		Eigen::Matrix<Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
		singularValuesInv.setZero();
		for (unsigned int i = 0; i < singularValues.size(); ++i) {
			if (singularValues(i) > tolerance)
			{
				singularValuesInv(i, i) = Scalar{ 1 } / singularValues(i);
			}
			else
			{
				singularValuesInv(i, i) = Scalar{ 0 };
			}
		}
		return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
	}




	static void write_obj_with_texture(
		const std::string& filename,
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXd& UV,
		const std::string& mtl
	) {
		std::ofstream obj_file(filename);

		if (!obj_file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return;
		}

		size_t pos = mtl.find_last_of("/\\");
		std::string mtl_filename = (pos == std::string::npos) ? mtl : mtl.substr(pos + 1);
		obj_file << "mtllib " << mtl_filename << std::endl;

		for (int i = 0; i < V.rows(); ++i) {
			obj_file << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
		}

		for (int i = 0; i < UV.rows(); ++i) {
			obj_file << "vt " << UV(i, 0) << " " << UV(i, 1) << std::endl;
		}

		obj_file << "usemtl mat_texture" << std::endl;

		for (int i = 0; i < F.rows(); ++i) {
			obj_file << "f "
				<< F(i, 0) + 1 << "/" << F(i, 0) + 1 << " " 
				<< F(i, 1) + 1 << "/" << F(i, 1) + 1 << " "
				<< F(i, 2) + 1 << "/" << F(i, 2) + 1 << std::endl;
		}

		obj_file.close();
	}


	static void write_mtl(const std::string& filename, const std::string& texture_filename) {
		std::ofstream mtl_file(filename);

		if (!mtl_file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return;
		}


		mtl_file << "newmtl mat_texture" << std::endl;
		mtl_file << "Ka 1.000 1.000 1.000" << std::endl;  
		mtl_file << "Kd 1.000 1.000 1.000" << std::endl;  
		mtl_file << "Ks 0.000 0.000 0.000" << std::endl; 
		mtl_file << "d 1.0" << std::endl;                 
		mtl_file << "map_Kd " << texture_filename << std::endl; 

		mtl_file.close();
	}
};


class MetaData {
	std::string contents;

	bool timerOn = false;
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	double allTime = 0;

public:
	void append(std::string key, double val) {
		contents += "# " + key + ": " + std::to_string(val) + "\n";
	}
	void append(std::string key, int val) {
		contents += "# " + key + ": " + std::to_string(val) + "\n";
	}


	void startTimer() {
		if (!timerOn) {
			start_time = std::chrono::high_resolution_clock::now();
			timerOn = true;
		}
		else
			std::cout << "start failed" << std::endl;
	}
	void endTimer(std::string key) {
		if (timerOn) {
			end_time = std::chrono::high_resolution_clock::now();
			timerOn = false;
			//Milliseconds
			double time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
			allTime += time;
			append(key, time);
		}
		else
			std::cout << "end failed" << std::endl;

	}
	void clearTimer(void) {
		timerOn = false;
	}
	void toFile(std::string filename) {
		std::cout << "output verbose: " << filename << std::endl;
		append("time---all", allTime);
		std::ofstream ofs(filename + ".txt");
		ofs << contents;
		ofs.close();
	}


};


template <typename T1, typename T2>
class DynamicRowMatrix {
private:
	T1 _mat;
	int _rowcnt, _cols;

public:
	DynamicRowMatrix(int cols) : _cols(cols) {
		_rowcnt = 0;
	}
	DynamicRowMatrix(T1& initmat) {
		_mat = initmat;
		_rowcnt = _mat.rows();
		_cols = _mat.cols();
	}

	int rowcnt() { return _rowcnt; }

	int insert(T2& row) {
		if (_rowcnt == _mat.rows()) {
			_mat.conservativeResize(_mat.rows() * 2 + 1, _cols);
		}
		//_mat.row(_rowcnt) = row;
		for (int i = 0; i < row.size(); i++)
			_mat(_rowcnt, i) = row(i);
		return _rowcnt++;
		//++_rowcnt;
	}

	void append(T1& appmat) {
		assert(_cols == appmat.cols());
		_mat.conservativeResize(_rowcnt + appmat.rows(), _cols);
		_mat.block(_rowcnt, 0, appmat.rows(), appmat.cols()) = appmat;
		_rowcnt = _mat.rows();
	}

	T1 matXd() {
		_mat.conservativeResize(_rowcnt, _cols);
		return _mat;
	}
};



class KNNSearch {
	PointCloud<double> cloud;
	using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
		PointCloud<double>, 3 /* dim */
	>;

	my_kd_tree_t* tree = nullptr;

public:
	KNNSearch();
	KNNSearch(const Eigen::MatrixXd& points) {
		init(points);
	}
	~KNNSearch() {
		delete tree;
	}
	void init(const Eigen::MatrixXd& points) {
		cloud.pts.clear();
		for (int i = 0; i < points.rows(); ++i) {
			cloud.pts.push_back(
				PointCloud<double>::Point{
					points(i, 0),
					points(i, 1),
					points(i, 2)
				}
			);
		}
		if (tree) delete tree;
		auto start_nano = std::chrono::high_resolution_clock::now();
		tree = new my_kd_tree_t(3 /*dim*/, cloud, { 10 /* max leaf */ });
		auto end_nano = std::chrono::high_resolution_clock::now();
	}


	void nnSearch(
		const Eigen::VectorXd& pt,
		Eigen::MatrixXd& PTS,
		Eigen::VectorXd& DIS,
		int k = 5
	) {
		double query_pt[3] = { pt(0), pt(1), pt(2) };

		size_t                num_results = k;
		std::vector<uint32_t> ret_index(num_results);
		std::vector<double>    out_dist_sqr(num_results);

		num_results = tree->knnSearch(
			&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]
		);

		PTS = Eigen::MatrixXd::Zero(num_results, 3);
		DIS = Eigen::VectorXd::Zero(num_results);

		for (int i = 0; i < ret_index.size(); i++) {
			int id = ret_index[i];
			double dis = out_dist_sqr[i];
			PTS(i, 0) = cloud.pts[id].x;
			PTS(i, 1) = cloud.pts[id].y;
			PTS(i, 2) = cloud.pts[id].z;
			DIS(i) = dis;
		}
	}


	void rSearch(
		const Eigen::VectorXd& pt,
		Eigen::VectorXi& ID,
		//Eigen::MatrixXd& PTS,
		Eigen::VectorXd& DIS,
		double search_radius
	) {
		double query_pt[3] = { pt(0), pt(1), pt(2) };
		std::vector<nanoflann::ResultItem<uint32_t, double>> ret_matches;

		const size_t nMatches =
			tree->radiusSearch(&query_pt[0], search_radius, ret_matches);

		ID.conservativeResize(nMatches);
		//PTS = Eigen::MatrixXd::Zero(nMatches, 3);
		DIS.conservativeResize(nMatches);
		//std::cout << nMatches << "  " << ret_matches.size() << std::endl;
		for (int i = 0; i < ret_matches.size(); i++) {
			int id = ret_matches[i].first;
			double dis = ret_matches[i].second;
			/*PTS(i, 0) = cloud.pts[id].x;
			PTS(i, 1) = cloud.pts[id].y;
			PTS(i, 2) = cloud.pts[id].z;*/
			ID(i) = id;
			DIS(i) = dis;
		}
	}


	std::pair<std::vector<int>, std::vector<double>> rSearch(
		const Eigen::VectorXd& pt,
		const double search_radius
	) {
		double query_pt[3] = { pt(0), pt(1), pt(2) };
		std::vector<nanoflann::ResultItem<uint32_t, double>> ret_matches;

		const size_t nMatches =
			tree->radiusSearch(&query_pt[0], search_radius, ret_matches);
		std::vector<int> ID;
		std::vector<double> DIS;

		ID.reserve(ret_matches.size());
		DIS.reserve(ret_matches.size());
#pragma omp parallel for		
		for (int i = 0; i < ret_matches.size(); i++) {
			int id = ret_matches[i].first;
			double dis = ret_matches[i].second;

			ID[i] = (id);
			DIS[i] = (dis);
		}

		return std::make_pair(ID, DIS);
	}


};

