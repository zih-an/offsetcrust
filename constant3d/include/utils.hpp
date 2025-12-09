#pragma once

#include <Eigen/Core>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>


class UTILS {
public:
	static std::string data_dir;
	static std::string input_dir, out_dir, sphere_file;

	static void normalization(Eigen::MatrixXd& M) {
		M = (M - M.minCoeff() * Eigen::MatrixXd::Ones(M.rows(), M.cols())) / (M.maxCoeff() - M.minCoeff());
	}
	static void normalization_2(Eigen::MatrixXd& M) {
		M = (M - M.minCoeff() * Eigen::MatrixXd::Ones(M.rows(), M.cols())) / (M.maxCoeff() - M.minCoeff());
		M = M * 2 - Eigen::MatrixXd::Ones(M.rows(), M.cols());  // (-1,1)
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
		std::cout << "output xyz..." << UTILS::out_dir + filename + "_normal.xyz" << std::endl;
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
	}

	static void outputXYZ(
		Eigen::MatrixXd& Vertices,
		std::string filename
	) {
		std::cout << "output xyz... verts" << std::endl;
		std::ofstream file(UTILS::data_dir + filename + "_verts.xyz");
		for (int i = 0; i < Vertices.rows(); i++) {
			file <<
				std::to_string(Vertices(i, 0)) + " "
				+ std::to_string(Vertices(i, 1)) + " "
				+ std::to_string(Vertices(i, 2)) + "\n";
		}
		file.close();
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
		for (int i = 0; i < row.size(); i++)
			_mat(_rowcnt, i) = row(i);
		return _rowcnt++;
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
	



