#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <cassert>

class OutputOBJ {

private:
	std::string res, faceMtl, vn_res;
	bool useFaceMtl = false;
	int colorid = 0;


public:
	void addVertices(Eigen::MatrixXd& V) {
		for (int i = 0; i < V.rows(); i++) {
			res += "v " + std::to_string(V(i, 0)) + " " + std::to_string(V(i, 1)) + " " + std::to_string(V(i, 2)) + "\n";
		}
	}
	void addFaces(Eigen::MatrixXi& F) {
		for (int i = 0; i < F.rows(); i++) {
			res += "f "
				+ std::to_string(F(i, 0) + 1)
				+ " " + std::to_string(F(i, 1) + 1)
				+ " " + std::to_string(F(i, 2) + 1)
				+ "\n";
		}
	}

	void addVertices(Eigen::MatrixXd& V, Eigen::MatrixXd& nV) {
		assert(V.rows() == nV.rows());
		for (int i = 0; i < V.rows(); i++) {
			res += "v " + std::to_string(V(i, 0)) + " " + std::to_string(V(i, 1)) + " " + std::to_string(V(i, 2)) + "\n";
		}

		for (int i = 0; i < nV.rows(); i++) {
			res += "vn " + std::to_string(nV(i, 0)) + " " + std::to_string(nV(i, 1)) + " " + std::to_string(nV(i, 2)) + "\n";
		}
	}

	void addFaces(std::vector<int>& face) {
		res += "f ";
		for (auto vid : face)
			res += std::to_string(vid + 1) + " ";
		if (!res.empty())
			res.pop_back();
		res += "\n";
	}


	void addVN(Eigen::VectorXd& face_normal) {
		res += "vn "
			+ std::to_string(face_normal(0))
			+ " " + std::to_string(face_normal(1))
			+ " " + std::to_string(face_normal(2)) + "\n";
	}

	void addFaces(std::vector<int>& face, int face_id) {
		/*res += "vn "
			+ std::to_string(face_normal(0))
			+ " " + std::to_string(face_normal(1))
			+ " " + std::to_string(face_normal(2)) + "\n";*/

		res += "f ";
		for (auto vid : face)
			res += std::to_string(vid + 1) + "//" + std::to_string(face_id + 1) + " ";
		if (!res.empty())
			res.pop_back();
		res += "\n";
	}

	void addFaces(std::set<int>& face) {
		res += "f ";
		for (auto vid : face)
			res += std::to_string(vid + 1) + " ";
		if (!res.empty())
			res.pop_back();
		res += "\n";
	}

	void addFaces(std::vector<int>& face, double r, double g, double b) {
		useFaceMtl = true;
		std::string mtlname = "facecolor" + std::to_string(colorid++) + "\n";
		res += "usemtl " + mtlname;
		addFaces(face);

		faceMtl += "newmtl " + mtlname;
		faceMtl += "Kd " + std::to_string(r) + " " + std::to_string(g) + " " + std::to_string(b) + "\n";
	}

	void addLine(std::set<int>& edge) {
		res += "l ";
		for (auto vid : edge)
			res += std::to_string(vid + 1) + " ";
		if (!res.empty())
			res.pop_back();
		res += "\n";
	}

	void addLine(int v1, int v2) {
		res += "l ";
		res += std::to_string(v1 + 1) + " " + std::to_string(v2 + 1);
		res += "\n";
	}

	//filename without .obj
	void toFile(std::string filename) {
		std::string mtlFilename = filename + ".mtl";
		filename += ".obj";
		std::ofstream fs(filename);

		if (useFaceMtl) {
			std::ofstream mtlfs(mtlFilename);
			std::cout << mtlFilename << std::endl;
			if (mtlfs.is_open()) {
				mtlfs << faceMtl;
				mtlfs.close();
			}
		}

		std::cout << filename << std::endl;
		if (fs.is_open()) {
			if (useFaceMtl) fs << "mtllib " + mtlFilename + "\n";
			fs << res;
			fs.close();
		}

	}

};

