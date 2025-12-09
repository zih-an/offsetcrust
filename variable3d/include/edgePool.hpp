#pragma once

#include <iostream>
#include <Eigen/Core>
#include <map>
#include <set>
#include <vector>


namespace Offset3D {
	struct EdgeInfo {
		int v0, v1;
		Eigen::VectorXd vv0, vv1;
		double len, corner_ratio;
		std::vector<int> adj_faces;

		int isConvex = -1;  //-1: planes; 0: concave; 1: convex
		void setEdgeInfo(
			int _v0, int _v1,
			Eigen::VectorXd& _vv0, Eigen::VectorXd& _vv1
		) {
			v0 = _v0;
			v1 = _v1;
			vv0 = _vv0;
			vv1 = _vv1;

			len = (_vv1 - _vv0).norm();
			double ideal_point_step = 0.001;  //need normalization scale
			//double ideal_point_step = 0.0001;
			corner_ratio = std::min(0.05, ideal_point_step / len);
			//corner_ratio = ratio;
		}
		bool isSetVal = false;
		bool operator==(const EdgeInfo& other) const {
			if ((this->v0 == other.v0 && this->v1 == other.v1) || (this->v1 == other.v0 && this->v0 == other.v1))
				return true;
			return false;
		}

		bool getEdgeRatio(int from, int to, Eigen::VectorXd& p) {
			bool inverse;
			if (from == v0 && to == v1) inverse = false;
			else if (to == v0 && from == v1) inverse = true;
			else return false;

			double ratio = inverse ? 1 - corner_ratio : corner_ratio;
			p = ratio * (vv1 - vv0) + vv0;
			return true;
		}

		Eigen::VectorXd getEdgeVec(int start) {
			if (start == v0) return vv1 - vv0;
			else if (start == v1) return vv0 - vv1;
			else std::cout << "wrong" << std::endl;
		}

	};

	class EdgePool {
	private:
		std::map<std::set<int>, EdgeInfo*> edges;
		double __ratio;

	public:
		EdgePool(double ratio) {
			__ratio = ratio;
		}

		~EdgePool() {
			for (auto& pair : edges) {
				delete pair.second;
			}
		}

		EdgeInfo* getEdge(int v0, int v1) {
			std::set<int> edge_key = { v0, v1 };
			auto it = edges.find(edge_key);
			if (it != edges.end()) {
				return it->second;
			}
			EdgeInfo* new_edge = new EdgeInfo();
			new_edge->v0 = v0;
			new_edge->v1 = v1;
			new_edge->corner_ratio = __ratio;
			edges[edge_key] = new_edge;
			return new_edge;
		}

		void forEachEdgeInfo(const std::function<void(EdgeInfo*)>& callback) const {
			for (const auto& [key, edgeInfo] : edges) {
				callback(edgeInfo);
			}
		}

	};
};
