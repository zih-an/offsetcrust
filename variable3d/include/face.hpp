#pragma once

#include <cassert>
#include <set>
#include <vector>
#include <iostream>
#include <cgal_types_3d.hpp>


namespace Offset3D {
	//now vorface struct
	struct VorFaceTri {
	public:
		std::vector<Point> points;
		int tri_pid, tri_qid;

		VorFaceTri() { }
		VorFaceTri(
			std::vector<Point> points,
			int tri_pid, int tri_qid
			//Eigen::VectorXd face_normal
		)
			:points(points), tri_pid(tri_pid), tri_qid(tri_qid)
		{ }

	};
}
