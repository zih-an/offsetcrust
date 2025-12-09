#include <offset_crust.hpp>
using namespace Offset3D;



OffsetCrust::OffsetCrust(
	std::string filename,
	int _num_angle,
	int blue_points,
	double corner_ratio,
	double _lambda
) {
	seg_angle = _num_angle;
	blue_noise_desired = blue_points;
	ratio = corner_ratio;
	lambda = _lambda;

	//0. read mesh (off file)
	igl::read_triangle_mesh(filename, OrgVertices, OrgFaces);

	if (VERBOSE_ON) {
		VERBOSE.append("original num V", int(OrgVertices.rows()));
		VERBOSE.append("original num F", int(OrgFaces.rows()));
	}
	// bouding box
	Eigen::AlignedBox<double, 3> bounding_box;
	bounding_box.setNull();
	for (int i = 0; i < OrgVertices.rows(); ++i) {
		Eigen::Vector3d p= OrgVertices.row(i);
		bounding_box.extend(p);
	}
	diagonal_length = bounding_box.diagonal().norm();

	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
	for (int i = 0; i < OrgFaces.rows(); i++) {
		Eigen::VectorXd 
			v1 = OrgVertices.row(OrgFaces(i, 0)),
			v2 = OrgVertices.row(OrgFaces(i, 1)),
			v3 = OrgVertices.row(OrgFaces(i, 2));

		simple_K::Point_3
			a(v1(0), v1(1), v1(2)),
			b(v2(0), v2(1), v2(2)),
			c(v3(0), v3(1), v3(2));
		triangles.push_back(Triangle(a,b,c));
	}
	std::cout << "AABB tree for orginal mesh..." << std::endl;
	tree = Tree(triangles.begin(), triangles.end());
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---AABB construction");
	}



	std::cout << "sphere: " << UTILS::sphere_file << std::endl;
	igl::readOBJ(UTILS::sphere_file, SPHV, SPHF);
	igl::per_vertex_normals(SPHV, SPHF, SPHVN);
	UTILS::normalization_2(SPHV);
}


void OffsetCrust::generateWeights(
	Eigen::MatrixXd& AllVertices,
	Eigen::VectorXd& AllWeights,
	double w_diff
) {
	for (int i = 0; i < AllVertices.rows(); i++) {
		if (i < PSize)
			AllWeights(i) = 0;
		else
			AllWeights(i) = w_diff;
	}
}


void OffsetCrust::compute_power_diagram(
	Eigen::MatrixXd& AllVertices,
	Eigen::MatrixXd& QNormals,
	Eigen::VectorXi& blue_validFlags,
	double w_diff,
	std::string vorFilename
) {
	//1. generate weights
	Eigen::VectorXd AllWeights(AllVertices.rows());
	generateWeights(
		AllVertices, AllWeights, w_diff
	);

	Eigen::VectorXi outFlag(PreQVerts.rows());
	std::cout << "check Q..." << std::endl;
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
#pragma omp parallel for
	for (int i = 0; i < PreQVerts.rows(); i++) {
		Eigen::VectorXd v = PreQVerts.row(i) + 100*q_offset * QNormals.row(i);
		simple_K::Point_3 a(v(0), v(1), v(2));
		// counts #intersections
		Eigen::Vector3d ref = Eigen::Vector3d::Random();
		simple_K::Point_3 REF(ref(0), ref(1), ref(2));
		simple_K::Ray_3 ray_query(a, REF);
		int num = tree.number_of_intersected_primitives(ray_query);
		int isOut = 0;
		if (num % 2 == 0) isOut = 1;
		outFlag(i) = isOut;
	}
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---AABB checkQ");
	}

	//2. compute
	Regular_triangulation rt;
	std::vector<std::pair<Weighted_point, unsigned>> w_points;
	for (int i = 0; i < PSize; i++) {
		std::vector<int> validQ;
		for (auto real_qid : pi2qi[i]) {
			if (!outFlag(real_qid) && __IS_INSIDE || outFlag(real_qid) && !__IS_INSIDE) {
				bool flag = true;
				if (real_qid < blue_validFlags.rows()) {
					flag = flag && blue_validFlags(i);
				}
				if(flag) validQ.push_back(real_qid);
			}
		}
		if (!validQ.empty()) {
			w_points.push_back(std::make_pair(
				Weighted_point(
					Point(AllVertices(i, 0), AllVertices(i, 1), AllVertices(i, 2)),
					AllWeights(i)
				),
				i
			));
			for (int real_qid : validQ) {
				int qid = real_qid + PSize;
				w_points.push_back(std::make_pair(
					Weighted_point(
						Point(AllVertices(qid, 0), AllVertices(qid, 1), AllVertices(qid, 2)),
						AllWeights(qid)
					),
					qid
				));
			}
		}
	}
	if (VERBOSE_ON) {
		VERBOSE.append("sample ALL", int(w_points.size()));
	}
	AllWeights.resize(0);
	outFlag.resize(0);

	std::cout << "power diagram computing... " << w_points.size() << "... ";
#ifdef CGAL_LINKED_WITH_TBB
	std::cout << "parallel" << std::endl;
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
	Regular_triangulation::Lock_data_structure locking_ds(bbox, 50);
	rt = Regular_triangulation(w_points.begin(), w_points.end(), &locking_ds);
	assert(rt.is_valid());
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---powerdiagram");
	}
#else
	std::cout << "sequential" << std::endl;
	rt.insert(w_points.begin(), w_points.end());
	assert(rt.is_valid());
#endif // CGAL_LINKED_WITH_TBB




	//3. extract voronoi.
	std::vector<VorFaceTri*> vor_diagram;
	extract_voronoi(rt, QNormals, vor_diagram);
	optim_vor_diagram(vor_diagram, QNormals, AllVertices, vorFilename);
}




//voronoi faces to triangle mesh.
void OffsetCrust::extract_voronoi(
	Regular_triangulation& rt,
	Eigen::MatrixXd& QNormals,
	std::vector<VorFaceTri*>& vor_diagram  //output (voronoi)
) {
	std::cout << "extracting voronoi information..." << std::endl;
	std::vector<Edge_iterator> edges_iter;
	edges_iter.reserve(rt.number_of_finite_edges());
	vor_diagram.reserve(rt.number_of_finite_edges());  //actual memory
	for (Edge_iterator eit = rt.finite_edges_begin(); eit != rt.finite_edges_end(); ++eit)
		edges_iter.push_back(eit);

	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
#pragma omp parallel for
	for (int i = 0; i < rt.number_of_finite_edges(); i++) {
		Edge_iterator eit = edges_iter[i];

		/* 0. preserve manifold ? */
		Vertex_handle source = eit->first->vertex(eit->second),
			target = eit->first->vertex(eit->third);
		int pid = std::min(source->info(), target->info()),
			qid = std::max(source->info(), target->info());
		if (pid < PSize && qid < PSize || pid >= PSize && qid >= PSize)
			continue;

		int real_qid = qid - PSize;

		/*1. get tet around the edge*/
		Cell_circulator
			cell_circulator = rt.incident_cells(*eit),
			done(cell_circulator);
		if (cell_circulator == nullptr) continue;

		std::vector<Point> vor_face;
		bool valid = true;
		do {
			if (rt.is_infinite(cell_circulator)) {
				valid = false;
				break;
			}

			Point p;
			p = rt.dual(cell_circulator);
			Eigen::VectorXd pvec(3);
			pvec <<
				CGAL::to_double(p.x()),
				CGAL::to_double(p.y()),
				CGAL::to_double(p.z());

			if (OFFSET > 1e-6) {
				simple_K::Point_3 a(pvec(0), pvec(1), pvec(2));
				 //counts #intersections
				Eigen::Vector3d ref = Eigen::Vector3d::Random();
				simple_K::Point_3 REF(ref(0), ref(1), ref(2));
				simple_K::Ray_3 ray_query(a, REF);
				int num = tree.number_of_intersected_primitives(ray_query);
				int isOut = 0;
				if (num % 2 == 0) isOut = 1;
				if (isOut && __IS_INSIDE || !isOut && !__IS_INSIDE) {
					valid = false;
					break;
				}
			}

			
			if (std::find(vor_face.begin(), vor_face.end(), p) == vor_face.end())
			{
				vor_face.push_back(p);
			}
		} while (++cell_circulator != done);
		if (!valid || vor_face.size() <= 2)
			continue;
		

#pragma omp critical
		{
			VorFaceTri* face = new VorFaceTri(vor_face, pid, qid);
			vor_diagram.push_back(face); 
		}
	}
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---extract face");
	}

}



void OffsetCrust::optim_vor_diagram(
	std::vector<VorFaceTri*>& vor_diagram,
	Eigen::MatrixXd& QNormals,
	Eigen::MatrixXd& AllVerts,  //QVerts.rows = QNormals.rows
	std::string vorFilename
) {
	std::cout << "collecting adjacent information... " << std::endl;
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
	std::map<Point, std::set<VorFaceTri*>> vor_vertices;
	for (int i = 0; i < vor_diagram.size(); i++) {
		for (Point p : vor_diagram[i]->points) {
			vor_vertices[p].insert(vor_diagram[i]);
		}
	}
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---pts adjinfo");
	}
	std::cout << "vor_vertices: " << vor_vertices.size() << std::endl;


	Eigen::MatrixXd
		vorV(vor_vertices.size(), 3),
		vorNewV(vor_vertices.size(), 3);
	std::map<Point, int> point_id_inVerts;

	//iterators for parallel
	std::vector< std::map<Point, std::set<VorFaceTri*>>::iterator> iters;
	for (auto it = vor_vertices.begin(); it != vor_vertices.end(); ++it)
		iters.push_back(it);

	std::cout << "QEM optimization..." << std::endl;
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}
#pragma omp parallel for
	for (int ii = 0; ii < iters.size(); ii++) {
		auto it = iters[ii];
		Point p = it->first;

		Eigen::VectorXd pV(3);
		pV <<
			CGAL::to_double(p.x()),
			CGAL::to_double(p.y()),
			CGAL::to_double(p.z());
		//1;

	//QEM
		Eigen::MatrixXd AllK = Eigen::MatrixXd::Zero(4, 4);
		for (auto face : it->second) {
			int pid = face->tri_pid, qid = face->tri_qid;
			int real_qid = qid - PSize;
			//Qid
			Eigen::VectorXd
				qnormal = QNormals.row(real_qid),
				q_onface = AllVerts.row(qid); 
			q_onface += (OFFSET - q_offset) * qnormal;

			//Pid
			if (pi2qi[pid].size() >= 2) {  //1 vs n sample. q_noraml only
				Eigen::MatrixXd qK(4, 4);
				QEM_K(qnormal, -qnormal.transpose() * q_onface, qK);
				AllK += qK;
				continue;
			}

			int pid_q = pi2qi[pid][0];
			Eigen::VectorXd
				pnormal = QNormals.row(pid_q),
				p_onface = AllVerts.row(pid);
			p_onface += OFFSET * pnormal;


			if (qnormal.dot(pnormal) + 1 < 1e-6) {
				//nearly 2pi(180��) must be plane.
				continue;
			}
			Eigen::MatrixXd qK(4, 4);
			QEM_K(qnormal, -qnormal.transpose() * q_onface, qK);
			AllK += qK;
			Eigen::MatrixXd pK(4, 4);
			QEM_K(pnormal, -pnormal.transpose() * p_onface, pK);
			AllK += pK;
		}

		Eigen::VectorXd base = lambda * pV - AllK.block(0, 3, 3, 1);
		Eigen::MatrixXd Iden = AllK.block(0, 0, 3, 3) + lambda * Eigen::MatrixXd::Identity(3, 3);
		Eigen::VectorXd tmp_newV = Iden.ldlt().solve(base);

#pragma omp critical
		{
			point_id_inVerts[p] = ii;
			vorV.row(ii) <<
				pV(0), pV(1), pV(2);

			vorNewV.row(ii) <<
				tmp_newV(0), tmp_newV(1), tmp_newV(2);

		}
	}
	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---QEM optim");
	}


	/*---------------------------output results---------------------------*/
	std::cout << "output poly obj..." << std::endl;
	OutputOBJ obj, objNew;
	obj.addVertices(vorV);
	objNew.addVertices(vorNewV);
	std::vector<std::vector<int>> polyF;
	for (int i = 0; i < vor_diagram.size(); i++) {
		std::vector<int> face;
		for (Point p : vor_diagram[i]->points) {
			int vid = point_id_inVerts[p];
			face.push_back(vid);
		}
		obj.addFaces(face);
		objNew.addFaces(face);
		polyF.push_back(face);
	}

	if (__outputPoly) {
		obj.toFile(vorFilename + "_poly");
		objNew.toFile(vorFilename + "_poly_optim");
	}

	std::cout << "postprocessing..." << std::endl;
	Eigen::MatrixXd cleanVorV, cleanVorNewV;
	Eigen::MatrixXi cleanVorF, cleanVorNewF;
	postproc_remove_dup(vorV, polyF, cleanVorV, cleanVorF);
	postproc_remove_dup(vorNewV, polyF, cleanVorNewV, cleanVorNewF);

	std::cout << "output tri obj (unduplicated)..." << std::endl;
	igl::writeOBJ(vorFilename + "_tri.obj", cleanVorV, cleanVorF);
	igl::writeOBJ(vorFilename + "_tri_optim.obj", cleanVorNewV, cleanVorNewF);


	/*delete all vorfaces*/
	for (int i = 0; i < vor_diagram.size(); i++) {
		delete vor_diagram[i];
	}
}



void OffsetCrust::postproc_remove_dup(
	Eigen::MatrixXd& polyV,
	std::vector<std::vector<int>>& polyF,
	Eigen::MatrixXd& newV,
	Eigen::MatrixXi& newF
) {
	std::cout << "remove duplicate & triangulate..." << std::endl;
	DynamicRowMatrix<Eigen::MatrixXi, Eigen::VectorXi> DyF(3);
	for (int i = 0; i < polyF.size(); i++) {
		for (int j = 2; j < polyF[i].size(); ++j) {
			Eigen::VectorXi fv(3);
			fv << polyF[i][0], polyF[i][j - 1], polyF[i][j];
			DyF.insert(fv);
		}
	}
	Eigen::MatrixXi F = DyF.matXd();

	Eigen::MatrixXi F1, F2;
	Eigen::VectorXi SVI, SVJ, J;
	igl::remove_duplicate_vertices(polyV, F, 1e-6, newV, SVI, SVJ, F1);
	igl::resolve_duplicated_faces(F1, F2, J);

	DynamicRowMatrix<Eigen::MatrixXi, Eigen::VectorXi> DyFnew(3);
	for (int i = 0; i < F2.rows(); i++) {
		if (F2(i, 0) == F2(i, 1) || F2(i, 0) == F2(i, 2) || F2(i, 1) == F2(i, 2))
			continue;
		Eigen::VectorXi tmp = F2.row(i);
		DyFnew.insert(tmp);
	}
	newF = DyFnew.matXd();
}



void OffsetCrust::QEM_K(
	Eigen::VectorXd& normal,
	double d,
	Eigen::MatrixXd& K
) {
	K.conservativeResize(4, 4);
	K << normal(0) * normal(0), normal(0)* normal(1), normal(0)* normal(2), normal(0)* d,
		normal(1)* normal(0), normal(1)* normal(1), normal(1)* normal(2), normal(1)* d,
		normal(2)* normal(0), normal(2)* normal(1), normal(2)* normal(2), normal(2)* d,
		d* normal(0), d* normal(1), d* normal(2), d* d;
}



void OffsetCrust::preprocessing_with_sphere(
	Eigen::MatrixXd& Vertices,
	Eigen::MatrixXi& Faces,
	Eigen::MatrixXd& OutAllVertices,
	Eigen::MatrixXd& OutQNormals,
	Eigen::VectorXi& blue_validFlags,
	double q_orient
) {
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}

	//1. presets
	igl::per_face_normals(Vertices, Faces, OrgFaceNormals);
	PQP = new ProjPoint(OrgVertices, OrgFaces, OrgFaceNormals, __IS_INSIDE);
	if (__IS_INSIDE) OrgFaceNormals = -OrgFaceNormals;


	// Heuristic to  determine radius from desired number 
	Eigen::MatrixXd P_blue, N_blue, Centroids;
	Eigen::MatrixXd outV, outVN;
	Eigen::VectorXi FI;
	std::map<int, std::vector<int>> f2blueP;

	if (__useCentroid) {
		Centroids.conservativeResize(Faces.rows(), 3);
		for (int fid = 0; fid < Faces.rows(); fid++) {
			int v0 = Faces(fid, 0),
				v1 = Faces(fid, 1),
				v2 = Faces(fid, 2);
			Eigen::VectorXd
				v0_p = Vertices.row(v0),
				v1_p = Vertices.row(v1),
				v2_p = Vertices.row(v2);
			Eigen::VectorXd cen = (v0_p + v1_p + v2_p) / 3;
			Centroids.row(fid) << cen(0), cen(1), cen(2);

			qi2pi[fid] = fid;
			pi2qi[fid].push_back(fid);
			f2blueP[fid].push_back(fid);
		}
		outV = Centroids, outVN = OrgFaceNormals;
		std::cout << "[Centroids] " << Centroids.rows() << std::endl;
	}
	else {
		// Heuristic to  determine radius from desired number 
		const double r = [&](const int n) {
			Eigen::VectorXd A;
			igl::doublearea(OrgVertices, OrgFaces, A);
			return sqrt(((A.sum() * 0.5 / (n * 0.6162910373)) / igl::PI));
			}(blue_noise_desired);

			printf("blue noise radius: %g\n", r);
			Eigen::MatrixXd B;
			igl::blue_noise(OrgVertices, OrgFaces, r, B, FI, P_blue);
			N_blue.conservativeResize(P_blue.rows(), P_blue.cols());
			for (int i = 0; i < FI.size(); i++) {
				int id = FI(i);
				N_blue.row(i) = OrgFaceNormals.row(id);
				qi2pi[i] = i;
				pi2qi[i].push_back(i);
				f2blueP[id].push_back(i);
			}
			outV = P_blue, outVN = N_blue;
			std::cout << "[Blue Noise] " << P_blue.rows() << std::endl;
	}
	
	Eigen::MatrixXd outVQ = outV + q_orient * q_offset * outVN;
	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd>
		DyPVertices(outV),
		DyQVertices(outVQ),
		DyQNormals(outVN);

	blue_validFlags = Eigen::VectorXi::Ones(outV.rows());

	//2. get all edges
	EdgePool edges(ratio);
	std::map<int, std::set<EdgeInfo*>> v2edges;
	std::map<int, std::vector<int>> v2faces;
	for (int fid = 0; fid < Faces.rows(); fid++) {
		for (int i = 0; i < 3; i++) {
			int v0 = Faces(fid, i),
				v1 = Faces(fid, (i + 1) % 3);
			EdgeInfo* e = edges.getEdge(v0, v1);
			if (!e->isSetVal) {
				Eigen::VectorXd
					vv0 = Vertices.row(v0), vv1 = Vertices.row(v1);
				e->setEdgeInfo(v0, v1, vv0, vv1);
			}
			e->adj_faces.push_back(fid);
			v2faces[v0].push_back(fid);
			v2edges[v0].insert(e);
		}
	}
	std::cout << "[ratio=] " << ratio << std::endl;
	std::cout << "[q_offset=] " << q_offset << std::endl;

	
	double slerpDegree = ((90. / seg_angle) * M_PI / 180.0);
	std::set<int> v_corners;
	edges.forEachEdgeInfo([&](EdgeInfo* edge) {
		if (edge->adj_faces.size() != 2) return;

		int f1 = edge->adj_faces[0], f2 = edge->adj_faces[1];
		Eigen::VectorXd
			n1 = OrgFaceNormals.row(f1), n2 = OrgFaceNormals.row(f2);
		double cosine = n1.transpose() * n2;
		cosine = std::clamp(cosine, -1.0, 1.0);
		if (cosine >= diheralBar) //not feature
			return;

		int v1 = edge->v0, v2 = edge->v1,
			f1v0, f2v0;
		Eigen::VectorXd
			p1 = edge->vv0, p2 = edge->vv1,
			uniq_f2v(3), uniq_f1v(3);

		//classify: convex / concave
		for (int i = 0; i < 3; i++) {
			int f1v = Faces(f1, i);
			if (f1v != v1 && f1v != v2) {
				f1v0 = f1v;
				uniq_f1v = Vertices.row(f1v);
			}

			int f2v = Faces(f2, i);
			if (f2v != v1 && f2v != v2) {
				f2v0 = f2v;
				uniq_f2v = Vertices.row(f2v);
			}
		}
		Eigen::VectorXd jud_vec = uniq_f2v - p1;
		double jud_val = jud_vec.transpose() * n1;
		edge->isConvex = (jud_val >= 0) ? 0 : 1;

		int num_points = (edge->corner_ratio > q_offset) ? 1. / (edge->corner_ratio) : 0;
		if (ratio > q_offset) {
			EdgeInfo
				* f1e1 = edges.getEdge(f1v0, v1), * f1e2 = edges.getEdge(f1v0, v2),
				* f2e1 = edges.getEdge(f2v0, v1), * f2e2 = edges.getEdge(f2v0, v2);
			Eigen::VectorXd
				f1sp1, f1sp2, f2sp1, f2sp2;
			f1e1->getEdgeRatio(v1, f1v0, f1sp1);
			f1e2->getEdgeRatio(v2, f1v0, f1sp2);
			f2e1->getEdgeRatio(v1, f2v0, f2sp1);
			f2e2->getEdgeRatio(v2, f2v0, f2sp2);

			for (auto bluep : f2blueP[f1]) {
				Eigen::VectorXd bluevp = outV.row(bluep);
				bool flag = UTILS::isPointInTriangle(bluevp, f1sp1, f1sp2, uniq_f1v);
				blue_validFlags(bluep) = blue_validFlags(bluep) && flag ? 1 : 0;
			}
			for (auto bluep : f2blueP[f2]) {
				Eigen::VectorXd bluevp = outV.row(bluep);
				bool flag = UTILS::isPointInTriangle(bluevp, f2sp1, f2sp2, uniq_f2v);
				blue_validFlags(bluep) = blue_validFlags(bluep) && flag ? 1 : 0;
			}


			K::FT point_step = 1. / num_points;
			//if (jud_val >= 0) { 
				/*1v1*/
			for (K::FT beta = point_step; ratio > q_offset && beta < 1 - q_offset; beta += point_step) {
				Eigen::VectorXd
					f1sp = CGAL::to_double(beta) * f1sp1 + CGAL::to_double(1 - beta) * f1sp2,
					f2sp = CGAL::to_double(beta) * f2sp1 + CGAL::to_double(1 - beta) * f2sp2;
				int spid1 = DyPVertices.insert(f1sp),
					spid2 = DyPVertices.insert(f2sp);
				DyQNormals.insert(n1);
				DyQNormals.insert(n2);
				Eigen::VectorXd
					f1sq = f1sp + q_offset * n1,
					f2sq = f2sp + q_offset * n2;
				int sqid1 = DyQVertices.insert(f1sq),
					sqid2 = DyQVertices.insert(f2sq);
				qi2pi[sqid1] = spid1;
				pi2qi[spid1].push_back(sqid1);  //1 v n
				qi2pi[sqid2] = spid2;
				pi2qi[spid2].push_back(sqid2);  //1 v n
			}

			//}
		}


		v_corners.insert(v1);
		v_corners.insert(v2);
		if (edge->isConvex != 1) return;


		K::FT point_step = 1. / (num_points + 2);
		for (K::FT beta = point_step; beta < 1 - q_offset; beta += point_step) {
			Eigen::VectorXd beta_p = CGAL::to_double(beta) * p1 + CGAL::to_double(1 - beta) * p2;
			int pi_rowid;
			bool isPushed = false;

			double Omega = std::acos(cosine);
			int num_angle = std::floor(Omega / slerpDegree);
			double
				alpha_step = num_angle <= 0 ? 2 : 1. / num_angle,
				sOmega = std::sin(Omega);

			for (double alpha = 0; cosine <= 0.99 && alpha <= 1; alpha += alpha_step) {
				Eigen::VectorXd alpha_n
					= std::sin((1 - alpha) * Omega) / sOmega * n1
					+ std::sin(alpha * Omega) / sOmega * n2;

				if (!isPushed) {
					pi_rowid = DyPVertices.insert(beta_p);
					isPushed = true;
				}
				alpha_n.normalize();
				DyQNormals.insert(alpha_n);

				Eigen::VectorXd beta_q = beta_p + q_offset * alpha_n;
				int qi_rowid = DyQVertices.insert(beta_q);
				qi2pi[qi_rowid] = pi_rowid;
				pi2qi[pi_rowid].push_back(qi_rowid);  //1 v n
			}
		}

		});

	K::FT point_step = 1. / num_corner_points;
	for (auto corner : v_corners) {
		Eigen::VectorXd beta_p = Vertices.row(corner);
		std::set<EdgeInfo*> adj_edges = v2edges[corner];
		std::vector<int> adjF = v2faces[corner];
		for (int fi : adjF) {
			if (ratio <= q_offset) break;
			int ep1, ep2, vid = corner;
			for (int i = 0; i < 3; i++) {
				int v = Faces(fi, i);
				if (v == vid) {
					ep1 = Faces(fi, (i + 1) % 3);
					ep2 = Faces(fi, (i + 2) % 3);
					break;
				}
			}
			Eigen::VectorXd
				uniq_fv = Vertices.row(vid),
				n = OrgFaceNormals.row(fi);

			//exclude some blue noise sample points.
			EdgeInfo
				* f1e1 = edges.getEdge(vid, ep1), * f1e2 = edges.getEdge(vid, ep2);
			Eigen::VectorXd
				f1sp1, f1sp2;
			f1e1->getEdgeRatio(vid, ep1, f1sp1);
			f1e2->getEdgeRatio(vid, ep2, f1sp2);

			for (auto bluep : f2blueP[fi]) {
				Eigen::VectorXd bluevp = outV.row(bluep);
				bool flag = UTILS::isPointInTriangle(bluevp, f1sp1, f1sp2, uniq_fv);
				blue_validFlags(bluep) = blue_validFlags(bluep) && !flag ? 1 : 0;
			}


			for (K::FT beta = point_step; beta < 1 - q_offset; beta += point_step) {
				Eigen::VectorXd
					f1sp = CGAL::to_double(beta) * f1sp1 + CGAL::to_double(1 - beta) * f1sp2;
				int spid1 = DyPVertices.insert(f1sp);
				DyQNormals.insert(n);
				Eigen::VectorXd
					f1sq = f1sp + q_offset * n;
				int sqid1 = DyQVertices.insert(f1sq);
				qi2pi[sqid1] = spid1;
				pi2qi[spid1].push_back(sqid1);  //1 v n
			}
		}

		std::vector<EdgeInfo*> convex_edges;
		for (auto edge : adj_edges) {
			if (edge->isConvex == 1) 
				convex_edges.push_back(edge);
		}
		bool isPushed = false;
		int pi_rowid;

		if (convex_edges.size() < 2)
			continue;
		
		/* 1vn */
		Eigen::VectorXd ones = Eigen::VectorXd::Ones(SPHV.rows());
		Eigen::MatrixXd trans_SPHV = (SPHV * q_offset + ones * beta_p.transpose());

		for (int svid = 0; svid < trans_SPHV.rows(); svid++) {
			Eigen::VectorXd
				sv = trans_SPHV.row(svid),
				alpha_n = SPHVN.row(svid);

			if (!isPushed) {
				pi_rowid = DyPVertices.insert(beta_p);
				isPushed = true;
			}
			DyQNormals.insert(alpha_n);
			int qi_rowid = DyQVertices.insert(sv);

			qi2pi[qi_rowid] = pi_rowid;
			pi2qi[pi_rowid].push_back(qi_rowid);  //1 v n
		}
		
	}


	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---preproc");
	}

	//4. final results
	PrePVerts = DyPVertices.matXd();
	PreQVerts = DyQVertices.matXd();
	PSize = DyPVertices.rowcnt();
	DyPVertices.append(PreQVerts);
	OutAllVertices = DyPVertices.matXd();
	OutQNormals = DyQNormals.matXd();


	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd>
		DyTestOut(3), DyTestNOut(3),
		DySharpOut(3), DySharpNOut(3);
	for (int i = 0; i < PreQVerts.rows(); i++) {
		int pid = qi2pi[i];
		Eigen::VectorXd
			qt = PreQVerts.row(i),
			qn = OutQNormals.row(i);
		if (pid < blue_validFlags.rows() && blue_validFlags(pid)==1) {
			DyTestOut.insert(qt);
			DyTestNOut.insert(qn);
		}
		else {
			DySharpOut.insert(qt);
			DySharpNOut.insert(qn);
		}
	}
	Eigen::MatrixXd
		testout = DyTestOut.matXd(),
		testoutn = DyTestNOut.matXd(),
		sharpout = DySharpOut.matXd(),
		sharpnout = DySharpNOut.matXd();
	UTILS::outputXYZ(testout, testoutn, "test");
	UTILS::outputXYZ(sharpout, sharpnout, "test_sharp");


	if (VERBOSE_ON) {
		VERBOSE.append("sample P", PSize);
		VERBOSE.append("sample Q", int(PreQVerts.rows()));
	}
}


void OffsetCrust::preprocessing_sharp(
	Eigen::MatrixXd& Vertices,
	Eigen::MatrixXi& Faces,
	Eigen::MatrixXd& OutAllVertices,
	Eigen::MatrixXd& OutQNormals,
	Eigen::VectorXi& blue_validFlags,
	double q_orient
) {
	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}

	//1. presets
	igl::per_face_normals(Vertices, Faces, OrgFaceNormals);
	PQP = new ProjPoint(OrgVertices, OrgFaces, OrgFaceNormals, __IS_INSIDE);
	if (__IS_INSIDE) OrgFaceNormals = -OrgFaceNormals;


	// Heuristic to  determine radius from desired number 
	Eigen::MatrixXd P_blue, N_blue, Centroids;
	Eigen::MatrixXd outV, outVN;
	Eigen::VectorXi FI;
	std::map<int, std::vector<int>> f2blueP;

	if (__useCentroid) {
		Centroids.conservativeResize(Faces.rows(), 3);
		for (int fid = 0; fid < Faces.rows(); fid++) {
			int v0 = Faces(fid, 0),
				v1 = Faces(fid, 1),
				v2 = Faces(fid, 2);
			Eigen::VectorXd
				v0_p = Vertices.row(v0),
				v1_p = Vertices.row(v1),
				v2_p = Vertices.row(v2);
			Eigen::VectorXd cen = (v0_p + v1_p + v2_p) / 3;
			Centroids.row(fid) << cen(0), cen(1), cen(2);

			qi2pi[fid] = fid;
			pi2qi[fid].push_back(fid);
			f2blueP[fid].push_back(fid);
		}
		outV = Centroids, outVN = OrgFaceNormals;
		std::cout << "[Centroids] " << Centroids.rows() << std::endl;
	}
	else {
		// Heuristic to  determine radius from desired number 
		const double r = [&](const int n) {
			Eigen::VectorXd A;
			igl::doublearea(OrgVertices, OrgFaces, A);
			return sqrt(((A.sum() * 0.5 / (n * 0.6162910373)) / igl::PI));
			}(blue_noise_desired);

			printf("blue noise radius: %g\n", r);
			Eigen::MatrixXd B;
			igl::blue_noise(OrgVertices, OrgFaces, r, B, FI, P_blue);
			N_blue.conservativeResize(P_blue.rows(), P_blue.cols());
			for (int i = 0; i < FI.size(); i++) {
				int id = FI(i);
				N_blue.row(i) = OrgFaceNormals.row(id);
				qi2pi[i] = i;
				pi2qi[i].push_back(i);
				f2blueP[id].push_back(i);
			}
			outV = P_blue, outVN = N_blue;
			std::cout << "[Blue Noise] " << P_blue.rows() << std::endl;
	}

	Eigen::MatrixXd outVQ = outV + q_orient * q_offset * outVN;
	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd>
		DyPVertices(outV),
		DyQVertices(outVQ),
		DyQNormals(outVN);

	blue_validFlags = Eigen::VectorXi::Ones(outV.rows());

	//2. get all edges
	EdgePool edges(ratio);
	std::map<int, std::set<EdgeInfo*>> v2edges;
	std::map<int, std::vector<int>> v2faces;
	for (int fid = 0; fid < Faces.rows(); fid++) {
		for (int i = 0; i < 3; i++) {
			int v0 = Faces(fid, i),
				v1 = Faces(fid, (i + 1) % 3);
			EdgeInfo* e = edges.getEdge(v0, v1);
			if (!e->isSetVal) {
				Eigen::VectorXd
					vv0 = Vertices.row(v0), vv1 = Vertices.row(v1);
				e->setEdgeInfo(v0, v1, vv0, vv1);
			}
			e->adj_faces.push_back(fid);
			v2faces[v0].push_back(fid);
			v2edges[v0].insert(e);
		}
	}
	std::cout << "[ratio=] " << ratio << std::endl;
	std::cout << "[q_offset=] " << q_offset << std::endl;

	std::set<int> v_corners;
	edges.forEachEdgeInfo([&](EdgeInfo* edge) {
		if (edge->adj_faces.size() != 2) return;

		int f1 = edge->adj_faces[0], f2 = edge->adj_faces[1];
		Eigen::VectorXd
			n1 = OrgFaceNormals.row(f1), n2 = OrgFaceNormals.row(f2);
		double cosine = n1.transpose() * n2;
		cosine = std::clamp(cosine, -1.0, 1.0);
		if (cosine >= diheralBar) //not feature
			return;

		int v1 = edge->v0, v2 = edge->v1,
			f1v0, f2v0;
		Eigen::VectorXd
			p1 = edge->vv0, p2 = edge->vv1,
			uniq_f2v(3), uniq_f1v(3);

		//classify: convex / concave
		for (int i = 0; i < 3; i++) {
			int f1v = Faces(f1, i);
			if (f1v != v1 && f1v != v2) {
				f1v0 = f1v;
				uniq_f1v = Vertices.row(f1v);
			}

			int f2v = Faces(f2, i);
			if (f2v != v1 && f2v != v2) {
				f2v0 = f2v;
				uniq_f2v = Vertices.row(f2v);
			}
		}
		Eigen::VectorXd jud_vec = uniq_f2v - p1;
		double jud_val = jud_vec.transpose() * n1;
		edge->isConvex = (jud_val >= 0) ? 0 : 1;

		//exclude some blue noise sample points.
		int num_points = (edge->corner_ratio > q_offset) ? 1. / (edge->corner_ratio) : 0;
		if (ratio > q_offset) {
			EdgeInfo
				* f1e1 = edges.getEdge(f1v0, v1), * f1e2 = edges.getEdge(f1v0, v2),
				* f2e1 = edges.getEdge(f2v0, v1), * f2e2 = edges.getEdge(f2v0, v2);
			Eigen::VectorXd
				f1sp1, f1sp2, f2sp1, f2sp2;
			f1e1->getEdgeRatio(v1, f1v0, f1sp1);
			f1e2->getEdgeRatio(v2, f1v0, f1sp2);
			f2e1->getEdgeRatio(v1, f2v0, f2sp1);
			f2e2->getEdgeRatio(v2, f2v0, f2sp2);

			for (auto bluep : f2blueP[f1]) {
				Eigen::VectorXd bluevp = outV.row(bluep);
				bool flag = UTILS::isPointInTriangle(bluevp, f1sp1, f1sp2, uniq_f1v);
				blue_validFlags(bluep) = blue_validFlags(bluep) && flag ? 1 : 0;
			}
			for (auto bluep : f2blueP[f2]) {
				Eigen::VectorXd bluevp = outV.row(bluep);
				bool flag = UTILS::isPointInTriangle(bluevp, f2sp1, f2sp2, uniq_f2v);
				blue_validFlags(bluep) = blue_validFlags(bluep) && flag ? 1 : 0;
			}


			K::FT point_step = 1. / num_points;

			for (K::FT beta = point_step; ratio > q_offset && beta < 1 - q_offset; beta += point_step) {
				Eigen::VectorXd
					f1sp = CGAL::to_double(beta) * f1sp1 + CGAL::to_double(1 - beta) * f1sp2,
					f2sp = CGAL::to_double(beta) * f2sp1 + CGAL::to_double(1 - beta) * f2sp2;
				int spid1 = DyPVertices.insert(f1sp),
					spid2 = DyPVertices.insert(f2sp);
				DyQNormals.insert(n1);
				DyQNormals.insert(n2);
				Eigen::VectorXd
					f1sq = f1sp + q_offset * n1,
					f2sq = f2sp + q_offset * n2;
				int sqid1 = DyQVertices.insert(f1sq),
					sqid2 = DyQVertices.insert(f2sq);
				qi2pi[sqid1] = spid1;
				pi2qi[spid1].push_back(sqid1);  //1 v n
				qi2pi[sqid2] = spid2;
				pi2qi[spid2].push_back(sqid2);  //1 v n
			}

			//}
		}


		v_corners.insert(v1);
		v_corners.insert(v2);

		
		K::FT point_step = 1. / (num_points + 2);
		for (K::FT beta = point_step; beta < 1 - q_offset; beta += point_step) {
			Eigen::VectorXd beta_p = CGAL::to_double(beta) * p1 + CGAL::to_double(1 - beta) * p2;
			int pi_rowid;
			bool isPushed = false;

			double
				alpha_step = 1;
			double Omega = std::acos(cosine),
				sOmega = std::sin(Omega);
			for (double alpha = 0; alpha <= 1+q_offset; alpha += alpha_step) {
				Eigen::VectorXd alpha_n
					= std::sin((1 - alpha) * Omega) / sOmega * n1
					+ std::sin(alpha * Omega) / sOmega * n2;

				if (!isPushed) {
					pi_rowid = DyPVertices.insert(beta_p);
					isPushed = true;
				}
				alpha_n.normalize();
				DyQNormals.insert(alpha_n);

				Eigen::VectorXd beta_q = beta_p + q_offset * alpha_n;
				int qi_rowid = DyQVertices.insert(beta_q);
				qi2pi[qi_rowid] = pi_rowid;
				pi2qi[pi_rowid].push_back(qi_rowid);  //1 v n
			}
		}

		});


	K::FT point_step = 1. / num_corner_points;
	//3.1 handle all corners
	for (auto corner : v_corners) {
		Eigen::VectorXd beta_p = Vertices.row(corner);
		std::set<EdgeInfo*> adj_edges = v2edges[corner];

		/*1v1*/
		std::vector<int> adjF = v2faces[corner];
		for (int fi : adjF) {
			int ep1, ep2, vid = corner;
			for (int i = 0; i < 3; i++) {
				int v = Faces(fi, i);
				if (v == vid) {
					ep1 = Faces(fi, (i + 1) % 3);
					ep2 = Faces(fi, (i + 2) % 3);
					break;
				}
			}
			Eigen::VectorXd
				uniq_fv = Vertices.row(vid),
				n = OrgFaceNormals.row(fi);


			if (ratio > q_offset) {
				
				//exclude some blue noise sample points.
				EdgeInfo
					* f1e1 = edges.getEdge(vid, ep1), * f1e2 = edges.getEdge(vid, ep2);
				Eigen::VectorXd
					f1sp1, f1sp2;
				f1e1->getEdgeRatio(vid, ep1, f1sp1);
				f1e2->getEdgeRatio(vid, ep2, f1sp2);

				for (auto bluep : f2blueP[fi]) {
					Eigen::VectorXd bluevp = outV.row(bluep);
					bool flag = UTILS::isPointInTriangle(bluevp, f1sp1, f1sp2, uniq_fv);
					blue_validFlags(bluep) = blue_validFlags(bluep) && !flag ? 1 : 0;
				}


				for (K::FT beta = point_step; beta < 1 - q_offset; beta += point_step) {
					Eigen::VectorXd
						f1sp = CGAL::to_double(beta) * f1sp1 + CGAL::to_double(1 - beta) * f1sp2;
					int spid1 = DyPVertices.insert(f1sp);
					DyQNormals.insert(n);
					Eigen::VectorXd
						f1sq = f1sp + q_offset * n;
					int sqid1 = DyQVertices.insert(f1sq);
					qi2pi[sqid1] = spid1;
					pi2qi[spid1].push_back(sqid1);  //1 v n
				}

			}

		}

		std::vector<EdgeInfo*> convex_edges;
		for (auto edge : adj_edges) {
			if (edge->isConvex == 1)
				convex_edges.push_back(edge);
		}
		
		if (convex_edges.size() < 2)
			continue;

		// add corner
		int pi_rowid = DyPVertices.insert(beta_p);
		for (int fi : adjF) {
			Eigen::VectorXd
				n = OrgFaceNormals.row(fi);
			Eigen::VectorXd qi_n = beta_p + q_orient * q_offset * n;
			DyQNormals.insert(n);
			int qi_rowid = DyQVertices.insert(qi_n);
			qi2pi[qi_rowid] = pi_rowid;
			pi2qi[pi_rowid].push_back(qi_rowid);  //1 v n
		}

	}


	if (VERBOSE_ON) {
		VERBOSE.endTimer("time---preproc");
	}

	//4. final results
	PrePVerts = DyPVertices.matXd();
	PreQVerts = DyQVertices.matXd();
	PSize = DyPVertices.rowcnt();
	DyPVertices.append(PreQVerts);
	OutAllVertices = DyPVertices.matXd();
	OutQNormals = DyQNormals.matXd();
	UTILS::outputXYZ(PreQVerts, OutQNormals, "test");

	if (VERBOSE_ON) {
		VERBOSE.append("sample P", PSize);
		VERBOSE.append("sample Q", int(PreQVerts.rows()));
	}
}





void OffsetCrust::core_main(
	double d, 
	bool is_inside, 
	std::string outOBJFile, 
	std::string outMetaFile,
	bool relative_dis,
	bool outPoly,
	bool useCentroid,
	bool do_sharp,
	double diheral
) {
	__outputPoly = outPoly;
	__useCentroid = useCentroid;
	diheralBar = std::cos(diheral);

	if (relative_dis) 
		OFFSET = d * diagonal_length;
	else 
		OFFSET = d;
	__IS_INSIDE = is_inside;
	VERBOSE.append("dis", d);
	VERBOSE.append("real_dis", OFFSET);
	
	Eigen::MatrixXd AllVertices, QNormals;
	Eigen::VectorXi blue_validFlags;
	if(do_sharp)
		preprocessing_sharp(OrgVertices, OrgFaces, AllVertices, QNormals, blue_validFlags);
	else
		preprocessing_with_sphere(OrgVertices, OrgFaces, AllVertices, QNormals, blue_validFlags);
	double w_diff = q_offset * q_offset - 2 * q_offset * OFFSET;
	compute_power_diagram(AllVertices, QNormals, blue_validFlags, w_diff, outOBJFile);

	VERBOSE.toFile(outMetaFile);


	return;
}



