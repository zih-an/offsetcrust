#include <MAT_recon.hpp>
using namespace Offset3D;


MATRecon::MATRecon(std::string filename) {

	if (VERBOSE_ON) {
		VERBOSE.startTimer();
	}

	read_mat(filename);

	if (VERBOSE_ON) {
		VERBOSE.endTimer("time-mat");
	}

	VERBOSE.toFile(filename + "_verbose");
}



void MATRecon::read_mat(std::string matfile) {
	std::ifstream infile(matfile);

	if (!infile.is_open()) {
		std::cerr << "Error: Cannot open file " << matfile << std::endl;
		return;
	}

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}

		std::istringstream iss(line);
		std::string type;
		iss >> type;

		if (type == "v") {
			double x, y, z, r;
			iss >> x >> y >> z >> r;
			MA_v.emplace_back(x, y, z);
			MA_r.push_back(r);
		}
		else if (type == "e") {
			int v1, v2;
			iss >> v1 >> v2;
			MA_e.push_back({ v1 , v2 });
		}
		else if (type == "f") {
			int v1, v2, v3;
			iss >> v1 >> v2 >> v3;
			MA_f.push_back({ v1 , v2 , v3 });
		}
	}

	infile.close();


	maV.conservativeResize(MA_v.size(), 3);
	__VDIS.conservativeResize(MA_v.size());
	for (size_t i = 0; i < MA_v.size(); ++i) {
		maV.row(i) << MA_v[i](0), MA_v[i](1), MA_v[i](2);
		__VDIS(i) = MA_r[i];
	}

	maF.conservativeResize(MA_f.size(), 3);
	for (size_t i = 0; i < MA_f.size(); ++i) {
		maF.row(i) = Eigen::Vector3i(MA_f[i][0], MA_f[i][1], MA_f[i][2]);
	}

	maE.conservativeResize(MA_e.size(), 2);
	for (size_t i = 0; i < MA_e.size(); ++i) {
		maE.row(i) = Eigen::Vector2i(MA_e[i][0], MA_e[i][1]);
	}


	Eigen::MatrixXd outP;
	std::vector<std::vector<Eigen::VectorXd>> outPNormals;
	Eigen::VectorXd outRadius;

	OffsetCrust crust;
	preprocessing_MA(
		maV, maF, maE, __VDIS,
		outP, outPNormals, outRadius
	);
	core_main(outP, outPNormals, outRadius, matfile + "_recon");

}


bool MATRecon::sampleTriangleBand(
	double desired_step,
	int fid,
	std::vector<PInfo>& output,
	double q_orient,
	std::vector<int>& blueids,
	Eigen::MatrixXd& P_blue,
	Eigen::VectorXi& blue_validFlags
) {
	//return true;
	Eigen::Vector3d
		v1 = OrgVertices.row(OrgFaces(fid, 0)),
		v2 = OrgVertices.row(OrgFaces(fid, 1)),
		v3 = OrgVertices.row(OrgFaces(fid, 2));
	Eigen::VectorXd normal = OrgFaceNormals.row(fid);

	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> tri_edges = {
		std::make_pair(v1, v2),
		std::make_pair(v2, v3),
		std::make_pair(v3, v1),
	};

	bool isValidTri = true;
	double b_step = desired_step;


	for (int i = 0; i < 3; i++) {
		std::pair<Eigen::Vector3d, Eigen::Vector3d>
			e1 = tri_edges[i],
			e2 = tri_edges[(i + 1) % 3],
			e3 = tri_edges[(i + 2) % 3];
		Eigen::Vector3d
			vec_e1 = e1.second - e1.first,
			vec_e2 = e2.first - e2.second,
			vec_e3 = e3.second - e3.first,
			n = vec_e2.cross(vec_e3);
		n.normalize();
		Eigen::Vector3d orth_e1 = n.cross(vec_e1).normalized();


		double
			area2 = vec_e2.cross(vec_e3).norm(),
			hh = area2 / vec_e1.norm(),
			ratio = (hh - b_step) / hh;

		Eigen::Vector3d
			pp1 = vec_e2 * ratio + e2.second,
			pp2 = vec_e3 * ratio + e3.first;

		double t_step = desired_step / (pp1 - pp2).norm();
		int n_step = std::ceil(1 / t_step);
		t_step = n_step <= 0 ? 2 : 1. / n_step;
		isValidTri = isValidTri && (t_step < 1);

		for (double t = t_step; t < 1 - eps; t += t_step) {
			Eigen::Vector3d
				p = (1.0 - t) * pp1 + t * pp2;
			Eigen::VectorXd
				vvp = p,
				vv1 = v1, vv2 = v2, vv3 = v3;
			bool isIn = UTILS::isPointInTriangle(p, vv1, vv2, vv3);
			if (!isIn) continue;
			Eigen::VectorXd
				p01 = p,
				q01 = p01 + q_orient * q_offset * normal;
			PInfo pp(p01, get_F(p01, fid), q01, normal);
			output.push_back(pp);
		}

		for (auto blueid : blueids) {
			Eigen::Vector3d
				bluevp = P_blue.row(blueid),
				e1_orient = vec_e1.normalized(),
				proj_e1 = e1.first + e1_orient * (bluevp - e1.first).dot(e1_orient),
				proj_vec = proj_e1 - bluevp;
			bool flag = proj_vec.norm() > b_step;
			blue_validFlags(blueid) = blue_validFlags(blueid) && flag ? 1 : 0;

		}
	}


	if (!isValidTri) {
		for (auto blueid : blueids) {
			Eigen::VectorXd bluevp = P_blue.row(blueid);
			blue_validFlags(blueid) = 0;
		}
		output.clear();
	}
	return isValidTri;
}



void MATRecon::preprocessing_MA(
	Eigen::MatrixXd& maV,
	Eigen::MatrixXi& maF,
	Eigen::MatrixXi& maE,
	Eigen::VectorXd VDIS,
	Eigen::MatrixXd& outP,
	std::vector<std::vector<Eigen::VectorXd>>& outPNormals,
	Eigen::VectorXd& outRadius
) {
	OrgVertices = maV;
	OrgFaces = maF;
	__VDIS = VDIS;

	Eigen::SparseMatrix<double> G;
	igl::grad(OrgVertices, OrgFaces, G);
	gradient = Eigen::Map<const Eigen::MatrixXd>((G * __VDIS).eval().data(), OrgFaces.rows(), 3);

	igl::per_face_normals(OrgVertices, OrgFaces, OrgFaceNormals);



	Eigen::MatrixXd P_blue;
	Eigen::VectorXi FI;
	std::map<int, std::vector<int>> f2blueP;
	{
		//int blue = blue_desired_low;
		int blue = 5000;
		double blue_radius_low = [&](const int n) {
			Eigen::VectorXd A;
			igl::doublearea(maV, maF, A);
			return sqrt(((A.sum() * 0.5 / (n * 0.6162910373)) / igl::PI));
			}
		(blue);
		double radius_lb = __VDIS.minCoeff() * 2 * std::sin(M_PI / (2 * seg_angle));
		blue_noise_radius = std::min(radius_lb, blue_radius_low);

		Eigen::MatrixXd B;
		igl::blue_noise(maV, maF, blue_noise_radius, B, FI, P_blue);
	}


	
	std::vector<PInfo> samples;
	//face
	for (int i = 0; i < P_blue.rows(); i++) {
		int id = FI(i);
		f2blueP[id].push_back(i);
		Eigen::VectorXd
			pt = P_blue.row(i),
			qn_pos = OrgFaceNormals.row(id),
			qn_neg = -qn_pos,
			qt = pt + q_offset * qn_pos,
			qt2 = pt + q_offset * qn_neg;
		Eigen::Vector3d grad = gradient.row(id);

		PInfo pp(pt, get_F(pt, id));
		pp.insert(qt, qn_pos);
		pp.insert(qt2, qn_neg);
		rotateN(pp.qn, grad);
		samples.push_back(pp);
	}
	//add band.
	Eigen::VectorXi blue_validFlags = Eigen::VectorXi::Ones(P_blue.rows());

	double desired_step = blue_noise_radius;
	for (int fid = 0; fid < OrgFaces.rows(); fid++) {
		std::vector<PInfo> edge_samples;
		std::vector<int> blueids = f2blueP[fid];
		sampleTriangleBand(
			desired_step, fid, edge_samples, 1,
			blueids, P_blue, blue_validFlags
		);
		samples.insert(samples.end(), edge_samples.begin(), edge_samples.end());
	}
	int intri_size = samples.size();


	//vert
	for (int i = 0; i < maV.rows(); i++) {
		Eigen::VectorXd pt = maV.row(i);
		samples.push_back(PInfo(pt, __VDIS(i)));
	}
	//edge
	for (int i = 0; i < maE.rows(); i++) {
		int pi1 = maE(i, 0),
			pi2 = maE(i, 1);
		Eigen::VectorXd
			p1 = OrgVertices.row(pi1),
			p2 = OrgVertices.row(pi2);

		double point_step = blue_noise_radius / (p1 - p2).norm();
		int n_step = std::ceil(1 / point_step);
		point_step = n_step <= 0 ? 2 : 1. / n_step;
		for (double beta = point_step; beta < 1 - eps; beta += point_step) {
			Eigen::VectorXd p = (1 - beta) * p1 + beta * p2;
			double rp = (1 - beta) * __VDIS(pi1) + beta * __VDIS(pi2);
			samples.push_back(PInfo(p, rp));
		}
	}


	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd>
		DyCloud(3), DyPVertices(3);
	std::map<int, int> cloud2smpl;
	for (int i = 0; i < samples.size(); i++) {
		PInfo pp = samples[i];
		int cid = DyCloud.insert(pp.pt);
		cloud2smpl[cid] = i;
	}
	Eigen::MatrixXd cloud = DyCloud.matXd();
	KNNSearch knn(cloud);
	double knn_radius = 1.5 * blue_noise_radius;
	double inter_t = 1;


	//1vN
	std::cout << "begin 1vn..." << std::endl;

	std::vector<double> PRadius;

#pragma omp parallel for 
	for (int i = 0; i < samples.size(); i++) {
		if (i < intri_size && i >= P_blue.rows()) continue;  //skip band points
		if (i < P_blue.rows() && blue_validFlags(i) == 0) continue;  //skip invalid p_blue


		bool tri_isValid = true;
		PInfo smplP = samples[i];
		if (smplP.qt.size() == 0) {
			Eigen::VectorXd pt = smplP.pt;
			double pwrdis_d = 0;  //on the sphere.

			//get knn
			Eigen::VectorXi nnID;
			Eigen::VectorXd nnDIS;
			knn.rSearch(pt, nnID, nnDIS, knn_radius);
			std::vector<Eigen::VectorXd> nn_itps;
			std::vector<double> rr_itps;

			for (int nn = 0; nn < nnID.rows(); nn++) {
				int id = cloud2smpl[nnID(nn)];
				Eigen::VectorXd nn_pt = samples[id].pt;
				if ((nn_pt - pt).norm() < eps) continue;
				double rr_itp =
					(1 - inter_t) * smplP.pr + inter_t * samples[id].pr;
				Eigen::VectorXd nn_itp =
					(1 - inter_t) * pt + inter_t * nn_pt;
				rr_itps.push_back(rr_itp);
				nn_itps.push_back(nn_itp);
			}
			//insert 1vn
			Eigen::VectorXd ones = Eigen::VectorXd::Ones(SPHV.rows());
			Eigen::MatrixXd trans_SPHV = (SPHV * smplP.pr + ones * pt.transpose());
			Eigen::MatrixXd trans_SPHV_eps = (SPHV * q_offset + ones * pt.transpose());
			tri_isValid = false;
			for (int svid = 0; svid < trans_SPHV.rows(); svid++) {
				Eigen::VectorXd
					sv = trans_SPHV.row(svid),
					alpha_n = SPHVN.row(svid);
				bool isValid = true;

				for (int nn = 0; nn < nn_itps.size(); nn++) {
					Eigen::VectorXd nn_pt = nn_itps[nn];
					double pwrdis_nn = (sv - nn_pt).norm() - rr_itps[nn];
					if (pwrdis_nn < pwrdis_d) {
						isValid = false;
						break;
					}
				}
				//insert
				tri_isValid = tri_isValid || isValid;
				Eigen::VectorXd ss_eps = trans_SPHV_eps.row(svid);
				smplP.insert(ss_eps, alpha_n);
			}
		}

		if (smplP.qt.size() != 0 && tri_isValid) {
#pragma omp critical 
			{
				int pid = DyPVertices.insert(smplP.pt);
				PRadius.push_back(smplP.pr);
				outPNormals.push_back(smplP.qn);
			}
		}
	}
	outP = DyPVertices.matXd();
	outRadius = Eigen::Map<Eigen::VectorXd>(PRadius.data(), PRadius.size());

}

