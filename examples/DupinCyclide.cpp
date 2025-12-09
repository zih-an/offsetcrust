#include <DupinCyclide.hpp>

double DupinCyclide::implicit_func(Eigen::Vector3d& pt) {
	double x = pt(0), y = pt(1), z = pt(2);
	double term = x * x + y * y + z * z + b * b - d * d;
	return term * term - 4 * (a * x - c * d) * (a * x - c * d) - 4 * b * b * y * y;
}

Eigen::Vector3d DupinCyclide::grad(Eigen::Vector3d& pt) {
	double x = pt(0), y = pt(1), z = pt(2);
	double term = x * x + y * y + z * z + b * b - d * d;

	double dx = 4 * x * term - 8 * a * (a * x - c * d);
	double dy = 4 * y * term - 8 * b * b * y;
	double dz = 4 * z * term;
	return Eigen::Vector3d(dx, dy, dz);
}


Eigen::Vector3d DupinCyclide::parametric_func(double u, double v) {
	double denominator = a - c * std::cos(u) * std::cos(v);
	if (std::abs(denominator) < 1e-6) { // ���������
		return Eigen::Vector3d(0, 0, 0); // ���׳��쳣
	}
	double x = (d * (c - a * std::cos(u) * std::cos(v)) + b * b * std::cos(u)) / denominator;
	double y = (b * std::sin(u) * (a - d * std::cos(v))) / denominator;
	double z = (b * std::sin(v) * (c * std::cos(u) - d)) / denominator;
	return Eigen::Vector3d(x, y, z);
}



double DupinCyclide::offset_radius(double u, double v) {
	if (isConstant) return (R_max + R_min) / 2;
	return R_min + (R_max - R_min) * 0.5 * (1 + std::cos(u)); 
}


Eigen::Vector3d DupinCyclide::surface_gradient(double u, double v) {
	// Evaluate S, Su, Sv
	auto S = parametric_func(u, v);

	double h = 1e-6; 
	Eigen::Vector3d Su = (parametric_func(u + h, v) - parametric_func(u - h, v)) / (2 * h);
	Eigen::Vector3d Sv = (parametric_func(u, v + h) - parametric_func(u, v - h)) / (2 * h);

	// Evaluate Ru, Rv
	double Ru = (offset_radius(u + h, v) - offset_radius(u - h, v)) / (2 * h);
	double Rv = (offset_radius(u, v + h) - offset_radius(u, v - h)) / (2 * h);


	// Compute metric tensor
	double E = Su.dot(Su);
	double F = Su.dot(Sv);
	double G = Sv.dot(Sv);

	Eigen::Matrix2d g;
	g << E, F,
		F, G;

	if (g.determinant() < 1e-12) {
		g += Eigen::Matrix2d::Identity() * 1e-6; 
	}

	Eigen::Vector2d grad_param(Ru, Rv);
	Eigen::Vector2d coeffs = g.inverse() * grad_param;


	Eigen::Vector3d grad3D = coeffs(0) * Su + coeffs(1) * Sv;
	return grad3D;
}




void DupinCyclide::mc(int res) {

	double L = std::max(
		a + d - c,
		std::abs(-a - (d + c))
	);
	Eigen::RowVector3d
		padding = Eigen::RowVector3d(0.1, 0.1, 0.1),
		min_corner = Eigen::RowVector3d(-L, -L, -L) - padding,
		max_corner = Eigen::RowVector3d(L, L, L) + padding;

	int nx = res + 1;
	int ny = res + 1;
	int nz = res + 1;


	Eigen::MatrixXd GV;
	igl::grid(Eigen::RowVector3i(nx, ny, nz), GV);
	Eigen::VectorXd Gf(GV.rows());
	igl::parallel_for(GV.rows(), [&](const int i)
		{
			GV.row(i).array() *= (max_corner - min_corner).array();
			GV.row(i) += min_corner;
			Eigen::Vector3d pt = GV.row(i);
			Gf(i) = implicit_func(pt); 
		}, 1000ul);

	Eigen::MatrixXd mcV;
	Eigen::MatrixXi mcF;
	igl::marching_cubes(Gf, GV, nx, ny, nz, 0, mcV, mcF);

	igl::writeOBJ(UTILS::data_dir + "cyclide_mc.obj", mcV, mcF);

}



void DupinCyclide::random_sample_para(int N) {
	int seed = 42;
	std::mt19937 rng(seed); 
	std::uniform_real_distribution<double> dist(0, 2 * M_PI);

	int num = N;
	V.conservativeResize(num, 3);
	VN.conservativeResize(num, 3);
	RGrad.conservativeResize(num, 3);
	Radius.conservativeResize(num);

	for (int i = 0; i < N; ++i) {
		double 
			u = dist(rng),
			v = dist(rng);

		Eigen::Vector3d pt = parametric_func(u, v);
		//vertices.push_back(pt);
		V.row(i) = pt;
		VN.row(i) = grad(pt).normalized();
		RGrad.row(i) = surface_gradient(u, v);
		Radius(i) = offset_radius(u, v);
	}

	
	UTILS::outputXYZ(V, VN, "cyclide_surf_normal");
	UTILS::outputXYZ(V, RGrad, "cyclide_normal");

}



void DupinCyclide::para(
	int u_steps,
	int v_steps
) {
	std::vector<Eigen::Vector3d> vertices;
	std::vector<Eigen::Vector3i> faces;

	int num = u_steps * v_steps;
	V.conservativeResize(num, 3);
	VN.conservativeResize(num, 3);
	RGrad.conservativeResize(num, 3);
	Radius.conservativeResize(num);

	for (int i = 0; i < u_steps; ++i) {
		double u = 2 * M_PI * i / (u_steps - 1);
		for (int j = 0; j < v_steps; ++j) {
			int id = i * v_steps + j;
			double v = 2 * M_PI * j / (v_steps - 1);
			Eigen::Vector3d pt = parametric_func(u, v);

			V.row(id) = pt;
			VN.row(id) = grad(pt).normalized();
			RGrad.row(id) = surface_gradient(u, v);
			Radius(id) = offset_radius(u, v);
		}
	}




	UTILS::outputXYZ(V, RGrad, "cyclide_normal");

	for (int i = 0; i < u_steps - 1; ++i) { 
		for (int j = 0; j < v_steps - 1; ++j) {
			int idx = i * v_steps + j;
			int next_i = (i + 1); 
			int next_j = (j + 1) % v_steps;

			int idx_next_i = next_i * v_steps + j;
			int idx_next_j = i * v_steps + next_j;
			int idx_next_ij = next_i * v_steps + next_j;


			faces.push_back(Eigen::Vector3i(idx, idx_next_i, idx_next_ij));
			faces.push_back(Eigen::Vector3i(idx, idx_next_ij, idx_next_j));
		}
	}

	F.conservativeResize(faces.size(), 3);
	for (int i = 0; i < faces.size(); ++i) {
		F.row(i) = faces[i].cast<int>();
	}

	igl::writeOBJ(UTILS::out_dir + std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "dupinc_cyclide.obj", V, F);
}



void DupinCyclide::sample(
	Eigen::MatrixXd& outV,
	Eigen::MatrixXd& outVN,
	Eigen::VectorXd& outR,
	Eigen::MatrixXd& outRGrad,
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
	outR.conservativeResize(outV.rows());
	outRGrad.conservativeResize(outV.rows(), outV.cols());

	for (int i = 0; i < FI.size(); i++) {
		int fid = FI(i);
		outVN.row(i) = FN.row(fid);
		Eigen::VectorXd pt = outV.row(i);
		Eigen::Vector3d grad;
		double r = get_F(pt, fid, grad);
		outR(i) = r;
		outRGrad.row(i) = grad.transpose();
	}

}



std::pair<double, double> DupinCyclide::estimate_initial_uv(const Eigen::Vector3d& p) {
	int sample_u = 20;
	int sample_v = 20;
	double min_dist = std::numeric_limits<double>::max();
	double best_u = 0.0, best_v = 0.0;

	for (int i = 0; i < sample_u; ++i) {
		double u = 2 * M_PI * i / sample_u;
		for (int j = 0; j < sample_v; ++j) {
			double v = 2 * M_PI * j / sample_v;
			Eigen::Vector3d q = parametric_func(u, v);
			double dist = (q - p).squaredNorm(); 
			if (dist < min_dist) {
				min_dist = dist;
				best_u = u;
				best_v = v;
			}
		}
	}
	return { best_u, best_v };
}



double DupinCyclide::pt2surfDis(Eigen::Vector3d& p) {
	auto para = estimate_initial_uv(p);
	double u = para.first, v = para.second;
	const int max_iters = 100;
	const double eps = 1e-6;
	const double h = 1e-6;

	for (int iter = 0; iter < max_iters; ++iter) {
		Eigen::Vector3d x = parametric_func(u, v);
		Eigen::Vector3d r = x - p;

		// Numerical derivatives
		Eigen::Vector3d xu = (parametric_func(u + h, v) - parametric_func(u - h, v)) / (2 * h);
		Eigen::Vector3d xv = (parametric_func(u, v + h) - parametric_func(u, v - h)) / (2 * h);

		Eigen::Matrix<double, 3, 2> J;
		J.col(0) = xu;
		J.col(1) = xv;

		Eigen::Matrix2d H = J.transpose() * J;
		H += 1e-8 * Eigen::Matrix2d::Identity(); // Damping
		Eigen::Vector2d delta = H.ldlt().solve(-J.transpose() * r);

		u += delta(0);
		v += delta(1);

		// Wrap u,v into [0,2pi]
		u = std::fmod(u, 2 * M_PI);
		if (u < 0) u += 2 * M_PI;
		v = std::fmod(v, 2 * M_PI);
		if (v < 0) v += 2 * M_PI;

		if (delta.norm() < eps)
			break;
	}

	Eigen::Vector3d q = parametric_func(u, v);
	double dis = (p - q).norm();
	std::cout << dis << std::endl;
	return dis;
}





void DupinCyclide::mat(int u_steps, int v_steps) {
	Eigen::MatrixXd SPHV, SPHVN;
	Eigen::MatrixXi SPHF;
	igl::readOBJ(UTILS::sphere_file, SPHV, SPHF);
	UTILS::normalization_2(SPHV);
	igl::per_vertex_normals(SPHV, SPHF, SPHVN);


	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd> verts(3),
		matV(3), matL(3), matQQ(3), matQQL(3);
	DynamicRowMatrix<Eigen::MatrixXi, Eigen::VectorXi> faces(3);

	std::vector<std::vector<Eigen::VectorXd>> PNormals;
	Eigen::VectorXd Radius(u_steps);
	for (int i = 0; i < u_steps; ++i) {
		double u = 2 * M_PI * i / (u_steps );

		double 
			x = a * std::cos(u),
			y = b * std::sin(u);
		Eigen::VectorXd pt = Eigen::Vector3d(x, y, 0);
		//double dis = pt2surfDis(pt);
		double dis = d - c * std::cos(u);
		Radius(i) = dis;

		Eigen::VectorXd ones = Eigen::VectorXd::Ones(SPHV.rows());
		Eigen::MatrixXd trans_SPHV = (SPHV*dis + ones * pt.transpose());
		Eigen::MatrixXi trans_SPHF = SPHF +
			i * SPHV.rows() * Eigen::MatrixXi::Ones(SPHF.rows(), SPHF.cols());

		verts.append(trans_SPHV);
		faces.append(trans_SPHF);
		matV.insert(pt);

		std::vector<Eigen::VectorXd> qns;
		
		for (int vv = 0; vv < v_steps; vv++) {
			double v = 2 * M_PI * vv / (v_steps );
			Eigen::Vector3d qt = parametric_func(u, v);
			Eigen::VectorXd qn = (qt - pt).normalized(),
				qtt = qt;
			qns.push_back(qn);
			matQQ.insert(qtt);
			
			matQQL.insert(pt);
			matQQL.insert(qtt);
		}
		PNormals.push_back(qns);

	}
	Eigen::MatrixXd 
		matVV = matV.matXd(),
		matQQV = matQQ.matXd(),
		matQQLL = matQQL.matXd();
	UTILS::outputXYZ(matQQV, std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "cyclide_mat_q");
	UTILS::outputOBJ_with_lines(matQQLL, std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "cyclide_ll");
	OffsetCrust crust;
	crust.core_main(matVV, PNormals, Radius, UTILS::out_dir + std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "cyclide_mat");



	/*output*/
	for (int i = 0; i < matVV.rows(); i++) {
		Eigen::VectorXd
			pt1 = matVV.row(i),
			pt2 = matVV.row((i + 1) % matVV.rows());
		matL.insert(pt1);
		matL.insert(pt2);
	}
	Eigen::MatrixXd matLL = matL.matXd();
	UTILS::outputOBJ_with_lines(matVV, matLL, std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "mat_directrix.obj");
	std::cout << UTILS::out_dir + std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "mat_directrix.obj" << std::endl;

	//para(100, 100);
	mc();
	Eigen::MatrixXd meshV = verts.matXd();
	Eigen::MatrixXi meshF = faces.matXd();
	igl::writeOBJ(UTILS::out_dir + std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_" + "mat_balls.obj", meshV, meshF);

}


//linear
double DupinCyclide::get_F(Eigen::VectorXd& pt, int fid, Eigen::Vector3d& grad) {
	int v0 = F(fid, 0),  //i
		v1 = F(fid, 1),  //j
		v2 = F(fid, 2);  //k
	Eigen::Vector3d
		vv0 = V.row(v0),
		vv1 = V.row(v1),
		vv2 = V.row(v2),
		ppt = pt;
	std::tuple<double, double, double> res = UTILS::tri_interpolation(vv0, vv1, vv2, fid, ppt);
	double val = std::get<0>(res) * Radius(v0) + std::get<1>(res) * Radius(v1) + std::get<2>(res) * Radius(v2);

	

	grad = UTILS::tri_grad(vv0, vv1, vv2, Radius(v0), Radius(v1), Radius(v2));
	

	return val;
}




void DupinCyclide::smooth_main() {
	OffsetCrust crust;
	isConstant = true;
	para(30, 15);

	crust.core_main(V, VN, Radius, RGrad, UTILS::out_dir + "smooth_constant", isConstant);
	UTILS::outputOBJ_with_normal(V, VN, "cyclide_sample", 0.2);

	isConstant = false;

	random_sample_para(50000);
	crust.core_main(V, VN, Radius, RGrad, UTILS::out_dir + "smooth_variable", true);
}


void DupinCyclide::mesh_main() {
	int n = 30 * 30;
	OffsetCrust crust;

	isConstant = true;
	para();
	Eigen::MatrixXd smplV, smplVN, smplGrad;
	Eigen::VectorXd smplR;
	sample(smplV, smplVN, smplR, smplGrad, n);
	crust.core_main(smplV, smplVN, smplR, smplGrad, 
		UTILS::out_dir + "mesh_constant", isConstant);

	UTILS::outputOBJ_with_normal(smplV, smplVN, "cyclide_mesh_sample", 0.2);

	isConstant = false;
	para();
	sample(smplV, smplVN, smplR, smplGrad, n);
	crust.core_main(smplV, smplVN, smplR, smplGrad,
		UTILS::out_dir + "mesh_variable", true);

	std::ofstream ifs(UTILS::out_dir + "mesh_variable_scalers.txt");
	for (int i = 0; i < V.rows(); i++)
		ifs << i << " " << Radius(i) << "\n";
	std::cout << UTILS::out_dir + "mesh_variable_scalers.txt" << std::endl;
	ifs.close();
}

