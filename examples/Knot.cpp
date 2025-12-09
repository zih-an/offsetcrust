#include <Knot.hpp>


Eigen::Vector3d Knot::curve_para(double u) {
	double
		x = (2 + std::cos(2 * u)) * std::cos(3 * u),
		y = (2 + std::cos(2 * u)) * std::sin(3 * u),
		z = std::sin(4 * u);

	return Eigen::Vector3d(x, y, z);
}


Eigen::Vector3d Knot::curve_deriv(double u) {
	double h = 1e-5;
	return (curve_para(u + h) - curve_para(u - h)) / (2 * h);
}

double Knot::radius_deriv(double u) {
	double h = 1e-5;
	return (radius_function(u + h) - radius_function(u - h)) / (2 * h);
}

Eigen::Vector3d Knot::radius_grad(double u) {
	Eigen::Vector3d dxdu = curve_deriv(u);
	double drdu = radius_deriv(u);
	double denom = dxdu.squaredNorm();  // |dx/du|^2
	return (drdu / denom) * dxdu;
}




void Knot::mat(int u_steps, int v_steps) {



	DynamicRowMatrix<Eigen::MatrixXd, Eigen::VectorXd> verts(3),
		matV(3), matL(3), matQQ(3), matQQL(3);
	DynamicRowMatrix<Eigen::MatrixXi, Eigen::VectorXi> faces(3);

	std::vector<std::vector<Eigen::VectorXd>> PNormals;
	Eigen::VectorXd Radius(u_steps);
	for (int i = 0; i < u_steps; ++i) {
		double u = 2 * M_PI * i / (u_steps);

		Eigen::VectorXd pt = curve_para(u);


		double dis = radius_function(u);
		Radius(i) = dis;

		matV.insert(pt);



	
		double dx_du = -(2 + cos(2 * u)) * 3 * sin(3 * u) - 2 * sin(2 * u) * cos(3 * u);
		double dy_du = (2 + cos(2 * u)) * 3 * cos(3 * u) - 2 * sin(2 * u) * sin(3 * u);
		double dz_du = 4 * cos(4 * u);

		Eigen::Vector3d T(dx_du, dy_du, dz_du);
		T.normalize();
		Eigen::Vector3d ref(0, 0, 1);
		if (fabs(T.dot(ref)) > 0.99) {
			ref = Eigen::Vector3d(1, 0, 0);
		}
		Eigen::Vector3d N = T.cross(ref).normalized();
		Eigen::Vector3d B = T.cross(N).normalized();
		std::vector<Eigen::VectorXd> qns;

	
		for (int vv = 0; vv < v_steps; vv++) {
			double v = 2 * M_PI * vv / (v_steps);
			Eigen::VectorXd
				qn = (std::cos(v) * N + std::sin(v) * B).normalized(),
				qtt = pt + dis * qn;
			qns.push_back(qn);
			matQQ.insert(qtt);

			matQQL.insert(pt);
			matQQL.insert(qtt);
		}

		Eigen::Vector3d grad = radius_grad(u);
		OffsetCrust::rotateN(qns, grad);
		PNormals.push_back(qns);

	}
	std::string prefix = std::to_string(rfid) + "_" + std::to_string(u_steps) + "_" + std::to_string(v_steps) + "_";
	Eigen::MatrixXd
		matVV = matV.matXd(),
		matQQV = matQQ.matXd(),
		matQQLL = matQQL.matXd();
	UTILS::outputXYZ(matQQV, prefix + "knot_mat_q");
	UTILS::outputOBJ_with_lines(matQQLL, prefix + "knot_ll");
	OffsetCrust crust;
	crust.core_main(matVV, PNormals, Radius, 
		UTILS::out_dir + prefix + "knot_mat");




	for (int i = 0; i < matVV.rows(); i++) {
		Eigen::VectorXd
			pt1 = matVV.row(i),
			pt2 = matVV.row((i + 1) % matVV.rows());
		matL.insert(pt1);
		matL.insert(pt2);
	}
	Eigen::MatrixXd matLL = matL.matXd();
	UTILS::outputOBJ_with_lines(matVV, matLL, prefix + "knot_mat_directrix.obj");
	


}