#pragma once

#include <offset_crust_neq.hpp>
#include <igl/sparse_voxel_grid.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <cmath>
#include <random>

using namespace Offset3D;

#pragma once
#include <functional>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

class RadiusFunctionFactory {
public:
    using RadiusFn = std::function<double(double)>;

    static RadiusFn get_by_index(int i) {
        const auto& functions = get_all();
        return functions[i % functions.size()];
    }

    static RadiusFn get_by_name(const std::string& name) {
        if (name == "fn1") return fn1;
        if (name == "fn2") return fn2;
        if (name == "fn3") return fn3;
        if (name == "fn4") return fn4;
        if (name == "fn5") return fn5;
        if (name == "fn6") return fn6;
        return fn1; // default
    }

    static const std::vector<RadiusFn>& get_all() {
        static std::vector<RadiusFn> functions = { fn1, fn2, fn3, fn4, fn5, fn6 };
        return functions;
    }

private:

    static double fn1(double u) {
        return 0.45;
    }

    static double fn2(double u) {
        return 0.4 + 0.3 * std::cos(u); // ¡Ê [0.1, 0.5]
    }

    static double fn3(double u) {
        return 0.4 + 0.2 * std::cos(3 * u); // ¡Ê [0.1, 0.5]
    }

    static double fn4(double u) {
        return 0.4 + 0.3 * std::exp(-5 * std::pow(std::sin(u), 2)); // ¡Ê ~[0.3, 0.5]
    }

    static double fn5(double u) {
        double base = 0.4 + 0.15 * std::cos(2 * u);
        double noise = 0.1 * std::sin(13 * u + std::cos(u)) + 0.05 * std::sin(20 * u); // ¸üÇ¿ÁÒµÄÔëÉù
        return base + noise; // ¡Ê [0.1, 0.5]
    }

    static double fn6(double u) {
        return 0.3 + 0.1 * std::cos(2 * u) + 0.1 * std::cos(8 * u); // ¡Ê [0.1, 0.5]
    }
};



class Knot {

private:
    Eigen::MatrixXd V, VN, RGrad;
    Eigen::VectorXd Radius;
    Eigen::MatrixXi F;

    bool isConstant = false;
    int rfid;

    Eigen::Vector3d curve_para(double u);
    Eigen::Vector3d radius_grad(double u);
    double radius_deriv(double u);
    Eigen::Vector3d curve_deriv(double u);

    std::function<double(double)> radius_function;

public:
    Knot(int rf = 0) {
        rfid = rf;
        radius_function = RadiusFunctionFactory::get_by_index(rfid);
    }


    void mat(int u_steps = 10, int v_steps = 10);
};


