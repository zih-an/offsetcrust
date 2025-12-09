#include <iostream>
#include <offset_crust_neq.hpp>
#include <utils.hpp>
#include <filesystem>  // C++17 filesystem library
#include <string>
#include <CLI/CLI.hpp>
namespace fs = std::filesystem;


#ifdef CGAL_LINKED_WITH_TBB
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>
#endif // CGAL_LINKED_WITH_TBB


std::string remove_suffix(const std::string& filename, const std::string& suffix) {
	if (filename.size() >= suffix.size() &&
		filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
		return filename.substr(0, filename.size() - suffix.size());
	}
	return filename;
}


int main(int argc, char** argv) {
	CLI::App app{ "Command Line Parser Example" };

	struct {
		std::string
			input_file,
			distance_file = "";

		bool inside = false,
			relative_dis = false,
			outPoly = false,
			useCentroid = false,
			doSharp = false,
			runGalleryDis = false;

		double dis = 0.02,
			corner_ratio = 0.05,
			lambda = 0.01,
			diheral_plane = 0.2;

		int _num_angle = 10,
			blue_points = 70000,
			parallel_n = 33;
	} args;

	app.add_option("--input", args.input_file, "Path to the input file")->required();
	app.add_option("-d", args.dis, "offset distance (relative), default=2%");
	app.add_flag("--diag", args.relative_dis, "relative offset distance");
	app.add_flag("--inside", args.inside, "Set to true if -inside flag is provided");
	app.add_option("--sphere", UTILS::sphere_file, "discrete spherical surface file");
	app.add_option("--outdir", UTILS::out_dir, "discrete spherical surface file");
	app.add_option("--slerp", args._num_angle, "slerp number. default=10");
	app.add_option("--blue", args.blue_points, "desired number of blue noise sample. default=70000");
	app.add_option("--corner", args.corner_ratio, "protected non-smooth region. default=0.05");
	app.add_option("--lambda", args.lambda, "optimization parameter. default=0.01");
	app.add_flag("--poly", args.outPoly, "output polygon results");
	app.add_flag("--centroid", args.useCentroid, "use centroid instead of blue noise sample. (recommend for dense triangle mesh)");
	app.add_flag("--sharp", args.doSharp, "all sharp preserved offset");
	app.add_option("--dis", args.distance_file, "user defined distance file. Every line format: vid dis. NOTE: please consider the norm of the gradient.");
	app.add_flag("--gallery", args.runGalleryDis, "all sharp preserved offset");
	app.add_option("--diheral", args.diheral_plane, "give a threshold (radian) to get sharp features (default=0.2)");
	app.add_option("-t", args.parallel_n, "num of threads. default = 33");


	CLI11_PARSE(app, argc, argv);


#ifdef CGAL_LINKED_WITH_TBB
	// Get the default number of threads
	int num_threads = oneapi::tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	std::cout << "num_threads: " << num_threads << std::endl;
	oneapi::tbb::global_control global_limit(
		oneapi::tbb::global_control::max_allowed_parallelism,
		args.parallel_n
	);
	num_threads = oneapi::tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	std::cout << "num_threads: " << num_threads << std::endl;
#endif // CGAL_LINKED_WITH_TBB

	Eigen::initParallel();
	omp_set_num_threads(args.parallel_n);
	Eigen::setNbThreads(args.parallel_n);
	std::cout << "Eigen THREADS = " << Eigen::nbThreads() << std::endl;

	std::filesystem::path filepath(args.input_file);
	std::string
		file = filepath.string(),
		filename = filepath.stem().filename().string(),
		outOBJfile = UTILS::out_dir + remove_suffix(filename, "_sf.off__sf"),
		outMetafile = UTILS::out_dir + remove_suffix(filename, "_sf.off__sf");

	int side = args.inside ? -1 : 1;
	std::cout << file << std::endl;
	std::cout << "===================" << filename << "   " << args.dis * side << "===================" << std::endl;

	std::string
		tmpobj = outOBJfile + "_" + std::to_string(args.dis * side) + "_" + std::to_string(args._num_angle),
		tmpmeta = outMetafile + "_" + std::to_string(args.dis * side) + "_" + std::to_string(args._num_angle);

	std::vector<Offset3D::Source> sources;
	if (args.distance_file != "") {
		std::ifstream ifs(args.distance_file);
		int vid;
		double dis;
		while (ifs >> vid >> dis) {
			sources.push_back(Offset3D::Source(vid, dis));
		}
	}
	else if (!args.runGalleryDis) {
		std::cout << "No distance." << std::endl;
		return -1;
	}

	Offset3D::OffsetCrust crust(file, args._num_angle, args.blue_points, args.corner_ratio);
	crust.core_main(
		sources,
		args.inside,
		tmpobj, tmpmeta,
		args.relative_dis,
		args.outPoly,
		args.useCentroid,
		args.doSharp,
		args.runGalleryDis,
		args.diheral_plane
	);


	return 0;
}



