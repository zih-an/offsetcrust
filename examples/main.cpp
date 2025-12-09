#include <DupinCyclide.hpp>
#include <Knot.hpp>
#include <MAT_recon.hpp>
using namespace Offset3D;

int main(int argc, char** argv) {
	UTILS::sphere_file = R"(.\sphere.obj)";
	UTILS::out_dir = R"(.\)";

/*-----------Cyclide-----------*/
	//DupinCyclide cyclide;
	//cyclide.mc();
	//cyclide.para();
	//cyclide.smooth_main();
	//cyclide.mesh_main();
	//cyclide.mat();

	/*std::vector<int>
		usteps = { 4, 10, 40 },
		vsteps = { 4, 8, 16, 32 };
	for (int i = 0; i < usteps.size(); i++) {
		for (int j = 0; j < vsteps.size(); j++) {
			cyclide.mat(usteps[i], vsteps[j]);
		}
	}*/


/*-----------Knot-----------*/
	//for (int i = 0; i < 6; i++) {
	//	Knot knot(i);
	//	//knot.mat(300, 5);
	//	knot.mat(300, 150);
	//}



/*-----------MAT2Surf-----------*/
	//MATRecon mat1(UTILS::out_dir + "elephant_b.ma");
	MATRecon mat2(UTILS::out_dir + "bug.ma");
	//MATRecon mat3(UTILS::out_dir + "chair.ma");
	//MATRecon mat4(UTILS::out_dir + "elk.ma");
	//MATRecon mat5(UTILS::out_dir + "vase.ma");
	//MATRecon mat6(UTILS::out_dir + "plane.ma");


	return 0;
}
