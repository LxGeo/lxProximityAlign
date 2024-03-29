#include "parameters.h"
//#include "optim_proximity_aligner.h"
#include "optim_proximity_aligner_linear.h"
#include "optim_proximity_aligner.h"



using namespace LxGeo::lxProximityAlign;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	//GDALAllRegister();
	KGDAL2CV* kgdal2cv = new KGDAL2CV();

	// Reads command-line parameters

	params = new Parameters(argc, argv);
	if (!params->initialized()) {
		delete params;
		return 1;
	}

	// Runs process
	OptimProximityAligner d_t= OptimProximityAligner();
	if (d_t.pre_check())
		d_t.run();

	// Quits
	delete kgdal2cv;
	delete params;

	clock_t t_end = clock();
	//std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}