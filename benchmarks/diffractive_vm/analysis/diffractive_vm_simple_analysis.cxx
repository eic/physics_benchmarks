#include "pleaseIncludeMe.h"
int diffractive_vm_simple_analysis(TString config_name)
{
	// read our configuration
	TString name_for_input = "/gpfs02/eic/ztu/EPIC/physics/Simulation_Campaign_Oct2022/physics_benchmarks/local_data/tmp/18on110/rec-" + config_name + ".root";

	auto file=new TFile(name_for_input);
	auto tree = (TTree *) file->Get("events");

	tree->Draw("MCParticles.momentum.z");


	return 0;
}