#include "pleaseIncludeMe.h"
int diffractive_vm_simple_analysis(const std::string& config_name, const int vm_type=1, const int mc_type=0)
{
	// read our configuration
	std::ifstream  config_file{config_name};
	nlohmann::json config;
	config_file >> config;

	const std::string rec_file      = config["rec_file"];
	const std::string detector      = config["detector"];
	const std::string output_prefix = config["output_prefix"];
	const std::string test_tag      = config["test_tag"];

	fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
	         "Running diffractive phi analysis...\n");
	fmt::print(" - Detector package: {}\n", detector);
	fmt::print(" - input file: {}\n", rec_file);
	fmt::print(" - output prefix: {}\n", output_prefix);
	fmt::print(" - test tag: {}\n", test_tag);

	TFile* file=new TFile(rec_file);
	TTree* tree = new TTree("events");

	tree->Draw("MCParticles.momentum.z");


	return 0;
}