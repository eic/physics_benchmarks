#include "pleaseIncludeMe.h"
int diffractive_vm_simple_analysis(const std::string& config_name)
{	

	// read our configuration	
	std::ifstream  config_file{config_name};
	nlohmann::json config;
	config_file >> config;

	const std::string rec_file      = config["rec_file"];
	const std::string detector      = config["detector"];
	const std::string output_prefix = config["output_prefix"];
	const std::string test_tag      = config["test_tag"];

	TString name_for_input = "/gpfs02/eic/ztu/EPIC/physics/Simulation_Campaign_Oct2022/physics_benchmarks/local_data/tmp/18on110/rec-" + config_name + ".root";

	auto file=new TFile(name_for_input);
	auto tree = (TTree *) file->Get("events");
    TTreeReader tree_reader(tree);       // !the tree reader
    tree->Print();

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    // MC particle pz array for each MC particle
    TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

        // Number of mc particles, reco particles and associations may differ
        fmt::print("New event. N reco particles: {}, N mc particles: {}, N assoc: {}\n",
                   reco_pz_array.GetSize(), mc_pz_array.GetSize(), rec_id.GetSize());

        // Iterate over associations
        for(unsigned int i=0; i < rec_id.GetSize(); i++) {

            // For each association pull index of reco and MC array
            auto reco_array_index = rec_id[i];
            auto mc_array_index = sim_id[i];

            float reco_pz = reco_pz_array[reco_array_index];
            float mc_pz = mc_pz_array[mc_array_index];
            fmt::print("   reco={:>10.4f} mc={:>10.4f}\n", reco_pz, mc_pz);
        }
    }
	

	return 0;
}