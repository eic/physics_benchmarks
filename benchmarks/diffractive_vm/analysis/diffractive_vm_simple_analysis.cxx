#include "pleaseIncludeMe.h"
int diffractive_vm_simple_analysis(TString config_name)
{	

	// read our configuration	
	// std::ifstream  config_file{config_name};
	// nlohmann::json config;
	// config_file >> config;

	// const std::string rec_file      = config["rec_file"];
	// const std::string detector      = config["detector"];
	// const std::string output_prefix = config["output_prefix"];
	// const std::string test_tag      = config["test_tag"];

	TString name_for_input = "/gpfs02/eic/ztu/EPIC/physics/Simulation_Campaign_Oct2022/physics_benchmarks/local_data/tmp/18on110/rec-" + config_name + ".root";

	auto file=new TFile(name_for_input);
	auto tree = (TTree *) file->Get("events");
    TTreeReader tree_reader(tree);       // !the tree reader
    tree->Print();

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    // MC particle pz array for each MC particle
    TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

    TFile* output=new TFile("rec-"+config_name+"simple-output.root");
    TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

    	for(int itrk=0;itrk<reco_pz_array.size();itrk++){
    		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
    		h_eta->Fill(trk.Eta());
    	}
    }
	output->Write();
	output->Close();

	return 0;
}