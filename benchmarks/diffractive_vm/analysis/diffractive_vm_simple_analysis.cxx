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

	TString name_of_input = (TString) rec_file;
	std::cout << "what is this rec_file = " << name_of_input << endl;
	
	auto file=new TFile(name_of_input);
	auto tree = (TTree *) file->Get("events");
    TTreeReader tree_reader(tree);       // !the tree reader
    
    // MC particle pz array for each MC particle
    TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<float> mc_mass_array = {tree_reader, "MCParticles.momentum.mass"};
    TTreeReaderArray<float> mc_charge_array = {tree_reader, "MCParticles.momentum.charge"};
    TTreeReaderArray<float> mc_pdg_array = {tree_reader, "MCParticles.momentum.PDG"};

    //Reconstructed EcalEndCapNClusters
    TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndCapNClusters.energy"};
    TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndCapNClusters.position.x"};
    TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndCapNClusters.position.y"};

    TTreeReaderArray<float> em_rec_id_array = {tree_reader, "EcalEndCapNClustersAssociations.recID"};
    TTreeReaderArray<float> em_sim_id_array = {tree_reader, "EcalEndCapNClustersAssociations.simID"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

    TString output_name_dir = output_prefix.c_str();
  	TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");

    TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
    TH1D* h_energy_MC = new TH1D("h_energy_MC",";E (GeV)",100,0,20);
    TH1D* h_energy_REC = new TH1D("h_energy_REC",";E (GeV)",100,0,20);
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-200,200,400,-200,200);

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

    	//MCParticles
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
    		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
    		if(mc_pdg_array[imc]==11 	
    			&& mctrk.Perp()>maxPt){
    			maxPt=mctrk.Perp();
    			mc_elect_index=imc;
    			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
    		}
    	}
    	h_energy_MC->Fill(scatMC.E());

    	//association to mc;
    	int rec_cluster_id=-1;
    	for(int iasso=0;iasso<em_rec_id_array.GetSize();iasso++){
    		if(em_sim_id_array[iasso]==mc_elect_index){
    			rec_cluster_id = em_rec_id_array[iasso];
    		}
    	}
    	
    	double energy = em_energy_array[rec_cluster_id];
    	double xpos = em_x_array[rec_cluster_id];
    	double ypos = em_y_array[rec_cluster_id];

    	h_energy_REC->Fill(energy);
    	h_emClus_position_REC->Fill(xpos,ypos);

    	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
    		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
    		h_eta->Fill(trk.Eta());
    	}
    }
	output->Write();
	output->Close();

	return 0;
}