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
    
    TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};

    // MC particle pz array for each MC particle
    TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

    //Reconstructed EcalEndcapNClusters
    TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
    TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
    TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};

    TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClustersAssociations.recID"};
    TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClustersAssociations.simID"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

    TString output_name_dir = output_prefix.c_str();
  	TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");

  	//events
    TH1D* h_Q2_e = new TH1D("h_Q2_e",";#eta",100,0,20);
    TH1D* h_y_e = new TH1D("h_y_e",";#eta",100,0,1);
 	TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
    TH1D* h_t_MC = new TH1D("h_t_MC",";t; counts",100,0,0.2);

    TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";#eta",100,0,20);
    TH1D* h_yREC_e = new TH1D("h_yREC_e",";#eta",100,0,1);
    TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);
    TH1D* h_trk_energy_REC = new TH1D("h_trk_energy_REC",";E_{REC} (GeV)",100,0,20);
    
    //track
    TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
    TH2D* h_trk_energy_res = new TH2D("h_trk_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} track ",100,0,20,1000,-1,1);

   	//energy clus
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
    TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} emcal",100,0,20,1000,-1,1);

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

		/*
	    Beam particles
	    */
    	TLorentzVector ebeam(0,0,0,0);
    	TLorentzVector pbeam(0,0,0,0);

    	TLorentzVector vmMC(0,0,0,0);
    	TLorentzVector kplusMC(0,0,0,0);
    	TLorentzVector kminusMC(0,0,0,0);

    	//MC level
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
    		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
    		if(mc_genStatus_array[imc]==4){
    			if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
  				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
    		}
    		if(mc_pdg_array[imc]==11 	
    			&& mctrk.Perp()>maxPt){
    			maxPt=mctrk.Perp();
    			mc_elect_index=imc;
    			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
    		}
    		if(mc_pdg_array[imc]==321
    			&& mc_genStatus_array[imc]==1) kplusMC.SetVectM(mctrk,MASS_KAON);
			if(mc_pdg_array[imc]==-321
    			&& mc_genStatus_array[imc]==1) kminusMC.SetVectM(mctrk,MASS_KAON);

    	}
    	vmMC=kplusMC+kminusMC;
    	//protection.
    	if(ebeam.E()==pbeam.E() && ebeam.E()==0) {
    		std::cout << "problem with MC incoming beams" << std::endl;
    		continue;
    	}
    	TLorentzVector qbeam=ebeam-scatMC;
    	double Q2=-(qbeam).Mag2();  
		double pq=pbeam.Dot(qbeam);
		double y= pq/pbeam.Dot(ebeam);
		if(Q2<1.||Q2>10.) continue;
		if(y<0.01||y>0.95) continue;

		h_Q2_e->Fill(Q2);
		h_y_e->Fill(y);
    	h_energy_MC->Fill(scatMC.E());

    	if(vmMC.E()!=0){
    		double method_E = -(qbeam-vmMC).Mag2();
    		h_t_MC->Fill( method_E );
    	}


    	//rec level
    	double maxEnergy=-99.;
    	double xpos=-999.;
    	double ypos=-999.;
    	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++){
    		if(em_energy_array[iclus]>maxEnergy){
    			maxEnergy=em_energy_array[iclus];
    			xpos=em_x_array[iclus];
    			ypos=em_y_array[iclus];
    		}
    	}

    	double res= (scatMC.E()-maxEnergy)/scatMC.E();
		h_energy_res->Fill(scatMC.E(), res);
    	
		h_energy_REC->Fill(maxEnergy);
		h_emClus_position_REC->Fill(xpos,ypos);

		//association of rec level scat' e
		int rec_elecct_index=-1;
		for(int i=0;i<sim_id.GetSize();i++){
			if(sim_id[i]==mc_elect_index){
				//find the rec_id
				rec_elecct_index = rec_id[i];
			}
		}
	    
	    TLorentzVector scatMCmatchREC(0,0,0,0);
	    TLorentzVector scatREC(0,0,0,0);
	    double maxP=-1.;
	    //track loop
    	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
    		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
    		if(rec_elecct_index!=-1
    			&& itrk==rec_elecct_index){
    			scatMCmatchREC.SetVectM(trk,MASS_ELECTRON);
    		}
    		if(trk.Mag()>maxP){
    			maxP=trk.Mag();
    			scatREC.SetVectM(trk,MASS_ELECTRON);
    		}

    		h_eta->Fill(trk.Eta());
    	}
    	if(scatREC.E()==0) continue;

		//track-base DIS kine;
		TLorentzVector qbeamREC=ebeam-scatREC;
    	double Q2REC=-(qbeamREC).Mag2();  
		double pqREC=pbeam.Dot(qbeamREC);
		double yREC= pqREC/pbeam.Dot(ebeam);
		h_Q2REC_e->Fill(Q2REC);
		h_yREC_e->Fill(yREC);
		h_trk_energy_REC->Fill(scatREC.E());

		res= (scatMC.E()-scatREC.E())/scatMC.E();
		h_trk_energy_res->Fill(scatMC.E(), res);
    }
	output->Write();
	output->Close();

	return 0;
}