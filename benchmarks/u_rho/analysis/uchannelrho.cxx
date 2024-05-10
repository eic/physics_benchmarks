#include "pleaseIncludeMe.h"

auto giveme_t_method_L(TLorentzVector eIn, 
					   TLorentzVector eOut, 
					   TLorentzVector pIn, 
					   TLorentzVector vmOut)
{
	TLorentzVector aInVec(pIn.Px(),pIn.Py(),pIn.Pz(),sqrt(pIn.Px()*pIn.Px() + pIn.Py()*pIn.Py() + pIn.Pz()*pIn.Pz() + MASS_PROTON*MASS_PROTON) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_PROTON*MASS_PROTON + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();

	return method_L;
}

auto giveme_u(TLorentzVector pIn, TLorentzVector vmOut){
	double uvalue = -(vmOut-pIn).Mag2();
	return uvalue;
}


int uchannelrho(TString rec_file="input.root", TString outputfile="output.root")
{	
if (gSystem->AccessPathName(rec_file.Data()) != 0) {
   // File does not exist
   cout<<Form("File %s does not exist.", rec_file.Data())<<endl;
   return 0;
}

void setbenchstatus(double eff){
	///////////// Set benchmark status!
        // create our test definition
        // test_tag
        common_bench::Test rho_reco_eff_test{
          {
            {"name", "rho_reconstruction_efficiency"},
            {"title", "rho Reconstruction Efficiency for rho -> pi+pi- in the B0"},
            {"description", "u-channel rho->pi+pi- reconstruction efficiency "
                       "when both pions should be within B0 acceptance"},
            {"quantity", "efficiency"},
            {"target", "0.9"}
          }
        };      //these 2 need to be consistent 
        double eff_target = 0.9;    //going to find a way to use the same variable

        if(eff<0 || eff>1){
          rho_reco_eff_test.error(-1);
        }else if(eff > eff_target){
          rho_reco_eff_test.pass(eff);
        }else{
          rho_reco_eff_test.fail(eff);
        }


        // write out our test data
        common_bench::write_test(rho_reco_eff_test, "rhorecoeff.json");
}

// read our configuration	
auto tree = new TChain("events");
TString name_of_input = (TString) rec_file;
std::cout << "Input file = " << name_of_input << endl;
tree->Add(name_of_input);
TTreeReader tree_reader(tree);       // !the tree reader

TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
// MC particle pz array for each MC particle
TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
TTreeReaderArray<double> mc_endx_array = {tree_reader, "MCParticles.endpoint.x"};
TTreeReaderArray<double> mc_endy_array = {tree_reader, "MCParticles.endpoint.y"};
TTreeReaderArray<double> mc_endz_array = {tree_reader, "MCParticles.endpoint.z"};
TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

//Reconstructed EcalEndcapNClusters
TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};
TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};

TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClusterAssociations.recID"};
TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClusterAssociations.simID"};

//B0 Tracker hits
TTreeReaderArray<double> reco_x_array = {tree_reader, "MCParticles.endpoint.x"};//= {tree_reader, "B0TrackerHits.position.x"};
TTreeReaderArray<double> reco_y_array = {tree_reader, "MCParticles.endpoint.y"};//= {tree_reader, "B0TrackerHits.position.y"};
TTreeReaderArray<double> reco_z_array = {tree_reader, "MCParticles.endpoint.z"};//= {tree_reader, "B0TrackerHits.position.z"};

// Reconstructed particles pz array for each reconstructed particle
TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};

TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

TTreeReaderArray<int> reco_PDG = {tree_reader,"ReconstructedChargedParticles.PDG"};

TString output_name_dir = outputfile;
cout << "Output file = " << output_name_dir << endl;
TFile* output = new TFile(output_name_dir,"RECREATE");

//events
TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",100,0,20);
TH1D* h_y_e = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
TH1D* h_t_MC = new TH1D("h_t_MC",";t_{MC}; counts",1000,0,5);

TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,20);
TH1D* h_yREC_e = new TH1D("h_yREC_e",";y_{e,REC}",100,0,1);
TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);
TH1D* h_trk_energy_REC = new TH1D("h_trk_energy_REC",";E_{REC} (GeV)",100,0,20);
TH1D* h_trk_Epz_REC = new TH1D("h_trk_Epz_REC",";E - p_{z} (GeV)",200,0,50);

//track
TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
TH2D* h_trk_energy_res = new TH2D("h_trk_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} track-base ",100,0,20,1000,-1,1);
TH2D* h_trk_Pt_res = new TH2D("h_trk_Pt_res",";p_{T,MC} (GeV); P_{T,MC}-P_{T,REC}/P_{T,MC} track-base ",100,0,15,1000,-1,1);
TH1D* h_Epz_REC = new TH1D("h_Epz_REC",";E - p_{z} (GeV)",200,0,50);
TH2D* h_Acceptance_REC = new TH2D("h_Acceptance_REC",";#eta; p_{T}(GeV/c REC)",500,0,10,500,0,10);
TH2D* h_Acceptance_angular_REC = new TH2D("h_Acceptance_angular_REC",";#phi (radians); #theta (mrad)",100,0,6.4,200,0,100);
TH2D* h_Acceptance_angular_RECPI = new TH2D("h_Acceptance_angular_RECPI",";#phi (radians); #theta (mrad)",100,0,6.4,200,0,100);
TH2D* h_Acceptance_angular_RECPIP = new TH2D("h_Acceptance_angular_RECPIP",";#phi (radians); #theta (mrad)",100,0,6.4,200,0,100);
TH2D* h_Acceptance_angular_RECPIM = new TH2D("h_Acceptance_angular_RECPIM",";#phi (radians); #theta (mrad)",100,0,6.4,200,0,100);
TH2D* h_Acceptance_xy_RECPI = new TH2D("h_Acceptance_xy_RECPI",";x(mm); y(mm)",1000,-1000,1000,1000,-1000,1000);
TH2D* h_Acceptance_angular_MC  = new TH2D("h_Acceptance_angular_MC" ,";#phi (radians); #theta (mrad)",500,0,6.4,500,0,100);
TH2D* h_etavseta_MC = new TH2D("h_etavseta_MC",";#eta(#pi^{+} MC); #eta(#pi^{-} MC)",250,0,10,250,0,10);
TH2D* h_etavseta_REC = new TH2D("h_etavseta_REC",";#eta(#pi^{+} REC); #eta(#pi^{-} REC)",250,0,10,250,0,10);

TProfile2D* h_effEtaPtPi  = new TProfile2D("h_effEtaPtPi",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);
TProfile2D* h_effEtaPtPip = new TProfile2D("h_effEtaPtPip",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);
TProfile2D* h_effEtaPtPim = new TProfile2D("h_effEtaPtPim",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);

TProfile2D* h_effPhiEtaPi  = new TProfile2D("h_effPhiEtaPi",";#phi (rad);#eta",50,0,6.4,50,4,6);
TProfile2D* h_effPhiEtaPip = new TProfile2D("h_effPhiEtaPip",";#phi (rad);#eta",50,0,6.4,50,4,6);
TProfile2D* h_effPhiEtaPim = new TProfile2D("h_effPhiEtaPim",";#phi (rad);#eta",50,0,6.4,50,4,6);

TH2D* h_RecoMomPi  = new TH2D("h_RecoMomPi", ";p (GeV/c) MC;p (GeV/c) reco.",100,0,100,100,0,100);
TH2D* h_RecoMomPip = new TH2D("h_RecoMomPip",";p (GeV/c) MC;p (GeV/c) reco.",100,0,100,100,0,100);
TH2D* h_RecoMomPim = new TH2D("h_RecoMomPim",";p (GeV/c) MC;p (GeV/c) reco.",100,0,100,100,0,100);

TH2D* h_RecoTransMomPi  = new TH2D("h_RecoTransMomPi", ";p_{T} (GeV/c) MC;p_{T} (GeV/c) reco.",100,0,1.6,100,0,1.6);
TH2D* h_RecoTransMomPip = new TH2D("h_RecoTransMomPip",";p_{T} (GeV/c) MC;p_{T} (GeV/c) reco.",100,0,1.6,100,0,1.6);
TH2D* h_RecoTransMomPim = new TH2D("h_RecoTransMomPim",";p_{T} (GeV/c) MC;p_{T} (GeV/c) reco.",100,0,1.6,100,0,1.6);

TH1D* h_PiMomRecoQuality = new TH1D("h_PiMomRecoQuality",";(p_{reco}-p_{MC})/p_{MC};counts",500,-1,1);
TH1D* h_PipMomRecoQuality = new TH1D("h_PipMomRecoQuality",";(p_{reco}-p_{MC})/p_{MC};counts",500,-1,1);
TH1D* h_PimMomRecoQuality = new TH1D("h_PimMomRecoQuality",";(p_{reco}-p_{MC})/p_{MC};counts",500,-1,1);
TH1D* h_PiTransMomRecoQuality = new TH1D("h_PiTransMomRecoQuality",";(p_{T,reco}-p_{T,MC})/p_{T,MC};counts",500,-1,1);
TH1D* h_PipTransMomRecoQuality = new TH1D("h_PipTransMomRecoQuality",";(p_{T,reco}-p_{T,MC})/p_{T,MC};counts",500,-1,1);
TH1D* h_PimTransMomRecoQuality = new TH1D("h_PimTransMomRecoQuality",";(p_{T,reco}-p_{T,MC})/p_{T,MC};counts",500,-1,1);

//VM & t
TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_MC_etacut = new TH1D("h_VM_mass_MC_etacut",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_REC_etacut = new TH1D("h_VM_mass_REC_etacut",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_REC_justpions = new TH1D("h_VM_mass_REC_justpions",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_REC_justpionsB0 = new TH1D("h_VM_mass_REC_justpionsB0",";mass (GeV)",200,0,4);
TH1D* h_VM_mass_REC_notjustpionsB0 = new TH1D("h_VM_mass_REC_notjustpionsB0",";mass (GeV)",200,0,4);
TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T} (GeV/c)",200,0,2);
TH2D* h_VM_res = new TH2D("h_VM_res",";p_{T,MC} (GeV); p_{T,MC}-E_{T,REC}/p_{T,MC}",100,0,2,1000,-1,1);
TH1D* h_t_REC = new TH1D("h_t_REC",";t_{REC} (GeV^{2}); counts",1000,0,5);
TH1D* h_t_trk_REC = new TH1D("h_t_trk_REC",";t_{REC}(GeV^{2}) track-base; counts",1000,0,5);
TH1D* h_t_combo_REC = new TH1D("h_t_combo_REC",";t_{combo,REC}(GeV^{2}); counts",1000,0,5);
TH2D* h_t_res = new TH2D("h_t_res",";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC}",1000,0,5,1000,-10,10);
TH2D* h_trk_t_res = new TH2D("h_trk_t_res",";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC} track-base",1000,0,5,1000,-10,10);
TH2D* h_t_2D = new TH2D("h_t_2D",";t_{MC} (GeV^{2}); t_{REC} (GeV^{2}) track-base",1000,0,5,1000,0,5);
TH2D* h_t_REC_2D = new TH2D("h_t_REC_2D",";t_{trk,REC} (GeV^{2}); t_{EEMC,REC} (GeV^{2})",1000,0,5,1000,0,5);
TH2D* h_t_RECMC_2D = new TH2D("h_t_RECMC_2D",";t_{MC} (GeV^{2}); t_{trk,REC} / t_{EEMC,REC} ",1000,0,5,200,-10,10);
TH2D* h_VM_endpointXY_MC = new TH2D("h_VM_endpointXY_MC",";x (cm); y (cm)",1000,-300,300,1000,-300,300);
TH1D* h_VM_endpointZ_MC = new TH1D("h_VM_endpointZ_MC",";z (cm); counts",1000,0,1000);
TH2D* h_u_REC_2D = new TH2D("h_u_REC_2D",";-#it{u}_{MC} (GeV^{2}); -#it{u}_{REC} (GeV^{2})",1000,0,5,1000,0,5);
double umin = -1.0;
double umax = 5.0;
int n_ubins = 100;
double u_binwidth = (umax-umin)/((double)n_ubins);
TH1D* h_u_REC = new TH1D("h_u_REC", ";-#it{u}_{REC} (GeV^{2}); dN/d|#it{u}| (GeV^{-2} scaled)",n_ubins,umin,umax);
TH1D* h_u_REC_justpions = new TH1D("h_u_REC_justpions", ";-#it{u}_{REC} (GeV^{2}); dN/d|#it{u}| (GeV^{-2} scaled)",n_ubins,umin,umax);
TH1D* h_u_REC_justpionsB0 = new TH1D("h_u_REC_justpionsB0", ";-#it{u}_{REC} (GeV^{2}); dN/d|#it{u}| (GeV^{-2} scaled)",n_ubins,umin,umax);
TH1D* h_u_REC_notjustpionsB0 = new TH1D("h_u_REC_notjustpionsB0", ";-#it{u}_{REC} (GeV^{2}); dN/d|#it{u}| (GeV^{-2} scaled)",n_ubins,umin,umax);
TH1D* h_u_MC = new TH1D("h_u_MC", ";-#it{u}_{MC} (GeV^{2}); dN/d|#it{u}| (GeV^{-2} scaled)",n_ubins,umin,umax);

//energy clus
TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
TH2D* h_emHits_position_REC = new TH2D("h_emHits_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} emcal",100,0,20,1000,-1,1);
TH1D* h_energy_calibration_REC = new TH1D("h_energy_calibration_REC",";E (GeV)",200,0,2);
TH1D* h_EoverP_REC = new TH1D("h_EoverP_REC",";E/p",200,0,2);
TH1D* h_ClusOverHit_REC = new TH1D("h_ClusOverHit_REC",";cluster energy / new cluster energy",200,0,2);

//
TH1D* h_numPositiveTracks = new TH1D("h_numPositiveTracks",";number of positive tracks;counts",10,-0.5,9.5);

cout<<"There are "<<tree->GetEntries()<<" events!!!!!!!"<<endl;
tree_reader.SetEntriesRange(0, tree->GetEntries());
while (tree_reader.Next()) {

	/*
	Beam particles
	*/
	TLorentzVector ebeam(0,0,0,0);
	TLorentzVector pbeam(0,0,0,0);

	TLorentzVector vmMC(0,0,0,0);
	TLorentzVector piplusMC(0,0,0,0);
	TLorentzVector piminusMC(0,0,0,0);

	//MC level
	TLorentzVector scatMC(0,0,0,0);
	unsigned int mc_elect_index=-1;
	double maxPt=-99.;
	for(unsigned int imc=0;imc<mc_px_array.GetSize();imc++){
		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
		if(mc_pdg_array[imc]!=11) mctrk.RotateY(0.025);//Rotate all non-electrons to hadron beam coordinate system
		if(mc_genStatus_array[imc]==4){//4 is Sartre.
			if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
		}
		if(mc_genStatus_array[imc]!=1) continue;
		if(mc_pdg_array[imc]==11 	
			&& mctrk.Perp()>maxPt){
			maxPt=mctrk.Perp();
			mc_elect_index=imc;
			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
		}
		if(mc_pdg_array[imc]==211 && mc_genStatus_array[imc]==1){ 
		  piplusMC.SetVectM(mctrk,MASS_PION);  
		  h_VM_endpointXY_MC->Fill(mc_endx_array[imc],mc_endy_array[imc]);  
		  h_VM_endpointZ_MC->Fill(mc_endz_array[imc]);
		  double phi = mctrk.Phi()>0 ? mctrk.Phi() : mctrk.Phi()+6.2831853;
                  h_Acceptance_angular_MC->Fill(phi,1000.0*mctrk.Theta());
		}
                if(mc_pdg_array[imc]==-211 && mc_genStatus_array[imc]==1){ 
		  piminusMC.SetVectM(mctrk,MASS_PION); 
		  h_VM_endpointXY_MC->Fill(mc_endx_array[imc],mc_endy_array[imc]);  
		  h_VM_endpointZ_MC->Fill(mc_endz_array[imc]);
                  double phi = mctrk.Phi()>0 ? mctrk.Phi() : mctrk.Phi()+6.2831853;
                  h_Acceptance_angular_MC->Fill(phi,1000.0*mctrk.Theta());
		}
	}
	vmMC=piplusMC+piminusMC;
	h_etavseta_MC->Fill(piplusMC.Eta(),piminusMC.Eta());
	//protection.
	if(ebeam.E()==pbeam.E() && ebeam.E()==0) {
		std::cout << "problem with MC incoming beams" << std::endl;
		continue;
	}
	TLorentzVector qbeam=ebeam-scatMC;
	double Q2=-(qbeam).Mag2();  
	double pq=pbeam.Dot(qbeam);
	double y= pq/pbeam.Dot(ebeam);
	
	//MC level phase space cut
	//if(Q2<1.||Q2>10.) continue;
	//if(y<0.01||y>0.85) continue;

	h_Q2_e->Fill(Q2);
	h_y_e->Fill(y);
	h_energy_MC->Fill(scatMC.E());

	double t_MC=0.;
	if(vmMC.E()!=0)
	{
		double method_E = -(qbeam-vmMC).Mag2();
		t_MC=method_E;
		h_t_MC->Fill( method_E );
		h_VM_mass_MC->Fill(vmMC.M());
		if(piplusMC.Theta()>0.009  && piplusMC.Theta()<0.013 && 
		    piminusMC.Theta()>0.009 && piminusMC.Theta()<0.013 ) h_VM_mass_MC_etacut->Fill(vmMC.M());
	}

	//rec level
	//leading cluster
	double maxEnergy=-99.;
	double xpos=-999.;
	double ypos=-999.;
	for(unsigned int iclus=0;iclus<em_energy_array.GetSize();iclus++){
		if(em_energy_array[iclus]>maxEnergy){
			maxEnergy=em_energy_array[iclus];
			xpos=em_x_array[iclus];
			ypos=em_y_array[iclus];
		}
	}
	//leading hit energy
	double maxHitEnergy=0.01;//threshold 10 MeV
	double xhitpos=-999.;
	double yhitpos=-999.;
	int hit_index=-1;
	for(unsigned int ihit=0;ihit<emhits_energy_array.GetSize();ihit++){	
		if(emhits_energy_array[ihit]>maxHitEnergy){
			maxHitEnergy=emhits_energy_array[ihit];
			        xhitpos=emhits_x_array[ihit];
			        yhitpos=emhits_y_array[ihit];
			        hit_index=ihit;
		}
	}
	//sum over all 3x3 towers around the leading tower
	double xClus=xhitpos*maxHitEnergy;
	double yClus=yhitpos*maxHitEnergy;
	for(unsigned int ihit=0;ihit<emhits_energy_array.GetSize();ihit++){
		double hitenergy=emhits_energy_array[ihit];
		double x=emhits_x_array[ihit];
		double y=emhits_y_array[ihit];
		double d=sqrt( (x-xhitpos)*(x-xhitpos) + (y-yhitpos)*(y-yhitpos));
		if(d<70. && ihit!=hit_index && hitenergy>0.01)  {
			maxHitEnergy+=hitenergy;//clustering around leading tower 3 crystal = 60mm.
			xClus+=x*hitenergy;
			yClus+=y*hitenergy;
		}
	}
	
	h_ClusOverHit_REC->Fill( maxEnergy / maxHitEnergy );
	//weighted average cluster position.
	xClus = xClus/maxHitEnergy;
	yClus = yClus/maxHitEnergy;
	double radius=sqrt(xClus*xClus+yClus*yClus);
	//if( radius>550. ) continue; //geometric acceptance cut
	//4.4% energy calibration.
	double clusEnergy=1.044*maxHitEnergy; 

	h_energy_REC->Fill(clusEnergy);
	//ratio of reco / truth Energy
	h_energy_calibration_REC->Fill( clusEnergy / scatMC.E() );
	//energy resolution
	double res= (scatMC.E()-clusEnergy)/scatMC.E();
	h_energy_res->Fill(scatMC.E(), res);
	h_emClus_position_REC->Fill(xpos,ypos);//default clustering position
	h_emHits_position_REC->Fill(xClus,yClus); //self clustering position
	
	//association of rec level scat' e
	int rec_elect_index=-1;
	for(unsigned int i=0;i<sim_id.GetSize();i++){
		if(sim_id[i]==mc_elect_index){
			//find the rec_id
			rec_elect_index = rec_id[i];
		}
	}
    
    TLorentzVector scatMCmatchREC(0,0,0,0);
    TLorentzVector scatREC(0,0,0,0);
    TLorentzVector scatClusEREC(0,0,0,0);
    TLorentzVector hfs(0,0,0,0);
    TLorentzVector particle(0,0,0,0);
    TLorentzVector protonREC(0,0,0,0);
    TLorentzVector protonRECasifpion(0,0,0,0);
    TLorentzVector piplusREC(0,0,0,0);
    TLorentzVector piminusREC(0,0,0,0);
    TLorentzVector vmREC(0,0,0,0);
    TLorentzVector vmRECfromproton(0,0,0,0);

    bool isPiMinusFound = false;
    bool isPiPlusFound = false;
    bool isProtonFound = false;

    double maxP=-1.;
    //track loop
	for(unsigned int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
		if(rec_elect_index!=-1
			&& itrk==rec_elect_index){
			scatMCmatchREC.SetVectM(trk,MASS_ELECTRON);//Reserved to calculate t.
		}
		if(trk.Mag()>maxP){
			//track-base 4 vector
			maxP=trk.Mag();
			scatREC.SetVectM(trk,MASS_ELECTRON);

			//use emcal energy to define 4 vector
			double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
			double eta=scatREC.Eta();
			double phi=scatREC.Phi();
			double pt = TMath::Sin(scatREC.Theta())*p;
			scatClusEREC.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
		}
	}
	//loop over track again;
	int numpositivetracks = 0;
	for(unsigned int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
		trk.RotateY(0.025);	
		particle.SetVectM(trk,MASS_PION);//assume pions;
		if(itrk!=rec_elect_index) {
    		hfs += particle; //hfs 4vector sum.
    		//selecting rho->pi+pi- daughters;
    		h_eta->Fill(trk.Eta());
    		//if(fabs(trk.Eta())<3.0){
    			if(reco_charge_array[itrk]>0){ 
			  numpositivetracks++; 
			  if(reco_PDG[itrk]==211){
			    piplusREC.SetVectM(trk,MASS_PION); 
			    isPiPlusFound=true;
			  }
                          if(reco_PDG[itrk]==2212){
                            protonREC.SetVectM(trk,MASS_PROTON); 
                            protonRECasifpion.SetVectM(trk,MASS_PION);
                            isProtonFound=true; 
                          }
			}
    			if(reco_charge_array[itrk]<0){ piminusREC.SetVectM(trk,MASS_PION); if(reco_PDG[itrk]==-211)isPiMinusFound=true;}
    		//}
		double pt = sqrt(reco_px_array[itrk]*reco_px_array[itrk] + reco_py_array[itrk]*reco_py_array[itrk]);
		h_Acceptance_REC->Fill(fabs(trk.Eta()),pt);
	        double phi = trk.Phi()>0 ? trk.Phi() : trk.Phi()+6.2831853;
		h_Acceptance_angular_REC->Fill(phi,1000.0*trk.Theta());
		if(reco_PDG[itrk]==211 || reco_PDG[itrk]==-211) h_Acceptance_angular_RECPI->Fill(phi,1000.0*trk.Theta());
                if(reco_PDG[itrk]==211) h_Acceptance_angular_RECPIP->Fill(phi,1000.0*trk.Theta());
                if(reco_PDG[itrk]==-211) h_Acceptance_angular_RECPIM->Fill(phi,1000.0*trk.Theta());
		}
	}
	h_numPositiveTracks->Fill(numpositivetracks);
        //loop over B0 hits
	for(unsigned int ihit=0; ihit<reco_x_array.GetSize();ihit++){
		TVector3 hit(reco_x_array[ihit],reco_y_array[ihit],reco_z_array[ihit]);
		hit.RotateY(0.025);
		h_Acceptance_xy_RECPI->Fill(hit.X(),hit.Y());
	}
	//4vector of VM;
	if(piplusREC.E()!=0. && piminusREC.E()!=0.){
		vmREC=piplusREC+piminusREC;
	}
        if(protonRECasifpion.E()!=0. && piminusREC.E()!=0.){
                vmRECfromproton=protonRECasifpion+piminusREC;
        }

	//pion reconstruction
	double thispipeff  = (isPiPlusFound) ? 1.0 : 0.0;
        double thispimeff = (isPiMinusFound) ? 1.0 : 0.0;
	h_effEtaPtPi ->Fill(piplusMC.Eta(), piplusMC.Pt(), thispipeff);
        h_effEtaPtPi ->Fill(piminusMC.Eta(),piminusMC.Pt(),thispimeff);
        h_effEtaPtPip->Fill(piplusMC.Eta(), piplusMC.Pt(), thispipeff);
        h_effEtaPtPim->Fill(piminusMC.Eta(),piminusMC.Pt(),thispimeff);
	//
	double thispipphi = piplusMC.Phi()>0  ? piplusMC.Phi()  : piplusMC.Phi()+6.2831853;
        double thispimphi = piminusMC.Phi()>0 ? piminusMC.Phi() : piminusMC.Phi()+6.2831853;
        if(piplusMC.Pt() >0.2) h_effPhiEtaPi ->Fill(thispipphi, piplusMC.Eta(), thispipeff);
        if(piminusMC.Pt()>0.2) h_effPhiEtaPi ->Fill(thispimphi,piminusMC.Eta(),thispimeff);
        if(piplusMC.Pt() >0.2) h_effPhiEtaPip->Fill(thispipphi, piplusMC.Eta(), thispipeff);
        if(piminusMC.Pt()>0.2) h_effPhiEtaPim->Fill(thispimphi,piminusMC.Eta(),thispimeff);
	//
	if(isPiPlusFound)  h_RecoMomPi ->Fill(piplusMC.P() ,piplusREC.P() );
        if(isPiMinusFound) h_RecoMomPi ->Fill(piminusMC.P(),piminusREC.P());
        if(isPiPlusFound)  h_RecoMomPip->Fill(piplusMC.P() ,piplusREC.P() );
        if(isPiMinusFound) h_RecoMomPim->Fill(piminusMC.P(),piminusREC.P());
        //
        if(isPiPlusFound)  h_RecoTransMomPi ->Fill(piplusMC.Pt() ,piplusREC.Pt() );
        if(isPiMinusFound) h_RecoTransMomPi ->Fill(piminusMC.Pt(),piminusREC.Pt());
        if(isPiPlusFound)  h_RecoTransMomPip->Fill(piplusMC.Pt() ,piplusREC.Pt() );
        if(isPiMinusFound) h_RecoTransMomPim->Fill(piminusMC.Pt(),piminusREC.Pt());
	//
	if(isPiPlusFound){
	  h_PiMomRecoQuality->Fill((piplusREC.P()-piplusMC.P())/piplusMC.P());
          h_PipMomRecoQuality->Fill((piplusREC.P()-piplusMC.P())/piplusMC.P());
          h_PiTransMomRecoQuality->Fill((piplusREC.Pt()-piplusMC.Pt())/piplusMC.Pt());
          h_PipTransMomRecoQuality->Fill((piplusREC.Pt()-piplusMC.Pt())/piplusMC.Pt());
	}
        if(isPiMinusFound){
          h_PiMomRecoQuality->Fill((piminusREC.P()-piminusMC.P())/piminusMC.P());
          h_PimMomRecoQuality->Fill((piminusREC.P()-piminusMC.P())/piminusMC.P());
          h_PiTransMomRecoQuality->Fill((piminusREC.Pt()-piminusMC.Pt())/piminusMC.Pt());
          h_PimTransMomRecoQuality->Fill((piminusREC.Pt()-piminusMC.Pt())/piminusMC.Pt());
        }

	//track-base e' energy REC;
	h_trk_energy_REC->Fill(scatMCmatchREC.E());
	
	//track-base e' energy resolution;
	res= (scatMC.E()-scatMCmatchREC.E())/scatMC.E();
	h_trk_energy_res->Fill(scatMC.E(), res);
	
	//track-base e' pt resolution;
	res= (scatMC.Pt()-scatMCmatchREC.Pt())/scatMC.Pt();
	h_trk_Pt_res->Fill(scatMC.Pt(), res);

	//track-base Epz scat' e
	double EpzREC= (scatMCmatchREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
	h_trk_Epz_REC->Fill( EpzREC );

	//EEMC cluster Epz scat' e
	EpzREC= (scatClusEREC+hfs).E() - (scatClusEREC+hfs).Pz();
	h_Epz_REC->Fill( EpzREC );

	//E over p
	double EoverP=scatClusEREC.E() / scatMCmatchREC.P();
	h_EoverP_REC->Fill( EoverP );

	//cluster-base DIS kine;
	TLorentzVector qbeamREC=ebeam-scatClusEREC;
	double Q2REC=-(qbeamREC).Mag2();  
	double pqREC=pbeam.Dot(qbeamREC);
	double yREC= pqREC/pbeam.Dot(ebeam);
	h_Q2REC_e->Fill(Q2REC);
	h_yREC_e->Fill(yREC);

	//Event selection:
	//if( EpzREC<27||EpzREC>40 ) continue;
	//if( EoverP<0.8||EoverP>1.18 ) continue;		

	//MC level phase space cut
	//if(Q2REC<1.||Q2REC>10.) continue;
	//if(yREC<0.01||yREC>0.85) continue;

	//VM rec
	if(vmREC.E()==0 && vmRECfromproton.E()==0) continue;
	double rho_mass = vmREC.M();
        double rho_mass_fromproton = vmRECfromproton.M();
	h_VM_mass_REC->Fill(rho_mass);
        h_VM_mass_REC->Fill(rho_mass_fromproton);
	if(piplusMC.Theta()>0.009  && piplusMC.Theta()<0.013 &&
                    piminusMC.Theta()>0.009 && piminusMC.Theta()<0.013 ) h_VM_mass_REC_etacut->Fill(vmREC.M());
        if(piplusREC.Eta()>3.9 && piminusREC.Eta()>3.9) h_VM_mass_REC_notjustpionsB0->Fill(rho_mass);
	if(isPiMinusFound && isPiPlusFound){ 
	  h_VM_mass_REC_justpions->Fill(rho_mass);
	  h_etavseta_REC->Fill(piplusREC.Eta(),piminusREC.Eta());
	  if(piplusREC.Eta()>3.9 && piminusREC.Eta()>3.9) h_VM_mass_REC_justpionsB0->Fill(rho_mass);
	}
	h_VM_pt_REC->Fill(vmREC.Pt());

    	//2 versions: track and energy cluster:
	double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,pbeam,vmREC);
    	double t_REC = giveme_t_method_L(ebeam,scatClusEREC,pbeam,vmREC);
    	h_t_trk_REC->Fill( t_trk_REC );
    	h_t_REC->Fill( t_REC );
    	h_t_REC_2D->Fill(t_trk_REC,t_REC);
	double u_REC = giveme_u(pbeam,vmREC);
        double u_REC_fromproton = giveme_u(pbeam,vmRECfromproton);
	double u_MC = giveme_u(pbeam,vmMC);
        //4vector of VM;
        if(piplusREC.E()!=0. && piminusREC.E()!=0.){
                h_u_REC->Fill(u_REC);
                h_u_REC_justpions->Fill(u_REC);
        }
        if(protonRECasifpion.E()!=0. && piminusREC.E()!=0.){
                h_u_REC->Fill(u_REC_fromproton);
        }
	h_u_REC_2D->Fill(u_MC,u_REC);
	/*h_u_REC->Fill(u_REC);
	if(piplusREC.Eta()>3.9 && piminusREC.Eta()>3.9) h_u_REC_notjustpionsB0->Fill(u_REC);
        if(isPiMinusFound && isPiPlusFound){ 
          h_u_REC_justpions->Fill(u_REC);
          if(piplusREC.Eta()>3.9 && piminusREC.Eta()>3.9) h_u_REC_justpionsB0->Fill(u_REC);
        }
	*/
	h_u_MC->Fill(u_MC);
    	if( (t_trk_REC/t_REC) > 0.5 && (t_trk_REC/t_REC) < 1.5 ){
    		h_t_combo_REC->Fill( (t_trk_REC+t_REC)/2. );//w=1./(fabs(1.0-(t_trk_REC/t_REC)))
    	}
    	h_t_RECMC_2D->Fill(t_MC,t_trk_REC/t_REC);

    	//t track resolution 
	res= (t_MC-t_trk_REC)/t_MC;
	h_trk_t_res->Fill(t_MC, res);
	
    	//t EEMC resolution;
	res= (t_MC-t_REC)/t_MC;
	h_t_res->Fill(t_MC, res);
	
	//2D t
	h_t_2D->Fill(t_MC,t_trk_REC);

    	//VM pt resolution;
    	res= (vmMC.Pt()-vmREC.Pt())/vmMC.Pt();
	h_VM_res->Fill(vmMC.Pt(), res);

}
h_u_REC->Scale(1.0/u_binwidth);
h_u_REC_justpions->Scale(1.0/u_binwidth);
h_u_REC_justpionsB0->Scale(1.0/u_binwidth);
h_u_REC_notjustpionsB0->Scale(1.0/u_binwidth);
h_u_MC->Scale(1.0/u_binwidth);

double minbineff = h_VM_mass_MC_etacut->FindBin(0.6);
double maxbineff = h_VM_mass_MC_etacut->FindBin(1.0);
double rhorecoeff = (1.0*h_VM_mass_REC_etacut->Integral(minbineff,maxbineff))/(1.0*h_VM_mass_MC_etacut->Integral(minbineff,maxbineff));

output->Write();
output->Close();


setbenchstatus(rhorecoeff);

return 0;
}
