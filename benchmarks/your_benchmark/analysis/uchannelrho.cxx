#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector2.h"
#include "TVector3.h"

#include "fmt/color.h"
#include "fmt/core.h"

#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

int uchannelrho(TString rec_file="input.root", TString outputfile="output.root"){	
	if (gSystem->AccessPathName(rec_file.Data()) != 0) {
	   // File does not exist
	   cout<<Form("File %s does not exist.", rec_file.Data())<<endl;
	   return 0;
	}
	
	// read our configuration	
	auto tree = new TChain("events");
	TString name_of_input = (TString) rec_file;
	std::cout << "Input file = " << name_of_input << endl;
	tree->Add(name_of_input);
	TTreeReader tree_reader(tree);       // !the tree reader
	
	//MC-level track attributes
	TTreeReaderArray<int>   mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
	TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
	TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
	TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
	TTreeReaderArray<int>   mc_pdg_array= {tree_reader, "MCParticles.PDG"};

	//Reco-level track attributes
	TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
	TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
	TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
	TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};
	TTreeReaderArray<int>   reco_type     = {tree_reader,"ReconstructedChargedParticles.type"};
	
	TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
	TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};
	
	TString output_name_dir = outputfile;
	cout << "Output file = " << output_name_dir << endl;
	TFile* output = new TFile(output_name_dir,"RECREATE");
	
	//Pion reconstruction efficiency histograms
	TProfile2D* h_effEtaPtPi  = new TProfile2D("h_effEtaPtPi",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);
	TProfile2D* h_effEtaPtPip = new TProfile2D("h_effEtaPtPip",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);
	TProfile2D* h_effEtaPtPim = new TProfile2D("h_effEtaPtPim",";#eta; p_{T}(GeV/c)",50,4,6,30,0,1.6);
	TProfile2D* h_effPhiEtaPi  = new TProfile2D("h_effPhiEtaPi",";#phi (rad);#eta",50,0,6.4,50,4,6);
	TProfile2D* h_effPhiEtaPip = new TProfile2D("h_effPhiEtaPip",";#phi (rad);#eta",50,0,6.4,50,4,6);
	TProfile2D* h_effPhiEtaPim = new TProfile2D("h_effPhiEtaPim",";#phi (rad);#eta",50,0,6.4,50,4,6);
	
	//rho vector meson mass
	TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";mass (GeV)",200,0,4);
	TH1D* h_VM_mass_MC_etacut = new TH1D("h_VM_mass_MC_etacut",";mass (GeV)",200,0,4);
	TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";mass (GeV)",200,0,4);
	TH1D* h_VM_mass_REC_etacut = new TH1D("h_VM_mass_REC_etacut",";mass (GeV)",200,0,4);
	TH1D* h_VM_mass_REC_justpions = new TH1D("h_VM_mass_REC_justpions",";mass (GeV)",200,0,4);
	
	cout<<"There are "<<tree->GetEntries()<<" events!!!!!!!"<<endl;
	tree_reader.SetEntriesRange(0, tree->GetEntries());
	while (tree_reader.Next()) {
	
		TLorentzVector vmMC(0,0,0,0);
		TLorentzVector piplusMC(0,0,0,0);
		TLorentzVector piminusMC(0,0,0,0);
		
		//MC level
		for(unsigned int imc=0;imc<mc_px_array.GetSize();imc++){
			TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
			if(mc_pdg_array[imc]!=11) mctrk.RotateY(0.025);//Rotate all non-electrons to hadron beam coordinate system
			if(mc_genStatus_array[imc]!=1) continue;
			if(mc_pdg_array[imc]==211 && mc_genStatus_array[imc]==1){ 
			  piplusMC.SetVectM(mctrk,MASS_PION);  
			}
		  if(mc_pdg_array[imc]==-211 && mc_genStatus_array[imc]==1){ 
				piminusMC.SetVectM(mctrk,MASS_PION); 
			}
		}
		
		vmMC=piplusMC+piminusMC;
		if(vmMC.E()!=0){
			h_VM_mass_MC->Fill(vmMC.M());
			if(piplusMC.Theta()>0.009  && piplusMC.Theta()<0.013 && 
			    piminusMC.Theta()>0.009 && piminusMC.Theta()<0.013 ) h_VM_mass_MC_etacut->Fill(vmMC.M());
		}
		
		//4-vector for proton reconstructed as if it were a pi+
		TLorentzVector protonRECasifpion(0,0,0,0);
		//pion 4-vectors
		TLorentzVector piplusREC(0,0,0,0);
		TLorentzVector piminusREC(0,0,0,0);
		//rho reconstruction 4-vector
		TLorentzVector vmREC(0,0,0,0);
		//rho reconstruction from mis-identified proton 4-vector
		TLorentzVector vmRECfromproton(0,0,0,0);
		
		bool isPiMinusFound = false;
		bool isPiPlusFound = false;
		bool isProtonFound = false;
		
		//track loop
		int numpositivetracks = 0;
		int failed = 0;
		for(unsigned int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
			TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
			
			//  Rotate in order to account for crossing angle 
			//  and express coordinates in hadron beam pipe frame
			//  This is just a patch, not a final solution.
			trk.RotateY(0.025);
	
			if(reco_type[itrk] == -1){ 
				failed++;
				continue;
			}
	

		  if(reco_charge_array[itrk]>0){ 
				numpositivetracks++; 
			  if ((sim_id[itrk - failed]==4 || sim_id[itrk - failed]==5) && reco_charge_array[itrk - failed]==1){
			    piplusREC.SetVectM(trk,MASS_PION); 
			    isPiPlusFound=true;
			  }
		     if(sim_id[itrk - failed]==6){
		     	protonRECasifpion.SetVectM(trk,MASS_PION);
		     	isProtonFound=true; 
		     }
			}
		  if(reco_charge_array[itrk]<0){ 
		  	piminusREC.SetVectM(trk,MASS_PION); 
		  	if((sim_id[itrk - failed]==4 || sim_id[itrk - failed]==5) && reco_charge_array[itrk - failed]==-1)	isPiMinusFound=true;
		  }
			
		}
		
		//4vector of VM;
		if(piplusREC.E()!=0. && piminusREC.E()!=0.){
			vmREC=piplusREC+piminusREC;
		}
		if(protonRECasifpion.E()!=0. && piminusREC.E()!=0.){
		  vmRECfromproton=protonRECasifpion+piminusREC;
		}
		
		//pion reconstruction efficiency
		double thispipeff = (isPiPlusFound) ? 1.0 : 0.0;
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
		
		//VM rec
		if(vmREC.E()==0 && vmRECfromproton.E()==0) continue;
		double rho_mass = vmREC.M();
		double rho_mass_fromproton = vmRECfromproton.M();
		h_VM_mass_REC->Fill(rho_mass);
		h_VM_mass_REC->Fill(rho_mass_fromproton);
		if(piplusMC.Theta()>0.009  && piplusMC.Theta()<0.013 &&
		                  piminusMC.Theta()>0.009 && piminusMC.Theta()<0.013 ) h_VM_mass_REC_etacut->Fill(vmREC.M());
		if(isPiMinusFound && isPiPlusFound){ 
		  h_VM_mass_REC_justpions->Fill(rho_mass);
		}
	}
	output->Write();
	output->Close();
	
	return 0;
}
