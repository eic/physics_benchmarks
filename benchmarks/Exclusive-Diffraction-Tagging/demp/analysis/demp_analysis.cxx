#include <fstream>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/RotationY.h>
#include <TMath.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

void demp_analysis(const std::string& config_name)
{
  
  //--------------------------------------------------------------------------------------------------------------------------------------------
  
  // Read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const int         ebeam         = config["ebeam"];
  const int         pbeam         = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DEMP analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix for histograms: {}\n", output_prefix);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);
   
  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  // Set output file for the histograms
  std::string output_name_hists = fmt::format("{}.root", output_prefix);
  cout << "Output file for histograms = " << output_name_hists << endl;
  TFile* ofile = new TFile(output_name_hists.c_str(), "RECREATE");

  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  // Set up input file chain
  TChain *mychain = new TChain("events");
  mychain->Add(rec_file.c_str());

    
  //--------------------------------------------------------------------------------------------------------------------------------------------
  
  // Initialize reader
  TTreeReader tree_reader(mychain);

  // Get weight information
  TTreeReaderArray<std::string> weight_keys(tree_reader,"GPStringKeys");
  TTreeReaderArray<std::vector<std::string>> weight_values(tree_reader,"GPStringValues");

  // Get Particle Information
  TTreeReaderArray<int>    partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int>    partPdg(tree_reader, "MCParticles.PDG");
  TTreeReaderArray<double> partMass(tree_reader,"MCParticles.mass");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader,"ReconstructedChargedParticles.momentum.x"); 
  TTreeReaderArray<float> trackMomY(tree_reader,"ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader,"ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<float> trackEng(tree_reader,"ReconstructedChargedParticles.energy");
  TTreeReaderArray<int>   trackPdg(tree_reader,"ReconstructedChargedParticles.PDG");
  TTreeReaderArray<float> trackMass(tree_reader,"ReconstructedChargedParticles.mass");
  TTreeReaderArray<float> trackCharge(tree_reader,"ReconstructedChargedParticles.charge");

  // ZDC Neutrons
  TTreeReaderArray<float> neutEng(tree_reader, "ReconstructedFarForwardZDCNeutrals.energy");
  TTreeReaderArray<float> neutMomX(tree_reader, "ReconstructedFarForwardZDCNeutrals.momentum.x");
  TTreeReaderArray<float> neutMomY(tree_reader, "ReconstructedFarForwardZDCNeutrals.momentum.y");
  TTreeReaderArray<float> neutMomZ(tree_reader, "ReconstructedFarForwardZDCNeutrals.momentum.z");
  TTreeReaderArray<unsigned int> neutClus(tree_reader, "ReconstructedFarForwardZDCNeutrals.clusters_end");

  // ZDC SiPM-on-tile HCal
  TTreeReaderArray<float> neutPosX_hcal(tree_reader, "HcalFarForwardZDCClusters.position.x");
  TTreeReaderArray<float> neutPosY_hcal(tree_reader, "HcalFarForwardZDCClusters.position.y");
  TTreeReaderArray<float> neutPosZ_hcal(tree_reader, "HcalFarForwardZDCClusters.position.z");
  TTreeReaderArray<float> neutEng_hcal(tree_reader, "HcalFarForwardZDCClusters.energy");
  
  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  // Define Histograms

  TH2* eTruthw_Thetap  = new TH2D("eTruthw_Thetap","e' truth #theta vs P; #theta (deg); P (GeV/c); Rate/bin (Hz)",100,120,165,100,4,6);
  TH2* piTruthw_Thetap = new TH2D("piTruthw_Thetap","#pi^{+} truth #theta vs P; #theta (deg); P (GeV/c); Rate/bin (Hz)",100,0,60,100,0,40);
  TH2* nTruthw_Thetap  = new TH2D("nTruthw_Thetap","n truth #theta vs P; #theta (Deg); P (GeV/c); Rate/bin (Hz)",100,0.0,3.5,100,0,45); 
  TH2* nTruthw_Thetaphi  = new TH2D("nTruthw_Thetaphi","n truth #theta vs #phi; #theta (mRad); #phi (deg); Rate/bin (Hz)",100,0.0,62.0,100,-200,200);
  TH2* nTruthw_rot_Thetap  = new TH2D("nTruthw_rot_Thetap","n truth #theta* vs P around p axis; #theta* (Deg); P (GeV/c); Rate/bin (Hz)",100,0.0,2.0,100,0,45);
  TH2* nTruthw_rot_Thetaphi  = new TH2D("nTruthw_rot_Thetaphi","n truth #theta* vs #phi* around p axis; #theta* (mRad); #phi* (deg); Rate/bin (Hz)",100,0.0,35.0,100,-200,200);

  TH2* eRecw_Thetap  = new TH2D("eRecw_Thetap","e' rec #theta vs P; #theta (deg); P (GeV/c); Rate/bin (Hz)",100,120,170,100,0,8);
  TH2* eRecw_Thetaphi  = new TH2D("eRecw_Thetaphi","e' rec #theta vs #phi; #theta (deg); P (GeV/c); Rate/bin (Hz)",100,135,170,100,-200,200);
  TH2* piRecw_Thetap = new TH2D("piRecw_Thetap","#pi^{+} rec #theta vs P; #theta (deg); P (GeV/c); Rate/bin (Hz)",100,0,60,100,0,40);
  TH2* piRecw_Thetaphi  = new TH2D("piRecw_Thetaphi","#pi^{+} rec #theta vs #phi; #theta (deg); Rate/bin (Hz); P (GeV/c)",100,0,55,100,-200,200);
  //TH2* eRecw_Thetap  = new TH2D("eRecw_Thetap","e^{+}' rec #theta vs P; #theta (deg); P (GeV/c)",100,0,60,100,0,15);  // positron
  //TH2* piRecw_Thetap = new TH2D("piRecw_Thetap","#pi^{-} rec #theta vs P; #theta (deg); P (GeV/c)",100,130,170,100,0,10); // pion -
  TH2* nRecw_Thetap  = new TH2D("nRecw_Thetap","n rec #theta vs P; #theta (Deg); P (GeV/c); Rate/bin (Hz)",100,0.8,2.0,100,0,60); 
  TH2* nRecw_Thetaphi  = new TH2D("nRecw_Thetaphi","n rec #theta vs #phi; #theta (mRad); #phi (deg); Rate/bin (Hz)",100,15.0,35.0,100,-200,200);
  TH2* nRecw_rot_Thetaphi  = new TH2D("nRecw_rot_Thetaphi","n rec #theta* vs #phi* around p axis; #theta* (mRad); #phi* (deg); Rate/bin (Hz)",100,0.0,10.0,100,-200,200);
  TH2* nRecw_rot_PosXY  = new TH2D("nRecw_rot_PosXY","n X vs Y around proton axis at Z = 35 m ( #theta* < 4.0 mRad, E > 10 GeV ); x (mm); y (mm); Rate/bin (Hz)",100,-200,200,100,-200,200);
  TH2* nRecw_rot_Thetap  = new TH2D("nRecw_rot_Thetap","n rec #theta* vs P around p axis ( #theta* < 4.0 mRad, E > 10 GeV ); #theta* (mRad); P (GeV/c); Rate/bin (Hz)",100,0.0,4.0,100,5,60);

  TH1* nRec_en = new TH1D("nRec_en", "n rec E for all clusters ( #theta* < 4.0 mRad ); E (GeV); Rate (Hz)", 100, 0.0, 60);
  nRec_en->SetLineWidth(2);
  TH1* nRec_clus = new TH1D("nRec_clus", "n clusters ( #theta* < 4.0 mRad )", 100, 0.0, 8.0);
  nRec_clus->SetLineWidth(2);

  // Neutron theta-phi plots 
  TH1* n_ThetaDiff = new TH1D("n_ThetaDiff", "#theta*_{pMiss_rec} - #theta*_{ZDC}; #theta*_{pMiss_rec} - #theta*_{ZDC}(Deg); Rate (Hz)", 100, -0.3, 1.5);
  n_ThetaDiff->SetLineWidth(2);
  TH1* n_PhiDiff = new TH1D("n_PhiDiff", " #phi*_{pMiss_rec} - #phi*_{ZDC}; #phi*_{pMiss_rec} - #phi*_{ZDC}(Deg); Rate (Hz)", 100, -50, 50);
  n_PhiDiff->SetLineWidth(2);
  TH2* n_ThetaPhiDiff = new TH2D("n_ThetaPhiDiff", "#theta*_{pMiss_rec} - #theta*_{ZDC} vs #phi*_{pMiss_rec} - #phi*_{ZDC}; #theta*_{pMiss_rec} - #theta*_{ZDC} (Deg); #phi*_{pMiss_rec} - #phi*_{ZDC} (Deg); Rate/bin (Hz)",100, -1.0, 1.0, 100, -75, 75);
  TH2* pMissRecw_Thetaphi = new TH2D("pMissRecw_Thetaphi", "pMiss rec #theta vs #phi; #theta (mRad); #phi (deg); Rate/bin (Hz)",100,15.0,35.0,100,-200,200);
  TH2* pMissRecw_rot_Thetaphi = new TH2D("pMissRecw_rot_Thetaphi", "pMiss rec  #theta* vs #phi* around p axis; #theta* (mRad); #phi* (deg); Rate/bin (Hz)",100,0.0,10.0,100,-200,200);
  TH2* n_TruthRecw_ThetaPhiDiff = new TH2D("n_TruthRecw_ThetaPhiDiff", " #theta*_{n_MC} - #theta*_{n_rec} vs #phi*_{n_MC} - #phi*_{n_rec}; #theta*_{n_MC} - #theta*_{n_rec} (Deg); #phi*_{n_MC} - #phi*_{n_rec} (Deg); Rate/bin (Hz)",100, -0.2, 0.2, 100, -25, 25);

  // Absolute difference -t plots
  TH1* htw_t1 = new TH1D("htw_t1", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}); Rate (Hz) ", 100, -0.1,0.1);
  htw_t1->SetLineColor(kBlue); htw_t1->SetLineWidth(2);
  TH1* htw_t2 = new TH1D("htw_t2", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}); Rate (Hz) ", 100, -0.1,0.1);
  htw_t2->SetLineColor(kRed);  htw_t2->SetLineWidth(2);
  TH1* htw_t3 = new TH1D("htw_t3", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}); Rate (Hz) ", 100, -0.1,0.1);
  htw_t3->SetLineColor(kMagenta); htw_t3->SetLineWidth(2);
  TH1* htw_t4 = new TH1D("htw_t4", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}); Rate (Hz) ", 100, -0.1,0.1);
  htw_t4->SetLineColor(kGreen); htw_t4->SetLineWidth(2);

  // -t reconstruction plots
  TH2* htw_rec1 = new TH2D("htw_rec1", "-t rec vs -t truth Distribution; -t_{rec} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec1->SetLineWidth(2);
  TH2* htwz_rec1 = new TH2D("htwz_rec1", "-t rec vs -t truth Distribution; -t_{rec} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec1->SetLineWidth(2); // zoomed version
  TH2* htw_rec2 = new TH2D("htw_rec2", "-t alt_rec vs -t truth Distribution; -t_{alt_rec} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec2->SetLineWidth(2);
  TH2* htwz_rec2 = new TH2D("htwz_rec2", "-t alt_rec vs -t truth Distribution; -t_{alt_rec} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec2->SetLineWidth(2); // zoomed version
  TH2* htw_rec3 = new TH2D("htw_rec3", "-t rec_pT vs -t truth Distribution; -t_{rec_pT} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec3->SetLineWidth(2);
  TH2* htwz_rec3 = new TH2D("htwz_rec3", "-t rec_pT vs -t truth Distribution; -t_{rec_pT} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec3->SetLineWidth(2); // zoomed version
  TH2* htw_rec4 = new TH2D("htw_rec4", "-t rec_corr vs -t truth Distribution; -t_{rec_corr} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec4->SetLineWidth(2);
  TH2* htwz_rec4 = new TH2D("htwz_rec4", "-t rec_corr vs -t truth Distribution; -t_{rec_corr} (GeV^{2});-t_{truth}(GeV^{2}); Rate/bin (Hz)", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec4->SetLineWidth(2); // zoomed version

  // Resolution plots
  TH1* htw_res_e = new TH1D("htw_res_e", "e' Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%); Rate (Hz)", 100, -20, 20);
  htw_res_e->SetLineWidth(2);
  TH1* htw_res_pi = new TH1D("htw_res_pi", "#pi^{+} Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%); Rate (Hz)", 100, -15, 15);
  htw_res_pi->SetLineWidth(2);
  TH1* htw_res_n1 = new TH1D("htw_res_n1", "n Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%); Rate (Hz)", 100, -60, 60);
  htw_res_n1->SetLineWidth(2);
  TH1* htw_res_n2 = new TH1D("htw_res_n2", "n Track #theta* Resolution Distribution (%); (#theta*_{rec} - #theta*_{MC})/#theta*_{MC} (%); Rate (Hz)", 100, -50, 50);
  htw_res_n2->SetLineWidth(2);
  TH1* htw_res_n3 = new TH1D("htw_res_n3", "n Track #phi* Resolution Distribution (%); (#phi*_{rec} - #phi*_{MC})/#phi*_{MC} (%); Rate (Hz)", 100, -50, 50);
  htw_res_n3->SetLineWidth(2);
  TH1* htw_res_n4 = new TH1D("htw_res_n4", "n Track Momentum Resolution Distribution (%); (P_{rec_corr} - P_{MC})/P_{MC} (%); Rate (Hz)", 100, -3, 3);
  htw_res_n4->SetLineWidth(2);

  TH1* htw_res1 = new TH1D("htw_res1", "-t Resolution Distribution (%); (t_{rec} - t_{truth})/t_{truth} (%); Rate (Hz)", 100, -200, 200);
  htw_res1->SetLineWidth(2);
  TH1* htw_res2 = new TH1D("htw_res2", "-t Resolution Distribution (%); (t_{alt_rec} - t_{truth})/t_{truth} (%); Rate (Hz)", 100, -110, 210);
  htw_res2->SetLineWidth(2);
  TH1* htw_res3 = new TH1D("htw_res3", "-t Resolution Distribution (%); (t_{rec_pT} - t_{truth})/t_{truth} (%); Rate (Hz)", 100, -110, 210);
  htw_res3->SetLineWidth(2);
  TH1* htw_res4 = new TH1D("htw_res4", "-t Resolution Distribution (%); (t_{rec_corr} - t_{truth})/t_{truth} (%); Rate (Hz)", 100, -100, 100);
  htw_res4->SetLineWidth(2);
  TH1* htw_res5 = new TH1D("htw_res5", "Q^{2} Resolution Distribution (%); (Q^{2}_{rec} - Q^{2}_{truth})/Q^{2}_{truth} (%); Rate (Hz)", 100, -20, 20);
  htw_res5->SetLineWidth(2);
  TH1* htw_res6 = new TH1D("htw_res6", "W Resolution Distribution (%); (W_{rec} - W_{truth})/W_{truth} (%); Rate (Hz)", 100, -60, 60);
  htw_res6->SetLineWidth(2);

  // Effeciency plots
  TH2* Q2_t_DetEff_Uncut = new TH2F("Q2_t_DetEff_Uncut", "Q^{2}_{truth} vs -t_{truth} for thrown events; Q^{2} (GeV^{2}); -t (GeV^{2}); Rate/bin (Hz)", 10, 0, 40, 10, 0, 1.5);
  TH2* Q2_t_DetEff_Cut = new TH2F("Q2_t_DetEff_Cut", "Q^{2}_{truth} vs -t_{truth} for detected events; Q^{2} (GeV^{2}); -t (GeV^{2}); Rate/bin (Hz)", 10, 0, 40, 10, 0, 1.5);
  TH2* Q2_t_DetEff = new TH2F("Q2_t_DetEff", "Q^{2}_{truth} vs -t_{truth} detected/thrown ratio; Q^{2} (GeV^{2}); -t (GeV^{2})", 10, 0, 40, 10, 0, 1.5);

  TH1* eTruthw_Eta_Uncut = new TH1D("eTruthw_Eta_Uncut", "e' #eta for thrown events; #eta; Rate (Hz)",100,-2,0);
  eTruthw_Eta_Uncut->SetLineWidth(2);
  TH1* eRecw_Eta_Cut = new TH1D("eRecw_Eta_Cut", "e' #eta for detected events; #eta; Rate (Hz)",100,-2,0);
  eRecw_Eta_Cut->SetLineWidth(2);
  TH1* eEff_Eta = new TH1D("eEff_Eta", "e' Tracking efficiency as fn of #eta; #eta; Eff", 100,-2,0);
  eEff_Eta->SetLineWidth(2);

  TH1* piTruthw_Eta_Uncut = new TH1D("piTruthw_Eta_Uncut", "#pi^{+} #eta for thrown events; #eta; Rate (Hz)",100,0,4);
  piTruthw_Eta_Uncut->SetLineWidth(2);
  TH1* piRecw_Eta_Cut = new TH1D("piRecw_Eta_Cut", "#pi^{+} #eta for detected events; #eta; Rate (Hz)",100,0,4);
  piRecw_Eta_Cut->SetLineWidth(2);
  TH1* piEff_Eta = new TH1D("piEff_Eta", "#pi^{+} Tracking efficiency as fn of #eta; #eta; Eff",100,0,4); 
  piEff_Eta->SetLineWidth(2);

  TH1* eTruthw_P_Uncut = new TH1D("eTruthw_P_Uncut", "e' P for thrown events; P (GeV/c); Rate (Hz)",100,4,7);
  eTruthw_P_Uncut->SetLineWidth(2);
  TH1* eRecw_P_Cut = new TH1D("eRecw_P_Cut", "e' P for detected events; P (GeV/c); Rate (Hz)",100,4,7);
  eRecw_P_Cut->SetLineWidth(2);
  TH1* eEff_P = new TH1D("eEff_P", "e' Tracking efficiency as fn of P; P (GeV/c); Eff", 100,4,7);
  eEff_P->SetLineWidth(2);

  TH1* piTruthw_P_Uncut = new TH1D("piTruthw_P_Uncut", "#pi^{+} P for thrown events; P (GeV/c); Rate (Hz)",100,0,30);
  piTruthw_P_Uncut->SetLineWidth(2);
  TH1* piRecw_P_Cut = new TH1D("piRecw_P_Cut", "#pi^{+} P for detected events; P (GeV/c); Rate (Hz)",100,0,30);
  piRecw_P_Cut->SetLineWidth(2);
  TH1* piEff_P = new TH1D("piEff_P", "#pi^{+} Tracking efficiency as fn of P; P (GeV/c); Eff",100,0,30); 
  piEff_P->SetLineWidth(2);

  // Neutrons HCal
  TH2* nRecw_Thetap_hcal  = new TH2D("nRecw_Thetap_hcal","n rec #theta vs P for 1 cluster events; #theta (Deg); P (GeV/c); Rate/bin (Hz)",100,0.8,2.0,100,0,40);
  TH2* nRecw_rot_Thetap_hcal  = new TH2D("nRecw_rot_Thetap_hcal","n rec #theta* vs P around p axis for 1 cluster events; #theta* (mRad); P (GeV/c); Rate/bin (Hz)",100,0,4.0,100,5,40);
  
  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  //Defining the four vectors
  ROOT::Math::PxPyPzEVector elec_beam; // initialized the 4 vector for electron beam
  ROOT::Math::PxPyPzEVector prot_beam; // initialized the 4 vector for proton beam

  ROOT::Math::PxPyPzEVector elec_mc; // initialized the 4 vector for truth electron
  ROOT::Math::PxPyPzEVector pi_mc; // initialized the 4 vector for truth pion
  ROOT::Math::PxPyPzEVector neut_mc; // initialized the 4 vector for truth neutron
  ROOT::Math::PxPyPzEVector neut_rot_mc; // initialized the 4 vector for truth neutron with a rotation of 25 mrad

  ROOT::Math::RotationY rot; // initialized rotation vector
  rot.SetAngle(0.025);

  ROOT::Math::PxPyPzEVector elec_rec; // initialized the 4 vector for reconstructed electron
  ROOT::Math::PxPyPzEVector pi_rec; // initialized the 4 vector for reconstructed pion
  ROOT::Math::PxPyPzEVector neut_rec; // initialized the 4 vector for reconstructed neutron
  ROOT::Math::PxPyPzEVector neut_rot_rec; // initialized the 4 vector for reconstructed neutron with a rotation of 25 mrad

  ROOT::Math::PxPyPzEVector virtphoton_truth; // intialized the 4 vector for truth virtual photon
  ROOT::Math::PxPyPzEVector ttruth; // intialized the 4 vector for ttruth (-t)from first loop
  ROOT::Math::PxPyPzEVector talttruth; // intialized the 4 vector for talttruth(-t) from second loop

  ROOT::Math::PxPyPzEVector virtphoton_rec; //intialized the 4 vector for reconstructed virtual photon
  ROOT::Math::PxPyPzEVector trec; // intialized the 4 vector for trec (-t)from first loop
  ROOT::Math::PxPyPzEVector taltrec; // intialized the 4 vector for taltrec(-t) from second loop
  ROOT::Math::PxPyPzEVector trecpT; // intialized the 4 vector for trecpT(-t)
  ROOT::Math::PxPyPzEVector trecpT_rot; // intialized the 4 vector for trecpT(-t) with a rotation of 25 mrad
  ROOT::Math::PxPyPzEVector p_miss_rec;  //intialized the 4 vector for missing momentum
  ROOT::Math::PxPyPzEVector p_miss_rot_rec; //intialized the 4 vector for missing momentum with a rotation of 25 mrad
  ROOT::Math::PxPyPzMVector neut_corr; // intialized the 4 vector for reconstructed corrected neutron
  ROOT::Math::PxPyPzEVector treccorr; // intialized the 4 vector for trecpT(-t)

  ROOT::Math::XYZVector neut_pos_hcal; // initialized the 3 vector for zdc position
  ROOT::Math::PxPyPzEVector neut_rec_hcal; // initialized the 4 vector for reconstructed neutorn in hcal
  ROOT::Math::PxPyPzEVector neut_rot_rec_hcal; // initialized the 4 vector for reconstructed neutron with a rotation of 25 mrad in hcal

  unsigned int count2 = 0; // counter on neutrons within 4 mrad
  int hcal_clus_size;
  double neut_rec_p_hcal;

  double neutMass = 0.93965420;
  double weight, partEng; // weight and energy of the particles
  double Q2_truth, W_truth, y_truth, t_truth, t_alttruth; // Truth kinematic variables
  double Q2_rec, W_rec, y_rec, t_rec, t_altrec, t_recpT, t_reccorr; // Reconstructed kinematic variables
  double neutPosX, neutPosY; // neutron position
  
  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  // Defining initial colliding beams
  double eMass = 0.000510998950; //electron beam
  double eEng = ebeam;
  double e_pmag = sqrt(pow(eEng,2)-pow(eMass,2));
  double e_p1 = 0.;
  double e_p2 = 0.;
  double e_p3 = -1*e_pmag;
  elec_beam.SetPxPyPzE(e_p1, e_p2, e_p3, eEng); 
 
  double pMass = 0.93827208816; // proton beam
  double pEng = pbeam;
  double p_pmag = sqrt(pow(pEng,2)-pow(pMass,2));
  double c_a = 0.025;
  double p_p1 = -p_pmag*sin(c_a);
  double p_p2 = 0.;
  double p_p3 = p_pmag*cos(c_a);
  prot_beam.SetPxPyPzE(p_p1, p_p2, p_p3, pEng);

  //--------------------------------------------------------------------------------------------------------------------------------------------
 
  bool x,y,z; // x,y, and z are for reconstructed electron, pion, and neutron
 
  while(tree_reader.Next()) { // Loop over events
 
    x = false, y = false, z = false;
 
    std::string weight_name = weight_keys[0]; // accessing weights of particles
    std::vector<std::string> weight_value = weight_values[0];
    weight = std::stod( *(weight_value.begin()) );

    
    //--------------------------------------------------------------------------------------------------------------------------------------------
    
    for(unsigned int i=0; i<partGenStat.GetSize(); i++) { // Loop over thrown particles
      partEng = sqrt(pow(partMomX[i],2) + pow(partMomY[i],2) + pow(partMomZ[i],2) + pow(partMass[i],2)); // Energy of all Monte Carlo particles
		
      if(partGenStat[i] == 1 && partPdg[i] == 11) { // Select stable thrown particles and look at electron
	elec_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
	eTruthw_Thetap -> Fill(elec_mc.Theta()*TMath::RadToDeg(), elec_mc.P(), weight);
	eTruthw_Eta_Uncut -> Fill(elec_mc.Eta(), weight);
	eTruthw_P_Uncut -> Fill(elec_mc.P(), weight);
      }
 
      if(partGenStat[i] == 1 && partPdg[i] == 211) { // Look at pion
	pi_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
	piTruthw_Thetap -> Fill(pi_mc.Theta()*TMath::RadToDeg(), pi_mc.P(), weight);
	piTruthw_Eta_Uncut -> Fill(pi_mc.Eta(), weight);
	piTruthw_P_Uncut -> Fill(pi_mc.P(), weight);
      }
 
      if(partGenStat[i] == 1 && partPdg[i] == 2112) { // Look at neutron
	neut_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], partEng);
	nTruthw_Thetap -> Fill(neut_mc.Theta()*TMath::RadToDeg(), neut_mc.P(), weight);
	nTruthw_Thetaphi -> Fill(neut_mc.Theta()*1000., neut_mc.Phi()*TMath::RadToDeg(), weight);
 
	neut_rot_mc = rot*neut_mc;  // rotate w.r.t to proton axis
	nTruthw_rot_Thetap -> Fill(neut_rot_mc.Theta()*TMath::RadToDeg(), neut_rot_mc.P(), weight);
	nTruthw_rot_Thetaphi -> Fill(neut_rot_mc.Theta()*1000., neut_rot_mc.Phi()*TMath::RadToDeg(), weight);
      }			           
 
    } // for over thrown particles
 
    
    //--------------------------------------------------------------------------------------------------------------------------------------------
    
    for(unsigned int i=0; i<trackPdg.GetSize(); i++) { // Loop over reconstructed particles 
      // if(trackPdg[i] == 11) { // Look at electron
      if(trackCharge[i] == -1 && trackMomZ[i] < 0) { 
	x = true;
	elec_rec.SetPxPyPzE(trackMomX[i],trackMomY[i],trackMomZ[i], trackEng[i]);
	eRecw_Thetap -> Fill(elec_rec.Theta()*TMath::RadToDeg(), elec_rec.P(), weight);
	eRecw_Thetaphi -> Fill(elec_rec.Theta()*TMath::RadToDeg(), elec_rec.Phi()*TMath::RadToDeg(), weight);
	eRecw_Eta_Cut -> Fill(elec_mc.Eta(), weight);
	eRecw_P_Cut -> Fill(elec_mc.P(), weight);
      }
 
      // if(trackPdg[i] == 211) { // Look at pion
      if(trackCharge[i] == +1 && trackMomZ[i] > 0) {
	y = true;
	pi_rec.SetPxPyPzE(trackMomX[i],trackMomY[i],trackMomZ[i], trackEng[i]);
	piRecw_Thetap -> Fill(pi_rec.Theta()*TMath::RadToDeg(), pi_rec.P(), weight);
	piRecw_Thetaphi ->  Fill(pi_rec.Theta()*TMath::RadToDeg(), pi_rec.Phi()*TMath::RadToDeg(), weight);
	piRecw_Eta_Cut -> Fill(pi_mc.Eta(), weight);
	piRecw_P_Cut -> Fill(pi_mc.P(), weight);
      }
 
    }// for over reconstructed particles

    //--------------------------------------------------------------------------------------------------------------------------------------------
    
    for(unsigned int i=0; i<neutEng.GetSize(); i++) { // Loop over zdc neutrons
 
      neut_rec.SetPxPyPzE(neutMomX[i],neutMomY[i],neutMomZ[i], neutEng[i]);
      nRecw_Thetap -> Fill(neut_rec.Theta()*TMath::RadToDeg(), neut_rec.P(), weight);
      nRecw_Thetaphi -> Fill(neut_rec.Theta()*1000., neut_rec.Phi()*TMath::RadToDeg(), weight);
    
      neut_rot_rec = rot*neut_rec; // rotate w.r.t to proton axis
      nRecw_rot_Thetaphi -> Fill(neut_rot_rec.Theta()*1000., neut_rot_rec.Phi()*TMath::RadToDeg(), weight);
      
      if(neut_rot_rec.Theta()*1000. < 4.0){ // acceptance of the zdc
        nRec_clus -> Fill(neutClus[i]);
	nRec_en -> Fill(neut_rot_rec.E(), weight);
	
   	if(neut_rot_rec.E()>10.0){ // neutron energy cut
	  z = true;
	  nRecw_rot_Thetap -> Fill(neut_rot_rec.Theta()*1000., neut_rot_rec.P(), weight);

	  neutPosX = 35000 * sin(neut_rot_rec.Theta()) * cos(neut_rot_rec.Phi()); // neutron position at r = z = 35.0 m
      	  neutPosY = 35000 * sin(neut_rot_rec.Theta()) * sin(neut_rot_rec.Phi());
          nRecw_rot_PosXY -> Fill(neutPosX, neutPosY, weight);
	}
 
      }
 
    }// for over zdc neutrons
 
    //--------------------------------------------------------------------------------------------------------------------------------------------
 
    for(unsigned int i=0; i<neutEng.GetSize(); i++) { // Loop over zdc neutrons in HCal
 
      hcal_clus_size = neutEng_hcal.GetSize(); //ZDC HCal cluster size -> No. of clusters in ZDC
 
      if(hcal_clus_size >0 ){ // Selected the events correspond to on clusters
 
	neut_pos_hcal.SetXYZ(neutPosX_hcal[0], neutPosY_hcal[0], neutPosZ_hcal[0]);
                        
	neut_rec_p_hcal = std::sqrt(pow(neutEng_hcal[0],2)- pow(neutMass,2)); // neutrons momentum
						
	neut_rec_hcal.SetPxPyPzE(neut_rec_p_hcal * sin(neut_pos_hcal.Theta()) * cos(neut_pos_hcal.Phi()), 
				 neut_rec_p_hcal * sin(neut_pos_hcal.Theta()) * sin(neut_pos_hcal.Phi()),   
				 neut_rec_p_hcal * cos(neut_pos_hcal.Theta()),
				 neutEng_hcal[i]);
	nRecw_Thetap_hcal -> Fill(neut_rec_hcal.Theta()*TMath::RadToDeg(), neut_rec_hcal.P(), weight);					    

	neut_rot_rec_hcal = rot*neut_rec_hcal; // rotate w.r.t to proton axis						
 
	if(neut_rot_rec_hcal.Theta()*1000. < 4.0 && neut_rot_rec_hcal.E()> 10.0){ 	
	  nRecw_rot_Thetap_hcal -> Fill(neut_rot_rec_hcal.Theta()*1000., neut_rot_rec_hcal.P(), weight);
	}
      }
 
    }// for over zdc neutrons in HCal

    //--------------------------------------------------------------------------------------------------------------------------------------------
    
    // Truth kinematic variables
    virtphoton_truth = (elec_beam - elec_mc);
    Q2_truth = -1*(virtphoton_truth.mag2());
    W_truth = ((virtphoton_truth + prot_beam).mag());
    y_truth = (prot_beam.Dot(virtphoton_truth))/(prot_beam.Dot(elec_beam)); // Energy Loss y 
    
    ttruth = (virtphoton_truth - pi_mc); 
    t_truth = -1*(ttruth.mag2()); // ttruth is the -t from the first loop
                                 
    talttruth = (prot_beam - neut_mc); 
    t_alttruth = -1*(talttruth.mag2()); // t_alttruth is the -t from the second loop
 
    // Efficiency plots
    if(neut_rot_mc.Theta()*1000. < 4.0){
      count2++;
      if(Q2_truth > 5 && Q2_truth < 35){
	Q2_t_DetEff_Uncut -> Fill(Q2_truth, t_truth, weight);
      }
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------------------
   
    // Reconstructed kinematic variables
    if (x == true && y == true && z == true ){ // if e', pi, and neutron are in coincidence
 
      virtphoton_rec = (elec_beam - elec_rec);
      Q2_rec = -1*(virtphoton_rec.mag2()); 
      W_rec = ((virtphoton_rec + prot_beam).mag()); 
      y_rec =  (prot_beam.Dot(virtphoton_rec))/(prot_beam.Dot(elec_beam)); // Energy Loss y

      if(Q2_rec > 5 && Q2_rec < 35){ // Q2 Cut 
 
	// t-method plots
	trec = (virtphoton_rec - pi_rec); // First method to reconstruct -t // No change in values after rotation.
	t_rec = -1*(trec.mag2()); // t_rec is the -t from the first loop 
	htw_rec1 -> Fill(t_rec, t_truth, weight); 
	htwz_rec1 -> Fill(t_rec, t_truth, weight); // zoomed version
                                 
	taltrec = (prot_beam - neut_rec); // Second method to reconstruct -t // No change in values after rotation.
	t_altrec = -1*(taltrec.mag2()); // t_altrec is the -t from the second loop
	htw_rec2 -> Fill(t_altrec, t_truth, weight); 
	htwz_rec2 -> Fill(t_altrec, t_truth, weight); // zoomed version
 
	trecpT = (elec_rec +  pi_rec); // Third method to reconstruct -t // Values changed after rotation.
	trecpT_rot = rot*trecpT;
	t_recpT = trecpT_rot.Perp2();
	htw_rec3 -> Fill(t_recpT, t_truth, weight); 
	htwz_rec3 -> Fill(t_recpT, t_truth, weight); // zoomed version
 
	p_miss_rec = (elec_beam + prot_beam) - (elec_rec + pi_rec) ; // Defined missing momentum information -> Fourth method to reconstruct -t
	pMissRecw_Thetaphi -> Fill(p_miss_rec.Theta()*1000., p_miss_rec.Phi()*TMath::RadToDeg(), weight);
	
	p_miss_rot_rec = rot*p_miss_rec; // rotate p_miss_rec w.r.t to proton axis
	pMissRecw_rot_Thetaphi -> Fill(p_miss_rot_rec.Theta()*1000., p_miss_rot_rec.Phi()*TMath::RadToDeg(), weight);
	
        neut_corr.SetCoordinates( p_miss_rec.P() * sin(neut_rec.Theta()) * cos(neut_rec.Phi()), 
				  p_miss_rec.P() * sin(neut_rec.Theta()) * sin(neut_rec.Phi()), 
				  p_miss_rec.P() * cos(neut_rec.Theta()),
				  neutMass );
 
	treccorr = (prot_beam - neut_corr); // No change in values after rotation.
	t_reccorr = -1*(treccorr.mag2());
	htw_rec4 -> Fill(t_reccorr, t_truth, weight); 
	htwz_rec4 -> Fill(t_reccorr, t_truth, weight); // zoomed version
 
	// Absolute difference -t plots
	htw_t1 -> Fill(t_rec - t_truth);
	htw_t2 -> Fill(t_altrec - t_truth);
	htw_t3 -> Fill(t_recpT - t_truth);
	htw_t4 -> Fill(t_reccorr - t_truth);
	
	// Neutron theta-phi plots 
	n_ThetaDiff -> Fill((p_miss_rot_rec.Theta() - neut_rot_rec.Theta())*TMath::RadToDeg(), weight);
	n_PhiDiff -> Fill((p_miss_rot_rec.Phi() - neut_rot_rec.Phi())*TMath::RadToDeg(), weight);
	n_ThetaPhiDiff -> Fill((p_miss_rot_rec.Theta() - neut_rot_rec.Theta())*TMath::RadToDeg(), (p_miss_rot_rec.Phi() - neut_rot_rec.Phi())*TMath::RadToDeg(), weight);
	n_TruthRecw_ThetaPhiDiff -> Fill((neut_rot_mc.Theta() - neut_rot_rec.Theta())*TMath::RadToDeg(), (neut_rot_mc.Phi() - neut_rot_rec.Phi())*TMath::RadToDeg(), weight); 
 
	// Resolution plots
	htw_res_e -> Fill((elec_rec.P() - elec_mc.P())/(elec_mc.P())*100, weight); 
	htw_res_pi -> Fill((pi_rec.P() - pi_mc.P())/(pi_mc.P())*100, weight);
	htw_res_n1 -> Fill((neut_rot_rec.P() - neut_rot_mc.P())/(neut_rot_mc.P())*100, weight);
	htw_res_n4 -> Fill((neut_corr.P() - neut_mc.P())/(neut_mc.P())*100, weight);
	htw_res_n2 -> Fill((neut_rot_rec.Theta() - neut_rot_mc.Theta())/(neut_rot_mc.Theta())*100, weight);
	htw_res_n3 -> Fill((neut_rot_rec.Phi() - neut_rot_mc.Phi())/(neut_rot_mc.Phi())*100, weight);
	
	htw_res1 -> Fill((t_rec - t_truth)/(t_truth)*100, weight); 
	htw_res2 -> Fill((t_altrec - t_truth)/(t_truth)*100, weight);
	htw_res3 -> Fill((t_recpT - t_truth)/(t_truth)*100, weight);
	htw_res4 -> Fill((t_reccorr - t_truth)/(t_truth)*100, weight);
	htw_res5 -> Fill((Q2_rec - Q2_truth)/(Q2_truth)*100, weight);
	htw_res6 -> Fill((W_rec - W_truth)/(W_truth)*100, weight);
 
	// Efficiency plots
	Q2_t_DetEff_Cut -> Fill(Q2_truth, t_truth, weight);
	
      } //Q2 cut
    } // if over x,y, and z
  } // End of event loop (while loop)
 
  // Efficiency plots
  Q2_t_DetEff -> Divide(Q2_t_DetEff_Cut, Q2_t_DetEff_Uncut, 1, 1, "b");
  eEff_Eta -> Divide(eRecw_Eta_Cut, eTruthw_Eta_Uncut, 1, 1, "b");
  piEff_Eta -> Divide(piRecw_Eta_Cut, piTruthw_Eta_Uncut, 1, 1, "b");
  eEff_P -> Divide(eRecw_P_Cut, eTruthw_P_Uncut, 1, 1, "b");
  piEff_P -> Divide(piRecw_P_Cut, piTruthw_P_Uncut, 1, 1, "b");
 
  cout<<"truth_neutron =  "<<count2<<endl;

  //--------------------------------------------------------------------------------------------------------------------------------------------

  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
    
}

