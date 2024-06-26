#include <TCanvas.h>
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
#include <TStyle.h>
#include <TLegend.h>

void demp_analysis()
{

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Set output file for the histograms
  TFile *ofile = TFile::Open("pion.hist.root","RECREATE");

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Set up input file chain
  TChain *mychain = new TChain("events");
  unsigned int count1 = 0; // counter on no. of files
  unsigned int count2 = 0; // counter on neutrons within 4 mrad

  unsigned int nFiles1 = 1;
  unsigned int nFiles2 = 2;
  unsigned int nFiles =  nFiles1 * nFiles2; // Total number of files to assign weights (or weight correction)

  for(unsigned int i = 0; i < nFiles1; i++) {
    for(unsigned int j = 0; j < nFiles2; j++) {
      count1++;
      mychain->Add(Form("Test_5on41_pi+_Reco_AllEvents_%i_%i.root", 1 + i, 5000 + j * 5000));
     // cout << "  " << Form("Test_5on41_pi+_Reco_AllEvents_%i_%i.root", 1 + i, 5000 + j * 5000) << endl;
    }
  }
 // cout<<"count1 =  "<<count1<<", nFiles =  "<<nFiles<<endl;

  if (count1 != nFiles) {
    std::cout << "Error in counting the files" << std::endl;
    std::exit(EXIT_FAILURE);  // Exit with a non-zero status to indicate failure
  }

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Initialize reader
  TTreeReader tree_reader(mychain);

  // Get weight information
  TTreeReaderArray<pair<string,vector<string>>> weight_map(tree_reader,"_stringMap"); 

  // Get Particle Information
  TTreeReaderArray<int>    partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<float>  partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<float>  partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<float>  partMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int>    partPdg(tree_reader, "MCParticles.PDG");
  TTreeReaderArray<double> partMass(tree_reader,"MCParticles.mass");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader,"ReconstructedChargedParticles.momentum.x"); 
  TTreeReaderArray<float> trackMomY(tree_reader,"ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader,"ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<float> trackEng(tree_reader,"ReconstructedChargedParticles.energy");
  TTreeReaderArray<int>   trackPdg(tree_reader,"ReconstructedChargedParticles.PDG");
  TTreeReaderArray<float> trackMass(tree_reader,"ReconstructedChargedParticles.mass");

  //ZDC Neutrons
  TTreeReaderArray<float> neutEng(tree_reader, "ReconstructedFarForwardZDCNeutrons.energy");
  TTreeReaderArray<float> neutMomX(tree_reader, "ReconstructedFarForwardZDCNeutrons.momentum.x");
  TTreeReaderArray<float> neutMomY(tree_reader, "ReconstructedFarForwardZDCNeutrons.momentum.y");
  TTreeReaderArray<float> neutMomZ(tree_reader, "ReconstructedFarForwardZDCNeutrons.momentum.z");
  TTreeReaderArray<unsigned int> neutClus(tree_reader, "ReconstructedFarForwardZDCNeutrons.clusters_end");

  //ZDC SiPM-on-tile HCal
  TTreeReaderArray<float> neutPosX_hcal(tree_reader, "HcalFarForwardZDCClusters.position.x");
  TTreeReaderArray<float> neutPosY_hcal(tree_reader, "HcalFarForwardZDCClusters.position.y");
  TTreeReaderArray<float> neutPosZ_hcal(tree_reader, "HcalFarForwardZDCClusters.position.z");
  TTreeReaderArray<float> neutEng_hcal(tree_reader, "HcalFarForwardZDCClusters.energy");

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Define Histograms

  TH2* eTruthw_Thetap  = new TH2D("eTruthw_Thethap","e' truth #theta vs P; #theta (deg); P (GeV/c)",100,120,165,100,4,6);
  TH2* piTruthw_Thetap = new TH2D("piTruthw_Thethap","#pi^{+} truth #theta vs P; #theta (deg); P (GeV/c)",100,0,60,100,0,40);
  TH2* nTruthw_Thetap  = new TH2D("nTruthw_Thethap","n truth #theta vs P; #theta (Deg); P (GeV/c)",100,0.0,3.5,100,0,45); 
  TH2* nTruthw_Thetaphi  = new TH2D("nTruthw_Thethaphi","n truth #theta vs #phi; #theta (mRad); #phi (deg)",100,0.0,62.0,100,-200,200);
  TH2* nTruthw_rot_Thetap  = new TH2D("nTruthw_rot_Thethap","n truth #theta* vs P around p axis; #theta* (Deg); P (GeV/c)",100,0.0,2.0,100,0,45);
  TH2* nTruthw_rot_Thetaphi  = new TH2D("nTruthw_rot_Thethaphi","n truth #theta* vs #phi* around p axis; #theta* (mRad); #phi* (deg)",100,0.0,35.0,100,-200,200);

  TH2* eRecw_Thetap  = new TH2D("eRecw_Thethap","e' rec #theta vs P; #theta (deg); P (GeV/c)",100,120,170,100,0,8);
  TH2* piRecw_Thetap = new TH2D("piRecw_Thethap","#pi^{+} rec #theta vs P; #theta (deg); P (GeV/c)",100,0,60,100,0,40);
  //TH2* eRecw_Thetap  = new TH2D("eRecw_Thethap","e^{+}' rec #theta vs P; #theta (deg); P (GeV/c)",100,0,60,100,0,15);  // positron
  //TH2* piRecw_Thetap = new TH2D("piRecw_Thethap","#pi^{-} rec #theta vs P; #theta (deg); P (GeV/c)",100,130,170,100,0,10); // pion -
  TH2* nRecw_Thetap  = new TH2D("nRecw_Thethap","n rec #theta vs P; #theta (Deg); P (GeV/c)",100,0.8,2.0,100,0,60); 
  TH2* nRecw_Thetaphi  = new TH2D("nRecw_Thethaphi","n rec #theta vs #phi; #theta (mRad); #phi (deg)",100,15.0,35.0,100,-200,200);
  TH2* nRecw_rot_Thetap  = new TH2D("nRecw_rot_Thethap","n rec #theta* vs P around p axis ( #theta* < 4.0 mRad, E > 10 GeV ); #theta* (mRad); P (GeV/c)",100,0.0,4.0,100,5,60);
  TH2* nRecw_rot_Thetaphi  = new TH2D("nRecw_rot_Thethaphi","n rec #theta* vs #phi* around p axis; #theta* (mRad); #phi* (deg)",100,0.0,10.0,100,-200,200);

  TH1* nRec_en = new TH1D("nRec_en", "n rec E ( #theta* < 4.0 mRad ); E (GeV)", 100, 0.0, 60);
  nRec_en->SetLineWidth(2);
  TH1* nRec_clus = new TH1D("nRec_clus", "n clusters ( #theta* < 4.0 mRad )", 100, 0.0, 8.0);
  nRec_clus->SetLineWidth(2);

  // Neutron theta-phi plots 
  TH1* n_ThetaDiff = new TH1D("n_ThetaDiff", "#theta_{pMiss_rec} - #theta_{ZDC}; #theta_{pMiss_rec} - #theta_{ZDC}(Deg)", 100, -1, 1);
  n_ThetaDiff->SetLineWidth(2);
  TH1* n_PhiDiff = new TH1D("n_PhiDiff", " #phi_{pMiss_rec} - #phi_{ZDC}; #phi_{pMiss_rec} - #phi_{ZDC}(Deg)", 100, -50, 50);
  n_PhiDiff->SetLineWidth(2);
  TH2* n_ThetaPhiDiff = new TH2D("n_ThetaPhiDiff", "#theta_{pMiss_rec} - #theta_{ZDC} vs #phi_{pMiss_rec} - #phi_{ZDC}; #theta_{pMiss_rec} - #theta_{ZDC} (Deg); #phi_{pMiss_rec} - #phi_{ZDC}(Deg)",100, -1, 1, 100, -25, 25);

  // Absolute difference -t plots
  TH1* htw_t1 = new TH1D("htw_t1", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}) ", 100, -0.1,0.1);
  htw_t1->SetLineColor(kBlue); htw_t1->SetLineWidth(2);
  TH1* htw_t2 = new TH1D("htw_t2", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}) ", 100, -0.1,0.1);
  htw_t2->SetLineColor(kRed);  htw_t2->SetLineWidth(2);
  TH1* htw_t3 = new TH1D("htw_t3", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}) ", 100, -0.1,0.1);
  htw_t3->SetLineColor(kMagenta); htw_t3->SetLineWidth(2);
  TH1* htw_t4 = new TH1D("htw_t4", "-t_{rec, alt_rec, rec_pT, rec_corr} - -t_{truth} Distribution; #Delta -t (GeV^{2}) ", 100, -0.1,0.1);
  htw_t4->SetLineColor(kGreen); htw_t4->SetLineWidth(2);

  // -t reconstruction plots
  TH2* htw_rec1 = new TH2D("htw_rec1", "-t rec vs -t truth Distribution; -t_{rec} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec1->SetLineWidth(2);
  TH2* htwz_rec1 = new TH2D("htwz_rec1", "-t rec vs -t truth Distribution; -t_{rec} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec1->SetLineWidth(2); // zoomed version
  TH2* htw_rec2 = new TH2D("htw_rec2", "-t alt_rec vs -t truth Distribution; -t_{alt_rec} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec2->SetLineWidth(2);
  TH2* htwz_rec2 = new TH2D("htwz_rec2", "-t alt_rec vs -t truth Distribution; -t_{alt_rec} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec2->SetLineWidth(2); // zoomed version
  TH2* htw_rec3 = new TH2D("htw_rec3", "-t rec_pT vs -t truth Distribution; -t_{rec_pT} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec3->SetLineWidth(2);
  TH2* htwz_rec3 = new TH2D("htwz_rec3", "-t rec_pT vs -t truth Distribution; -t_{rec_pT} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec3->SetLineWidth(2); // zoomed version
  TH2* htw_rec4 = new TH2D("htw_rec4", "-t rec_corr vs -t truth Distribution; -t_{rec_corr} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,2.0,100, 0.0,1.5);
  htw_rec4->SetLineWidth(2);
  TH2* htwz_rec4 = new TH2D("htwz_rec4", "-t rec_corr vs -t truth Distribution; -t_{rec_corr} (GeV^{2});-t_{truth}(GeV^{2})", 100, 0.0,0.2,100, 0.0,0.2);
  htwz_rec4->SetLineWidth(2); // zoomed version

  // Resolution plots
  TH1* htw_res_e = new TH1D("htw_res_e", "e' Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%)", 100, -20, 20);
  htw_res_e->SetLineWidth(2);
  TH1* htw_res_pi = new TH1D("htw_res_pi", "#pi^{+} Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%)", 100, -15, 15);
  htw_res_pi->SetLineWidth(2);
  TH1* htw_res_n1 = new TH1D("htw_res_n1", "n Track Momentum Resolution Distribution (%); (P_{rec} - P_{MC})/P_{MC} (%)", 100, -60, 60);
  htw_res_n1->SetLineWidth(2);
  TH1* htw_res_n2 = new TH1D("htw_res_n2", "n Track #theta Resolution Distribution (%); (#theta*_{rec} - #theta*_{MC})/#theta*_{MC} (%)", 100, -50, 50);
  htw_res_n2->SetLineWidth(2);
  TH1* htw_res_n3 = new TH1D("htw_res_n3", "n Track #phi Resolution Distribution (%); (#phi_{rec} - #phi_{MC})/#phi_{MC} (%)", 100, -50, 50);
  htw_res_n3->SetLineWidth(2);
  TH1* htw_res_n4 = new TH1D("htw_res_n4", "n Track Momentum Resolution Distribution (%); (P_{rec_corr} - P_{MC})/P_{MC} (%)", 100, -3, 3);
  htw_res_n4->SetLineWidth(2);

  TH1* htw_res1 = new TH1D("htw_res1", "-t Resolution Distribution (%); (t_{rec} - t_{truth})/t_{truth} (%)", 100, -200, 200);
  htw_res1->SetLineWidth(2);
  TH1* htw_res2 = new TH1D("htw_res2", "-t Resolution Distribution (%); (t_{alt_rec} - t_{truth})/t_{truth} (%)", 100, -110, 210);
  htw_res2->SetLineWidth(2);
  TH1* htw_res3 = new TH1D("htw_res3", "-t Resolution Distribution (%); (t_{rec_pT} - t_{truth})/t_{truth} (%)", 100, -110, 210);
  htw_res3->SetLineWidth(2);
  TH1* htw_res4 = new TH1D("htw_res4", "-t Resolution Distribution (%); (t_{rec_corr} - t_{truth})/t_{truth} (%)", 100, -100, 100);
  htw_res4->SetLineWidth(2);
  TH1* htw_res5 = new TH1D("htw_res5", "Q^{2} Resolution Distribution (%); (Q^{2}_{rec} - Q^{2}_{truth})/Q^{2}_{truth} (%)", 100, -20, 20);
  htw_res5->SetLineWidth(2);
  TH1* htw_res6 = new TH1D("htw_res6", "W Resolution Distribution (%); (W_{rec} - W_{truth})/W_{truth} (%)", 100, -60, 60);
  htw_res6->SetLineWidth(2);

  // Effeciency plots
  TH2* Q2_t_DetEff_Uncut = new TH2F("Q2_t_DetEff_Uncut", "Q^{2}_{truth} vs -t_{truth} for thrown events; Q^{2} (GeV^{2}); -t (GeV^{2})", 10, 0, 40, 10, 0, 1.5);
  TH2* Q2_t_DetEff_Cut = new TH2F("Q2_t_DetEff_Cut", "Q^{2}_{truth} vs -t_{truth} for detected events; Q^{2} (GeV^{2}); -t (GeV^{2})", 10, 0, 40, 10, 0, 1.5);
  TH2* Q2_t_DetEff = new TH2F("Q2_t_DetEff", "Q^{2}_{truth} vs -t_{truth} detected/thrown ratio; Q^{2} (GeV^{2}); -t (GeV^{2})", 10, 0, 40, 10, 0, 1.5);

  TH1* eTruthw_Eta_Uncut = new TH1D("eTruthw_Eta_Uncut", "e' #eta for thrown events; #eta",100,-2,0);
  TH1* eRecw_Eta_Cut = new TH1D("eRecw_Eta_Cut", "e' #eta for detected events; #eta",100,-2,0);
  TH1* eEff_Eta = new TH1D("eEff_Eta", "e' Tracking efficiency as fn of #eta; #eta; Eff(%)", 100,-2,0);

  TH1* piTruthw_Eta_Uncut = new TH1D("piTruthw_Eta_Uncut", "#pi^{+} #eta for thrown events; #eta",100,0,4);
  TH1* piRecw_Eta_Cut = new TH1D("piRecw_Eta_Cut", "#pi^{+} #eta for detected events; #eta",100,0,4);
  TH1* piEff_Eta = new TH1D("piEff_Eta", "#pi^{+} Tracking efficiency as fn of #eta; #eta; Eff(%)",100,0,4); 

  TH1* eTruthw_P_Uncut = new TH1D("eTruthw_P_Uncut", "e' P for thrown events; P (GeV/c)",100,4,7);
  TH1* eRecw_P_Cut = new TH1D("eRecw_P_Cut", "e' P for detected events; P (GeV/c)",100,4,7);
  TH1* eEff_P = new TH1D("eEff_P", "e' Tracking efficiency as fn of P; P (GeV/c); Eff(%)", 100,4,7);

  TH1* piTruthw_P_Uncut = new TH1D("piTruthw_P_Uncut", "#pi^{+} P for thrown events; P (GeV/c)",100,0,30);
  TH1* piRecw_P_Cut = new TH1D("piRecw_P_Cut", "#pi^{+} P for detected events; P (GeV/c)",100,0,30);
  TH1* piEff_P = new TH1D("piEff_P", "#pi^{+} Tracking efficiency as fn of P; P (GeV/c); Eff(%)",100,0,30); 

  // Neutrons HCal
  TH2* nRecw_Thetap_hcal  = new TH2D("nRecw_Thethap_hcal","n rec #theta vs P for 1 cluster events; #theta (Deg); P (GeV/c)",100,0.8,2.0,100,0,40);
  TH2* nRecw_rot_Thetap_hcal  = new TH2D("nRecw_rot_Thethap_hcal","n rec #theta* vs P around p axis for 1 cluster events; #theta* (mRad); P (GeV/c)",100,0,4.0,100,5,40);

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //Defining the four vectors
  ROOT::Math::PxPyPzEVector elec_beam; // initialized the 4 vector for electron beam
  ROOT::Math::PxPyPzEVector prot_beam; // initialized the 4 vector for proton beam

  ROOT::Math::PxPyPzEVector elec_mc; // initialized the 4 vector for truth electron
  ROOT::Math::PxPyPzEVector pi_mc; // initialized the 4 vector for truth pion
  ROOT::Math::PxPyPzEVector neut_mc; // initialized the 4 vector for truth neutron
  ROOT::Math::PxPyPzEVector neut_rot_mc; // initialized the 4 vector for truth neutron with rotation of 25 mrad

  ROOT::Math::RotationY rot; // initialized rotation vector
  rot.SetAngle(0.025);

  ROOT::Math::PxPyPzEVector elec_rec; // initialized the 4 vector for reconstructed electron
  ROOT::Math::PxPyPzEVector pi_rec; // initialized the 4 vector for reconstructed pion
  ROOT::Math::PxPyPzEVector neut_rec; // initialized the 4 vector for reconstructed neutron
  ROOT::Math::PxPyPzEVector neut_rot_rec; // initialized the 4 vector for reconstructed neutron with rotation of 25 mrad

  ROOT::Math::PxPyPzEVector virtphoton_truth; // intialized the 4 vector for truth virtual photon
  ROOT::Math::PxPyPzEVector ttruth; // intialized the 4 vector for ttruth (-t)from first loop
  ROOT::Math::PxPyPzEVector talttruth; // intialized the 4 vector for talttruth(-t) from second loop

  ROOT::Math::PxPyPzEVector virtphoton_rec; //intialized the 4 vector for reconstructed virtual photon
  ROOT::Math::PxPyPzEVector trec; // intialized the 4 vector for trec (-t)from first loop
  ROOT::Math::PxPyPzEVector taltrec; // intialized the 4 vector for taltrec(-t) from second loop
  ROOT::Math::PxPyPzEVector trecpT; // intialized the 4 vector for trecpT(-t)
  ROOT::Math::PxPyPzEVector p_miss_rec; //intialized the 4 vector for missing momentum
  ROOT::Math::PxPyPzMVector neut_corr, neut_corr_rot; // intialized the 4 vector for reconstructed corrected neutron
  ROOT::Math::PxPyPzEVector treccorr; // intialized the 4 vector for trecpT(-t)

  ROOT::Math::XYZVector neut_pos_hcal; // initialized the 3 vector for zdc position
  ROOT::Math::PxPyPzEVector neut_rec_hcal; // initialized the 4 vector for reconstructed neutorn in hcal
  ROOT::Math::PxPyPzEVector neut_rot_rec_hcal; // initialized the 4 vector for reconstructed neutron with rotation of 25 mrad in hcal
  int hcal_clus_size;
  double neut_rec_p_hcal;

  double neutMass = 0.93965420;
  double weight, partEng; // weight and energy of the particles
  double Q2_truth, W_truth, y_truth, t_truth, t_alttruth; // Truth kinematic variables
  double Q2_rec, W_rec, y_rec, t_rec, t_altrec, t_recpT, t_reccorr; // Reconstructed kinematic variables

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Defining initial colliding beams
  double eMass = 0.000510998950; //electron beam
  double eEng = 5;
  double e_pmag = sqrt(pow(eEng,2)-pow(eMass,2));
  double e_p1 = 0.,e_p2 = 0.;
  double e_p3 = -1*e_pmag;
  elec_beam.SetPxPyPzE(e_p1, e_p2, e_p3, eEng); 
 
              
  double pMass = 0.93827208816; // proton beam
  double pEng = 41;
  double p_pmag = sqrt(pow(pEng,2)-pow(pMass,2));
  double c_a = 0.025;
  double p_p1 = -p_pmag*TMath::Sin(c_a);
  double p_p2 = 0.;
  double p_p3 = p_pmag*TMath::Cos(c_a);
  prot_beam.SetPxPyPzE(p_p1, p_p2, p_p3, pEng);

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  bool x,y,z; // x,y, and z are for reconstructed electron, pion, and neutron
 
  while(tree_reader.Next()) { // Loop over events
 
    x = false, y = false, z = false;
 
    std::string weight_name = weight_map[0].first; // accessing weights of particles
    vector<string> weight_value = weight_map[0].second;
    weight = std::stod( *(weight_value.begin()) );

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
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
 
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    for(unsigned int i=0; i<trackPdg.GetSize(); i++) { // Loop over reconstructed particles 
      if(trackPdg[i] == 11) { // Look at electron
	x = true;
	elec_rec.SetPxPyPzE(trackMomX[i],trackMomY[i],trackMomZ[i], trackEng[i]);
	eRecw_Thetap -> Fill(elec_rec.Theta()*TMath::RadToDeg(), elec_rec.P(), weight);
	eRecw_Eta_Cut -> Fill(elec_mc.Eta(), weight);
	eRecw_P_Cut -> Fill(elec_mc.P(), weight);
      }
 
      if(trackPdg[i] == 211) { // Look at pion
	y = true;
	pi_rec.SetPxPyPzE(trackMomX[i],trackMomY[i],trackMomZ[i], trackEng[i]);
	piRecw_Thetap -> Fill(pi_rec.Theta()*TMath::RadToDeg(), pi_rec.P(), weight);
	piRecw_Eta_Cut -> Fill(pi_mc.Eta(), weight);
	piRecw_P_Cut -> Fill(pi_mc.P(), weight);
      }
 
    }// for over reconstructed particles

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    for(unsigned int i=0; i<neutEng.GetSize(); i++) { // Loop over zdc neutrons
 
      neut_rec.SetPxPyPzE(neutMomX[i],neutMomY[i],neutMomZ[i], neutEng[i]);
      nRecw_Thetap -> Fill(neut_rec.Theta()*TMath::RadToDeg(), neut_rec.P(), weight);
      nRecw_Thetaphi -> Fill(neut_rec.Theta()*1000., neut_rec.Phi()*TMath::RadToDeg(), weight);
 
      neut_rot_rec = rot*neut_rec; // rotate w.r.t to proton axis
      nRecw_rot_Thetaphi -> Fill(neut_rot_rec.Theta()*1000., neut_rot_rec.Phi()*TMath::RadToDeg(), weight);
 
      if(neut_rot_rec.Theta()*1000. < 4.0){ // acceptance of the zdc
	nRec_en -> Fill(neut_rot_rec.E());
	nRec_clus -> Fill(neutClus[i]);
 
	if(neut_rot_rec.E()>10.0){ // neutron energy cut
	  z = true;
	  nRecw_rot_Thetap -> Fill(neut_rot_rec.Theta()*1000., neut_rot_rec.P(), weight);
	}
 
      }
 
    }// for over zdc neutrons
 
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
 
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Reconstructed kinematic variables

    if (x == true && y == true && z == true ){ // if e', pi, and neutron are in coincidence
 
      virtphoton_rec = (elec_beam - elec_rec);
      Q2_rec = -1*(virtphoton_rec.mag2()); 
      W_rec = ((virtphoton_rec + prot_beam).mag()); 
      y_rec =  (prot_beam.Dot(virtphoton_rec))/(prot_beam.Dot(elec_beam)); // Energy Loss y

      if(Q2_rec > 5 && Q2_rec < 35){ // Q2 Cut 
 
	// t-method plots
	trec = (virtphoton_rec - pi_rec); // First method to reconstruct -t
	t_rec = -1*(trec.mag2()); // t_rec is the -t from the first loop
	htw_rec1 -> Fill(t_rec, t_truth, weight); 
	htwz_rec1 -> Fill(t_rec, t_truth, weight); 
                                 
	taltrec = (prot_beam - neut_rec); // Second method to reconstruct -t
	t_altrec = -1*(taltrec.mag2()); // t_altrec is the -t from the second loop
	htw_rec2 -> Fill(t_altrec, t_truth, weight); 
	htwz_rec2 -> Fill(t_altrec, t_truth, weight);
 
	trecpT = (elec_rec +  pi_rec); // Third method to reconstruct -t
	t_recpT = trecpT.Perp2();
	htw_rec3 -> Fill(t_recpT, t_truth, weight); 
	htwz_rec3 -> Fill(t_recpT, t_truth, weight);
 
	p_miss_rec = (elec_beam + prot_beam) - (elec_rec + pi_rec) ; // Defined missing momentum information -> Fourth method to reconstruct -t
 
	neut_corr.SetCoordinates( p_miss_rec.P() * sin(neut_rec.Theta()) * cos(neut_rec.Phi()), 
				  p_miss_rec.P() * sin(neut_rec.Theta()) * sin(neut_rec.Phi()), 
				  p_miss_rec.P() * cos(neut_rec.Theta()),
				  neutMass );
 
	treccorr = (prot_beam - neut_corr);
	t_reccorr = -1*(treccorr.mag2());
	htw_rec4 -> Fill(t_reccorr, t_truth, weight); 
	htwz_rec4 -> Fill(t_reccorr, t_truth, weight); 
 
	// Absolute difference -t plots
	htw_t1 -> Fill(t_rec - t_truth);
	htw_t2 -> Fill(t_altrec - t_truth);
	htw_t3 -> Fill(t_recpT - t_truth);
	htw_t4 -> Fill(t_reccorr - t_truth);

	// Neutron theta-phi plots 
	n_ThetaDiff -> Fill((p_miss_rec.Theta() - neut_rec.Theta())*TMath::RadToDeg(), weight);
	n_PhiDiff -> Fill((p_miss_rec.Phi() - neut_rec.Phi())*TMath::RadToDeg(), weight);
	n_ThetaPhiDiff -> Fill((p_miss_rec.Theta() - neut_rec.Theta())*TMath::RadToDeg(), (p_miss_rec.Phi() - neut_rec.Phi())*TMath::RadToDeg(), weight);
 
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
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  gStyle->SetPadRightMargin(0.125); // left space on right side
  gStyle->SetStatX(0.895); // move the stat bax on left or right side
  gStyle->SetStatY(0.90); // move the stat bax on up or down side
  gStyle->SetPalette(55);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetHistLineWidth(2);
  
  TCanvas *c1 = new TCanvas("c1"); //truth information     
  c1->SetLogz(); 
  eTruthw_Thetap->Scale(1.0/nFiles);
  eTruthw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  eTruthw_Thetap->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogz();
  piTruthw_Thetap->Scale(1.0/nFiles);
  piTruthw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  piTruthw_Thetap->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogz();
  nTruthw_Thetaphi->Scale(1.0/nFiles);
  nTruthw_Thetaphi->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nTruthw_Thetaphi->Draw("colz");
  
  TCanvas *c3a = new TCanvas("c3a");
  c3a->SetLogz();
  nTruthw_rot_Thetaphi->Scale(1.0/nFiles);
  nTruthw_rot_Thetaphi->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nTruthw_rot_Thetaphi->Draw("colz");
  
  TCanvas *c3b = new TCanvas("c3b"); 
  c3b->SetLogz();
  nTruthw_Thetap->Scale(1.0/nFiles);
  nTruthw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nTruthw_Thetap->Draw("colz");  
  
  TCanvas *c3c = new TCanvas("c3c"); 
  c3c->SetLogz();
  nTruthw_rot_Thetap->Scale(1.0/nFiles);
  nTruthw_rot_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nTruthw_rot_Thetap->Draw("colz");  
  
  TCanvas *c4 = new TCanvas("c4"); //reconstructed information
  c4->SetLogz();
  eRecw_Thetap->Scale(1.0/nFiles);
  eRecw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  eRecw_Thetap->Draw("colz");
  
  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogz();
  piRecw_Thetap->Scale(1.0/nFiles);
  piRecw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  piRecw_Thetap->Draw("colz");
  
  TCanvas *c5a = new TCanvas("c5a");
  c5a->SetLogz();
  nRecw_Thetaphi->Scale(1.0/nFiles);
  nRecw_Thetaphi->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_Thetaphi->Draw("colz");
  
  TCanvas *c5b = new TCanvas("c5b");
  c5b->SetLogz();
  nRecw_rot_Thetaphi->Scale(1.0/nFiles);
  nRecw_rot_Thetaphi->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_rot_Thetaphi->Draw("colz");
  
  TCanvas *c5c = new TCanvas("c5c");
  c5c->SetLogz();
  nRecw_Thetap_hcal->Scale(1.0/nFiles);
  nRecw_Thetap_hcal->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_Thetap_hcal->Draw("colz");
  
  TCanvas *c5d = new TCanvas("c5d");
  c5d->SetLogz();
  nRecw_rot_Thetap_hcal->Scale(1.0/nFiles);
  nRecw_rot_Thetap_hcal->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_rot_Thetap_hcal->Draw("colz");
  
  TCanvas *c6a = new TCanvas("c6a");
  nRec_en->Draw();
  
  TCanvas *c6b = new TCanvas("c6b");
  nRec_clus->Draw();
  
  TCanvas *c7 = new TCanvas("c7");
  c7->SetLogz();
  nRecw_Thetap->Scale(1.0/nFiles);
  nRecw_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_Thetap->Draw("colz");
  
  TCanvas *c7a = new TCanvas("c7a");
  c7a->SetLogz();
  nRecw_rot_Thetap->Scale(1.0/nFiles);
  nRecw_rot_Thetap->GetZaxis()->SetTitle("Rate/bin (Hz)");
  nRecw_rot_Thetap->Draw("colz");
  
  TCanvas *c8 = new TCanvas("c8"); // -t reconstruction plots
  c8->SetLogz();
  htw_rec1->Scale(1.0/nFiles);
  htw_rec1->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htw_rec1->Draw("colz");
  
  TCanvas *c8a = new TCanvas("c8a");
  c8a->SetLogz();
  htwz_rec1->Scale(1.0/nFiles);
  htwz_rec1->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htwz_rec1->Draw("colz");
  
  TCanvas *c9 = new TCanvas("c9");
  c9->SetLogz();
  htw_rec2->Scale(1.0/nFiles);
  htw_rec2->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htw_rec2->Draw("colz");
  
  TCanvas *c9a = new TCanvas("c9a");
  c9a->SetLogz();
  htwz_rec2->Scale(1.0/nFiles);
  htwz_rec2->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htwz_rec2->Draw("colz");
  
  TCanvas *c10 = new TCanvas("c10");
  c10->SetLogz();
  htw_rec3->Scale(1.0/nFiles);
  htw_rec3->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htw_rec3->Draw("colz");
  
  TCanvas *c10a = new TCanvas("c10a");
  c10a->SetLogz();
  htwz_rec3->Scale(1.0/nFiles);
  htwz_rec3->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htwz_rec3->Draw("colz");
  
  TCanvas *c11 = new TCanvas("c11");
  c11->SetLogz();
  htw_rec4->Scale(1.0/nFiles);
  htw_rec4->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htw_rec4->Draw("colz");
  
  TCanvas *c11a = new TCanvas("c11a");
  c11a->SetLogz();
  htwz_rec4->Scale(1.0/nFiles);
  htwz_rec4->GetZaxis()->SetTitle("Rate/bin (Hz)");
  htwz_rec4->Draw("colz");
  
  TCanvas *c12 = new TCanvas("c12"); // Resolution plots
  htw_res_e->Scale(1.0/nFiles);
  htw_res_e->Draw("HIST");
  
  TCanvas *c13 = new TCanvas("c13");
  htw_res_pi->Scale(1.0/nFiles);
  htw_res_pi->Draw("HIST");
  
  TCanvas *c14a = new TCanvas("c14a");
  htw_res_n1->Scale(1.0/nFiles);
  htw_res_n1->Draw("HIST");
  
  TCanvas *c14b = new TCanvas("c14b");
  htw_res_n2->Scale(1.0/nFiles);
  htw_res_n2->Draw("HIST");
  
  TCanvas *c14c = new TCanvas("c14c");
  htw_res_n3->Scale(1.0/nFiles);
  htw_res_n3->Draw("HIST");
  
  TCanvas *c14d = new TCanvas("c14d");
  htw_res_n4->Scale(1.0/nFiles);
  htw_res_n4->GetYaxis()->SetTitle("Rate (Hz)");
  htw_res_n4->Draw("HIST");
  
  TCanvas *c15a = new TCanvas("c15a"); //t-resoutions
  htw_res1->Scale(1.0/nFiles);
  htw_res1->Draw("HIST");
  
  TCanvas *c15b = new TCanvas("c15b");
  htw_res2->Scale(1.0/nFiles);
  htw_res2->Draw("HIST");
  
  TCanvas *c15c = new TCanvas("c15c");
  htw_res3->Scale(1.0/nFiles);
  htw_res3->Draw("HIST");
  
  TCanvas *c15d = new TCanvas("c15d");
  htw_res4->Scale(1.0/nFiles);
  htw_res4->Draw("HIST");
  
  TCanvas *c15e = new TCanvas("c15e");
  htw_res5->Scale(1.0/nFiles);
  htw_res5->Draw("HIST");
  
  TCanvas *c15f = new TCanvas("c15f");
  htw_res6->Scale(1.0/nFiles);
  htw_res6->Draw("HIST");
  
  TCanvas *c16a = new TCanvas("c16a"); // Neutron theta-phi plots
  n_ThetaDiff->Scale(1.0/nFiles);
  n_ThetaDiff->Draw("HIST");
  
  TCanvas *c16b = new TCanvas("c16b");
  n_PhiDiff->Scale(1.0/nFiles);
  n_PhiDiff->Draw("HIST");
  
  TCanvas *c16c = new TCanvas("c16c");
  c16c->SetLogz();
  n_ThetaPhiDiff->Scale(1.0/nFiles);
  n_ThetaPhiDiff->GetZaxis()->SetTitle("Rate/bin (Hz)");
  n_ThetaPhiDiff->Draw("colz");
  
  TCanvas *c17 = new TCanvas("c17"); // Absolute difference -t plots
  htw_t1->Scale(1.0/nFiles), htw_t2->Scale(1.0/nFiles), htw_t3->Scale(1.0/nFiles), htw_t4->Scale(1.0/nFiles);
  htw_t4->Draw("HIST");
  htw_t2->Draw("HIST SAME");
  htw_t3->Draw("HIST SAME");
  htw_t1->Draw("HIST SAME");
  
  TLegend *leg17 = new TLegend (0.8,0.45,0.6,0.75);
  leg17->SetBorderSize(0);leg17->SetFillStyle(0); 
  leg17->AddEntry(htw_t1,"t_{rec} - t_{truth}","l");
  leg17->AddEntry(htw_t2,"t_{alt_rec} - t_{truth}","l");
  leg17->AddEntry(htw_t3,"t_{recpT} - t_{truth}","l");
  leg17->AddEntry(htw_t4,"t_{rec_corr} - t_{truth}","l");
  leg17->Draw();
  
  TCanvas *c18 = new TCanvas("c18"); // Efficiency plots
  Q2_t_DetEff->Draw("colz");
  
  TCanvas *c18a = new TCanvas("c18a"); 
  Q2_t_DetEff_Cut->Scale(1.0/nFiles);
  Q2_t_DetEff_Cut->GetZaxis()->SetTitle("Rate/bin (Hz)");
  Q2_t_DetEff_Cut->Draw("colz");
  
  TCanvas *c18b = new TCanvas("c18b"); 
  Q2_t_DetEff_Uncut->Scale(1.0/nFiles);
  Q2_t_DetEff_Cut->GetZaxis()->SetTitle("Rate/bin (Hz)");
  Q2_t_DetEff_Uncut->Draw("colz"); 
 
  TCanvas *c19 = new TCanvas("c19"); // elec eta eff.
  eEff_Eta->Draw("HIST");
  
  TCanvas *c19a = new TCanvas("c19a"); 
  eTruthw_Eta_Uncut->Scale(1.0/nFiles);
  eTruthw_Eta_Uncut->Draw("HIST");
  
  TCanvas *c19b = new TCanvas("c19b"); 
  eRecw_Eta_Cut->Scale(1.0/nFiles);
  eRecw_Eta_Cut->Draw("HIST");
  
  TCanvas *c20 = new TCanvas("c20");  // elec mom eff.
  eEff_P->Draw("HIST");
  
  TCanvas *c20a = new TCanvas("c20a"); 
  eTruthw_P_Uncut->Scale(1.0/nFiles);
  eTruthw_P_Uncut->Draw("HIST");
  
  TCanvas *c20b = new TCanvas("c20b"); 
  eRecw_P_Cut->Scale(1.0/nFiles);
  eRecw_P_Cut->Draw("HIST");
  
  TCanvas *c21 = new TCanvas("c21"); // pi eta eff.
  piEff_Eta->Draw("HIST");
  
  TCanvas *c21a = new TCanvas("c21a"); 
  piTruthw_Eta_Uncut->Scale(1.0/nFiles);
  piTruthw_Eta_Uncut->Draw("HIST");
  
  TCanvas *c21b = new TCanvas("c21b"); 
  piRecw_Eta_Cut->Scale(1.0/nFiles);
  piRecw_Eta_Cut->Draw("HIST");
  
  TCanvas *c22 = new TCanvas("c22"); //pi mom eff.
  piEff_P->Draw("HIST");
  
  TCanvas *c22a = new TCanvas("c22a"); 
  piTruthw_P_Uncut->Scale(1.0/nFiles);
  piTruthw_P_Uncut->Draw("HIST");
  
  TCanvas *c22b = new TCanvas("c22b"); 
  piRecw_P_Cut->Scale(1.0/nFiles);
  piRecw_P_Cut->Draw("HIST");
  
  
  c1->Print("DEMP_reco.pdf[");
  c1->Print("DEMP_reco.pdf");
  c2->Print("DEMP_reco.pdf");
  c3->Print("DEMP_reco.pdf");
  c3a->Print("DEMP_reco.pdf");
  c3b->Print("DEMP_reco.pdf");
  c3c->Print("DEMP_reco.pdf");
  c4->Print("DEMP_reco.pdf");
  c5->Print("DEMP_reco.pdf");
  c5a->Print("DEMP_reco.pdf");
  c5b->Print("DEMP_reco.pdf");
  c5c->Print("DEMP_reco.pdf");
  c5d->Print("DEMP_reco.pdf");
  c6a->Print("DEMP_reco.pdf");
  c6b->Print("DEMP_reco.pdf");
  c7->Print("DEMP_reco.pdf");
  c7a->Print("DEMP_reco.pdf");
  c8->Print("DEMP_reco.pdf");
  c8a->Print("DEMP_reco.pdf");
  c9->Print("DEMP_reco.pdf");
  c9a->Print("DEMP_reco.pdf");
  c10->Print("DEMP_reco.pdf");
  c10a->Print("DEMP_reco.pdf");
  c11->Print("DEMP_reco.pdf");
  c11a->Print("DEMP_reco.pdf");
  c12->Print("DEMP_reco.pdf");
  c13->Print("DEMP_reco.pdf");
  c14a->Print("DEMP_reco.pdf");
  c14b->Print("DEMP_reco.pdf");
  c14c->Print("DEMP_reco.pdf");
  c14d->Print("DEMP_reco.pdf");
  c15a->Print("DEMP_reco.pdf");
  c15b->Print("DEMP_reco.pdf");
  c15c->Print("DEMP_reco.pdf");
  c15d->Print("DEMP_reco.pdf");
  c15e->Print("DEMP_reco.pdf");
  c15f->Print("DEMP_reco.pdf");
  c16a->Print("DEMP_reco.pdf");
  c16b->Print("DEMP_reco.pdf");
  c16c->Print("DEMP_reco.pdf");
  c17->Print("DEMP_reco.pdf");
  c18b->Print("DEMP_reco.pdf");
  c18a->Print("DEMP_reco.pdf");
  c18->Print("DEMP_reco.pdf");
  c19a->Print("DEMP_reco.pdf");
  c19b->Print("DEMP_reco.pdf");
  c19->Print("DEMP_reco.pdf");
  c20a->Print("DEMP_reco.pdf");
  c20b->Print("DEMP_reco.pdf");
  c20->Print("DEMP_reco.pdf");
  c21a->Print("DEMP_reco.pdf");
  c21b->Print("DEMP_reco.pdf");
  c21->Print("DEMP_reco.pdf");
  c22a->Print("DEMP_reco.pdf");
  c22b->Print("DEMP_reco.pdf");
  c22->Print("DEMP_reco.pdf");
  c22->Print("DEMP_reco.pdf]");
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
   ofile->Write(); // Write histograms to file
   ofile->Close(); // Close output file
    
}

