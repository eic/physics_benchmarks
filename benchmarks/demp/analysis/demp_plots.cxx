#include <fstream>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

void demp_plots(const std::string& config_name)
{

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string hists_file    = config["hists_file"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const int         ebeam         = config["ebeam"];
  const int         pbeam         = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DEMP analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file for histograms: {}\n", hists_file);
  fmt::print(" - output prefix for plots: {}\n", output_prefix);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Read file with histograms
  TFile* file = new TFile(hists_file.c_str());

  std::cout<<"Reading histograms..."<<std::endl;

  TH2D* eTruthw_Thetap = (TH2D*) file->Get("eTruthw_Thetap");
  TH2D* piTruthw_Thetap = (TH2D*) file->Get("piTruthw_Thetap");
  TH2D* nTruthw_Thetaphi = (TH2D*) file->Get("nTruthw_Thetaphi");
  TH2D* nTruthw_rot_Thetaphi = (TH2D*) file->Get("nTruthw_rot_Thetaphi");
  TH2D* nTruthw_Thetap = (TH2D*) file->Get("nTruthw_Thetap");
  TH2D* nTruthw_rot_Thetap = (TH2D*) file->Get("nTruthw_rot_Thetap");
  TH2D* eRecw_Thetap = (TH2D*) file->Get("eRecw_Thetap");
  TH2D* eRecw_Thetaphi = (TH2D*) file->Get("eRecw_Thetaphi");
  TH2D* piRecw_Thetap = (TH2D*) file->Get("piRecw_Thetap");
  TH2D* piRecw_Thetaphi = (TH2D*) file->Get("piRecw_Thetaphi");
  TH2D* nRecw_Thetaphi = (TH2D*) file->Get("nRecw_Thetaphi");
  TH2D* nRecw_rot_Thetaphi = (TH2D*) file->Get("nRecw_rot_Thetaphi");
  TH2D* nRecw_Thetap_hcal = (TH2D*) file->Get("nRecw_Thetap_hcal");
  TH2D* nRecw_rot_Thetap_hcal = (TH2D*) file->Get("nRecw_rot_Thetap_hcal");
  TH1D* nRec_en = (TH1D*) file->Get("nRec_en");
  TH1D* nRec_clus = (TH1D*) file->Get("nRec_clus");
  TH2D* nRecw_Thetap = (TH2D*) file->Get("nRecw_Thetap");
  TH2D* nRecw_rot_Thetap = (TH2D*) file->Get("nRecw_rot_Thetap");
  TH2D* htw_rec1 = (TH2D*) file->Get("htw_rec1");
  TH2D* htwz_rec1 = (TH2D*) file->Get("htwz_rec1");
  TH2D* htw_rec2 = (TH2D*) file->Get("htw_rec2");
  TH2D* htwz_rec2 = (TH2D*) file->Get("htwz_rec2");
  TH2D* htw_rec3 = (TH2D*) file->Get("htw_rec3");
  TH2D* htwz_rec3 = (TH2D*) file->Get("htwz_rec3");
  TH2D* htw_rec4 = (TH2D*) file->Get("htw_rec4");
  TH2D* htwz_rec4 = (TH2D*) file->Get("htwz_rec4");
  TH1D* htw_res_e = (TH1D*) file->Get("htw_res_e");
  TH1D* htw_res_pi = (TH1D*) file->Get("htw_res_pi");
  TH1D* htw_res_n1 = (TH1D*) file->Get("htw_res_n1");
  TH1D* htw_res_n2 = (TH1D*) file->Get("htw_res_n2");  
  TH1D* htw_res_n3 = (TH1D*) file->Get("htw_res_n3");
  TH1D* htw_res_n4 = (TH1D*) file->Get("htw_res_n4");
  TH1D* htw_res1 = (TH1D*) file->Get("htw_res1");
  TH1D* htw_res2 = (TH1D*) file->Get("htw_res2");
  TH1D* htw_res3 = (TH1D*) file->Get("htw_res3");
  TH1D* htw_res4 = (TH1D*) file->Get("htw_res4");
  TH1D* htw_res5 = (TH1D*) file->Get("htw_res5");
  TH1D* htw_res6 = (TH1D*) file->Get("htw_res6");
  TH1D* n_ThetaDiff = (TH1D*) file->Get("n_ThetaDiff");
  TH1D* n_PhiDiff = (TH1D*) file->Get("n_PhiDiff");
  TH2D* n_ThetaPhiDiff = (TH2D*) file->Get("n_ThetaPhiDiff");
  TH2D* pMissRecw_Thetaphi = (TH2D*) file->Get("pMissRecw_Thetaphi");
  TH2D* pMissRecw_rot_Thetaphi = (TH2D*) file->Get("pMissRecw_rot_Thetaphi");  
  TH2D* n_TruthRecw_ThetaPhiDiff = (TH2D*) file->Get("n_TruthRecw_ThetaPhiDiff");
  TH1D* htw_t1 = (TH1D*) file->Get("htw_t1");
  TH1D* htw_t2 = (TH1D*) file->Get("htw_t2");
  TH1D* htw_t3 = (TH1D*) file->Get("htw_t3");
  TH1D* htw_t4 = (TH1D*) file->Get("htw_t4"); 
  TH2D* Q2_t_DetEff = (TH2D*) file->Get("Q2_t_DetEff");
  TH2D* Q2_t_DetEff_Cut = (TH2D*) file->Get("Q2_t_DetEff_Cut");
  TH2D* Q2_t_DetEff_Uncut = (TH2D*) file->Get("Q2_t_DetEff_Uncut");
  TH1D* eEff_Eta = (TH1D*) file->Get("eEff_Eta");
  TH1D* eTruthw_Eta_Uncut = (TH1D*) file->Get("eTruthw_Eta_Uncut");
  TH1D* eRecw_Eta_Cut = (TH1D*) file->Get("eRecw_Eta_Cut");
  TH1D* eEff_P = (TH1D*) file->Get("eEff_P");
  TH1D* eTruthw_P_Uncut = (TH1D*) file->Get("eTruthw_P_Uncut");
  TH1D* eRecw_P_Cut = (TH1D*) file->Get("eRecw_P_Cut");
  TH1D* piEff_Eta = (TH1D*) file->Get("piEff_Eta");
  TH1D* piTruthw_Eta_Uncut = (TH1D*) file->Get("piTruthw_Eta_Uncut");
  TH1D* piRecw_Eta_Cut = (TH1D*) file->Get("piRecw_Eta_Cut");
  TH1D* piEff_P = (TH1D*) file->Get("piEff_P");
  TH1D* piTruthw_P_Uncut = (TH1D*) file->Get("piTruthw_P_Uncut");
  TH1D* piRecw_P_Cut = (TH1D*) file->Get("piRecw_P_Cut");

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Make plots and save to PDF file
  std::cout<<"Making plots..."<<std::endl;

  gStyle->SetPadRightMargin(0.125); // left space on right side
  gStyle->SetStatX(0.895); // move the stat bax on left or right side
  gStyle->SetStatY(0.90); // move the stat bax on up or down side
  gStyle->SetPalette(55);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetHistLineWidth(2);
  
  TCanvas *c1 = new TCanvas("c1"); //truth information     
  c1->SetLogz(); 
  eTruthw_Thetap->Draw("colz");

  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogz();
  piTruthw_Thetap->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogz();
  nTruthw_Thetaphi->Draw("colz");
  
  TCanvas *c3a = new TCanvas("c3a");
  c3a->SetLogz();
  nTruthw_rot_Thetaphi->Draw("colz");
  
  TCanvas *c3b = new TCanvas("c3b"); 
  c3b->SetLogz();
  nTruthw_Thetap->Draw("colz");  
  
  TCanvas *c3c = new TCanvas("c3c"); 
  c3c->SetLogz();
  nTruthw_rot_Thetap->Draw("colz");  
  
  TCanvas *c4 = new TCanvas("c4"); //reconstructed information
  c4->SetLogz();
  eRecw_Thetap->Draw("colz");
  
  TCanvas *c4a = new TCanvas("c4a"); 
  c4a->SetLogz();
  eRecw_Thetaphi->Draw("colz");
  
  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogz();
  piRecw_Thetap->Draw("colz");
  
  TCanvas *c5_a = new TCanvas("c5_a");
  c5_a->SetLogz();
  piRecw_Thetaphi->Draw("colz");
  
  TCanvas *c5a = new TCanvas("c5a");
  c5a->SetLogz();
  nRecw_Thetaphi->Draw("colz");
  
  TCanvas *c5b = new TCanvas("c5b");
  c5b->SetLogz();
  nRecw_rot_Thetaphi->Draw("colz");
  
  TCanvas *c5c = new TCanvas("c5c");
  c5c->SetLogz();
  nRecw_Thetap_hcal->Draw("colz");
  
  TCanvas *c5d = new TCanvas("c5d");
  c5d->SetLogz();
  nRecw_rot_Thetap_hcal->Draw("colz");
  
  TCanvas *c6a = new TCanvas("c6a");
  nRec_en->Draw("HIST");
  
  TCanvas *c6b = new TCanvas("c6b");
  nRec_clus->Draw();
  
  TCanvas *c7 = new TCanvas("c7");
  c7->SetLogz();
  nRecw_Thetap->Draw("colz");
  
  TCanvas *c7a = new TCanvas("c7a");
  c7a->SetLogz();
  nRecw_rot_Thetap->Draw("colz");
  
  TCanvas *c8 = new TCanvas("c8"); // -t reconstruction plots
  c8->SetLogz();
  htw_rec1->Draw("colz");
  
  TCanvas *c8a = new TCanvas("c8a");
  c8a->SetLogz();
  htwz_rec1->Draw("colz");
  
  TCanvas *c9 = new TCanvas("c9");
  c9->SetLogz();
  htw_rec2->Draw("colz");
  
  TCanvas *c9a = new TCanvas("c9a");
  c9a->SetLogz();
  htwz_rec2->Draw("colz");
  
  TCanvas *c10 = new TCanvas("c10");
  c10->SetLogz();
  htw_rec3->Draw("colz");
  
  TCanvas *c10a = new TCanvas("c10a");
  c10a->SetLogz();
  htwz_rec3->Draw("colz");
  
  TCanvas *c11 = new TCanvas("c11");
  c11->SetLogz();
  htw_rec4->Draw("colz");
 
  TCanvas *c11a = new TCanvas("c11a");
  c11a->SetLogz();
  htwz_rec4->Draw("colz");
  
  TCanvas *c12 = new TCanvas("c12"); // Resolution plots
  htw_res_e->Draw("HIST");
  
  TCanvas *c13 = new TCanvas("c13");
  htw_res_pi->Draw("HIST");
  
  TCanvas *c14a = new TCanvas("c14a");
  htw_res_n1->Draw("HIST");
  
  TCanvas *c14b = new TCanvas("c14b");
  htw_res_n2->Draw("HIST");
  
  TCanvas *c14c = new TCanvas("c14c");
  htw_res_n3->Draw("HIST");
  
  TCanvas *c14d = new TCanvas("c14d");
  htw_res_n4->Draw("HIST");
  
  TCanvas *c15a = new TCanvas("c15a"); //t-resoutions
  htw_res1->Draw("HIST");
  
  TCanvas *c15b = new TCanvas("c15b");
  htw_res2->Draw("HIST");
  
  TCanvas *c15c = new TCanvas("c15c");
  htw_res3->Draw("HIST");
  
  TCanvas *c15d = new TCanvas("c15d");
  htw_res4->Draw("HIST");
  
  TCanvas *c15e = new TCanvas("c15e");
  htw_res5->Draw("HIST");
  
  TCanvas *c15f = new TCanvas("c15f");
  htw_res6->Draw("HIST");

  TCanvas *c16a = new TCanvas("c16a"); // Neutron theta-phi plots
  n_ThetaDiff->Draw("HIST");
  
  TCanvas *c16b = new TCanvas("c16b");
  n_PhiDiff->Draw("HIST");
  
  TCanvas *c16c = new TCanvas("c16c");
  c16c->SetLogz();
  n_ThetaPhiDiff->Draw("colz");
  
  TCanvas *c16d = new TCanvas("c16d");
  c16d->SetLogz();
  pMissRecw_Thetaphi->Draw("colz");
  
  TCanvas *c16e = new TCanvas("c16e");
  c16e->SetLogz();
  pMissRecw_rot_Thetaphi->Draw("colz");
  
  TCanvas *c16f = new TCanvas("c16f");
  c16f->SetLogz();
  n_TruthRecw_ThetaPhiDiff->Draw("colz");
  
  TCanvas *c17 = new TCanvas("c17"); // Absolute difference -t plots
  htw_t4->Draw("HIST");
  htw_t3->Draw("HIST SAME");
  htw_t2->Draw("HIST SAME");
  htw_t1->Draw("HIST SAME");
  
  TLegend *leg17 = new TLegend (0.8,0.45,0.6,0.75);
  leg17->SetBorderSize(0);leg17->SetFillStyle(0); 
  leg17->AddEntry(htw_t1,"t_{rec} - t_{truth}","l");
  leg17->AddEntry(htw_t2,"t_{alt_rec} - t_{truth}","l");
  leg17->AddEntry(htw_t3,"t_{recpT} - t_{truth}","l");
  leg17->AddEntry(htw_t4,"t_{rec_corr} - t_{truth}","l");
  leg17->Draw();
  
  TCanvas *c18 = new TCanvas("c18"); // Efficiency plots
  Q2_t_DetEff-> GetZaxis()->SetRangeUser(0.0,1.0);
  Q2_t_DetEff->Draw("colz");
  
  TCanvas *c18a = new TCanvas("c18a"); 
  Q2_t_DetEff_Cut->Draw("colz");
  
  TCanvas *c18b = new TCanvas("c18b"); 
  Q2_t_DetEff_Uncut->Draw("colz"); 
 
  TCanvas *c19 = new TCanvas("c19"); // elec eta eff.
  eEff_Eta->Draw("HIST");
  
  TCanvas *c19a = new TCanvas("c19a"); 
  eTruthw_Eta_Uncut->Draw("HIST");
  
  TCanvas *c19b = new TCanvas("c19b"); 
  eRecw_Eta_Cut->Draw("HIST");
  
  TCanvas *c20 = new TCanvas("c20");  // elec mom eff.
  eEff_P->Draw("HIST");
  
  TCanvas *c20a = new TCanvas("c20a"); 
  eTruthw_P_Uncut->Draw("HIST");
  
  TCanvas *c20b = new TCanvas("c20b"); 
  eRecw_P_Cut->Draw("HIST");
  
  TCanvas *c21 = new TCanvas("c21"); // pi eta eff.
  piEff_Eta->Draw("HIST");
  
  TCanvas *c21a = new TCanvas("c21a"); 
  piTruthw_Eta_Uncut->Draw("HIST");
  
  TCanvas *c21b = new TCanvas("c21b"); 
  piRecw_Eta_Cut->Draw("HIST");
  
  TCanvas *c22 = new TCanvas("c22"); //pi mom eff.
  piEff_P->Draw("HIST");
  
  TCanvas *c22a = new TCanvas("c22a"); 
  piTruthw_P_Uncut->Draw("HIST");
  
  TCanvas *c22b = new TCanvas("c22b"); 
  piRecw_P_Cut->Draw("HIST");

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  // Print plots to pdf file
  c1->Print(fmt::format("{}.pdf[", output_prefix).c_str());
  c1->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c2->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c3b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c3->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c3a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c3c->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c4->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c4a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5_a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5c->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c5d->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c6b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c6a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c7->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c7a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c8->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c8a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c9->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c9a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c10->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c10a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c11->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c11a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c12->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c13->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c14a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c14b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c14c->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c14d->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15c->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15d->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15e->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c15f->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16f->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16d->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16e->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c16c->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c17->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c18b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c18a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c18->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c19a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c19b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c19->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c20a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c20b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c20->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c21a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c21b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c21->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c22a->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c22b->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c22->Print(fmt::format("{}.pdf", output_prefix).c_str());
  c22->Print(fmt::format("{}.pdf]", output_prefix).c_str());
  
}

