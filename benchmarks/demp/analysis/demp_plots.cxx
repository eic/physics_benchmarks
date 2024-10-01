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

//--------------------------------------------------------------------------------------------------------------------------------------------

//Plotting style for histograms
void plots_2D(TH2D* hist, const int& nfiles = 1, const double font = 42, const double tsize = 0.04){
 
  hist->SetStats(1);
  
  hist->GetXaxis()->SetTitleFont(font); //x-axis
  hist->GetXaxis()->SetTitleSize(tsize);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetLabelFont(font);
  hist->GetXaxis()->SetLabelSize(tsize);
  hist->GetXaxis()->SetLabelOffset(0.01);
 
  hist->GetYaxis()->SetTitleFont(font); //y-axis
  hist->GetYaxis()->SetTitleSize(tsize);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetLabelFont(font);
  hist->GetYaxis()->SetLabelSize(tsize);
  hist->GetYaxis()->SetLabelOffset(0.01);

  hist->GetZaxis()->SetTitleFont(font); //z-axis
  hist->GetZaxis()->SetTitleSize(0.036);
  hist->GetZaxis()->SetTitleOffset(1.05);
  hist->GetZaxis()->SetLabelFont(font);
  hist->GetZaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetLabelOffset(0.01);
  
  hist->Scale(1.0/nfiles);
  hist->Draw("colz");

}

void plots_1D(TH1D* hist, const int& nfiles = 1, const double font = 42, const double tsize = 0.04){
   
  hist->SetStats(1);

  hist->GetXaxis()->SetTitleFont(font); //x-axis
  hist->GetXaxis()->SetTitleSize(tsize);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetLabelFont(font);
  hist->GetXaxis()->SetLabelSize(tsize);
  hist->GetXaxis()->SetLabelOffset(0.01);
 
  hist->GetYaxis()->SetTitleFont(font); //y-axis
  hist->GetYaxis()->SetTitleSize(tsize);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetLabelFont(font);
  hist->GetYaxis()->SetLabelSize(tsize);
  hist->GetYaxis()->SetLabelOffset(0.01);

  hist->Scale(1.0/nfiles);
  hist->Draw("HIST");

}

//--------------------------------------------------------------------------------------------------------------------------------------------

void demp_plots(const std::string& config_name)
{

  //--------------------------------------------------------------------------------------------------------------------------------------------
  
  // Read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string hists_file    = config["hists_file"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const int         ebeam         = config["ebeam"];
  const int         pbeam         = config["pbeam"];
  const int         nfiles        = config["nfiles"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DEMP analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file for histograms: {}\n", hists_file);
  fmt::print(" - output prefix for plots: {}\n", output_prefix);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);
  fmt::print(" - nfiles: {}\n", nfiles);

  //-------------------------------------------------------------------------------------------------------------------------------------------

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
  TH2D* nRecw_rot_PosXY = (TH2D*) file->Get("nRecw_rot_PosXY");
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

  //--------------------------------------------------------------------------------------------------------------------------------------------

  // Make plots and save to PDF file
  std::cout<<"Making plots..."<<std::endl;
  
  gStyle->SetPadRightMargin(0.125);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetOptTitle(1); //setting title
  gStyle->SetTitleAlign(1);
  gStyle->SetTitleX(0.35);
  gStyle->SetTitleY(0.92);
  gStyle->SetTitleW(0.32); 
  gStyle->SetTitleH(0.065);
  gStyle->SetOptStat(1); //setting stat box 
  gStyle->SetStatX(0.875);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.14);
  gStyle->SetStatH(0.11);
  gStyle->SetPalette(55);
  gStyle->SetHistLineWidth(2);
  
  TCanvas *c1 = new TCanvas("c1"); //truth information     
  c1->SetLogz();
  plots_2D(eTruthw_Thetap, nfiles);

  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogz();
  plots_2D(piTruthw_Thetap, nfiles);
   
  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogz();
  plots_2D(nTruthw_Thetaphi, nfiles);
 
  TCanvas *c3a = new TCanvas("c3a");
  c3a->SetLogz();
  plots_2D(nTruthw_rot_Thetaphi, nfiles);
  
  TCanvas *c3b = new TCanvas("c3b"); 
  c3b->SetLogz();
  plots_2D(nTruthw_Thetap, nfiles);   
  
  TCanvas *c3c = new TCanvas("c3c"); 
  c3c->SetLogz();
  plots_2D(nTruthw_rot_Thetap, nfiles);  
  
  TCanvas *c4 = new TCanvas("c4"); //reconstructed information
  c4->SetLogz();
  plots_2D(eRecw_Thetap, nfiles);
  
  TCanvas *c4a = new TCanvas("c4a"); 
  c4a->SetLogz();
  plots_2D(eRecw_Thetaphi, nfiles);
  
  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogz();
  plots_2D(piRecw_Thetap, nfiles);
  
  TCanvas *c5_a = new TCanvas("c5_a");
  c5_a->SetLogz();
  plots_2D(piRecw_Thetap, nfiles);
  
  TCanvas *c5a = new TCanvas("c5a");
  c5a->SetLogz();
  plots_2D(nRecw_Thetaphi, nfiles);
  
  TCanvas *c5b = new TCanvas("c5b");
  c5b->SetLogz();
  plots_2D(nRecw_rot_Thetaphi, nfiles);
  
  TCanvas *c5c = new TCanvas("c5c");
  c5c->SetLogz();
  plots_2D(nRecw_Thetap_hcal, nfiles);
  
  TCanvas *c5d = new TCanvas("c5d");
  c5d->SetLogz();
  plots_2D(nRecw_rot_Thetap_hcal, nfiles);
  
  TCanvas *c6a = new TCanvas("c6a");
  plots_1D(nRec_en, nfiles);
  
  TCanvas *c6b = new TCanvas("c6b");
  plots_1D(nRec_clus);
  
  TCanvas *c7 = new TCanvas("c7");
  c7->SetLogz();
  plots_2D(nRecw_Thetap, nfiles);
  
  TCanvas *c7a = new TCanvas("c7a");
  c7a->SetLogz();
  plots_2D(nRecw_rot_Thetap, nfiles);
 
  TCanvas *c7b = new TCanvas("c7b");
  c7b->SetLogz();
  plots_2D(nRecw_rot_PosXY, nfiles);

  TCanvas *c8 = new TCanvas("c8"); // -t reconstruction plots
  c8->SetLogz();
  plots_2D(htw_rec1, nfiles);
  
  TCanvas *c8a = new TCanvas("c8a");
  c8a->SetLogz();
  plots_2D(htwz_rec1, nfiles);
  
  TCanvas *c9 = new TCanvas("c9");
  c9->SetLogz();
  plots_2D(htw_rec2, nfiles);
  
  TCanvas *c9a = new TCanvas("c9a");
  c9a->SetLogz();
  plots_2D(htwz_rec2, nfiles);
  
  TCanvas *c10 = new TCanvas("c10");
  c10->SetLogz();
  plots_2D(htw_rec3, nfiles);
  
  TCanvas *c10a = new TCanvas("c10a");
  c10a->SetLogz();
  plots_2D(htwz_rec3, nfiles);
  
  TCanvas *c11 = new TCanvas("c11");
  c11->SetLogz();
  plots_2D(htw_rec4, nfiles);
 
  TCanvas *c11a = new TCanvas("c11a");
  c11a->SetLogz();
  plots_2D(htwz_rec4, nfiles);
  
  TCanvas *c12 = new TCanvas("c12"); // Resolution plots
  plots_1D(htw_res_e, nfiles);
  
  TCanvas *c13 = new TCanvas("c13");
  plots_1D(htw_res_pi, nfiles);
  
  TCanvas *c14a = new TCanvas("c14a");
  plots_1D(htw_res_n1, nfiles);
  
  TCanvas *c14b = new TCanvas("c14b");
  plots_1D(htw_res_n2, nfiles);
  
  TCanvas *c14c = new TCanvas("c14c");
  plots_1D(htw_res_n3, nfiles);
  
  TCanvas *c14d = new TCanvas("c14d");
  plots_1D(htw_res_n4, nfiles);
  
  TCanvas *c15a = new TCanvas("c15a"); //t-resoutions
  plots_1D(htw_res1, nfiles);
  
  TCanvas *c15b = new TCanvas("c15b");
  plots_1D(htw_res2, nfiles);
  
  TCanvas *c15c = new TCanvas("c15c");
  plots_1D(htw_res3, nfiles);
  
  TCanvas *c15d = new TCanvas("c15d");
  plots_1D(htw_res4, nfiles);
  
  TCanvas *c15e = new TCanvas("c15e");
  plots_1D(htw_res5, nfiles);
  
  TCanvas *c15f = new TCanvas("c15f");
  plots_1D(htw_res6, nfiles);

  TCanvas *c16a = new TCanvas("c16a"); // Neutron theta-phi plots
  plots_1D(n_ThetaDiff, nfiles);
  
  TCanvas *c16b = new TCanvas("c16b");
  plots_1D(n_PhiDiff, nfiles);
  
  TCanvas *c16c = new TCanvas("c16c");
  c16c->SetLogz();
  plots_2D(n_ThetaPhiDiff, nfiles);
  
  TCanvas *c16d = new TCanvas("c16d");
  c16d->SetLogz();
  plots_2D(pMissRecw_Thetaphi, nfiles);
  
  TCanvas *c16e = new TCanvas("c16e");
  c16e->SetLogz();
  plots_2D(pMissRecw_rot_Thetaphi, nfiles);
  
  TCanvas *c16f = new TCanvas("c16f");
  c16f->SetLogz();
  plots_2D(n_TruthRecw_ThetaPhiDiff, nfiles);
  
  TCanvas *c17 = new TCanvas("c17"); // Absolute difference -t plots
  htw_t1->Scale(1.0/nfiles), htw_t2->Scale(1.0/nfiles), htw_t3->Scale(1.0/nfiles), htw_t4->Scale(1.0/nfiles); // special case
  plots_1D(htw_t4);
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
  plots_2D(Q2_t_DetEff);
  
  TCanvas *c18a = new TCanvas("c18a"); 
  plots_2D(Q2_t_DetEff_Cut, nfiles);
  
  TCanvas *c18b = new TCanvas("c18b"); 
  plots_2D(Q2_t_DetEff_Uncut, nfiles); 
 
  TCanvas *c19 = new TCanvas("c19"); // elec eta eff.
  plots_1D(eEff_Eta);
  
  TCanvas *c19a = new TCanvas("c19a"); 
  plots_1D(eTruthw_Eta_Uncut, nfiles);
  
  TCanvas *c19b = new TCanvas("c19b"); 
  plots_1D(eRecw_Eta_Cut, nfiles);
  
  TCanvas *c20 = new TCanvas("c20");  // elec mom eff.
  plots_1D(eEff_P);
  
  TCanvas *c20a = new TCanvas("c20a"); 
  plots_1D(eTruthw_P_Uncut, nfiles);
  
  TCanvas *c20b = new TCanvas("c20b"); 
  plots_1D(eRecw_P_Cut, nfiles);
  
  TCanvas *c21 = new TCanvas("c21"); // pi eta eff.
  plots_1D(piEff_Eta);
  
  TCanvas *c21a = new TCanvas("c21a"); 
  plots_1D(piTruthw_Eta_Uncut, nfiles);
  
  TCanvas *c21b = new TCanvas("c21b"); 
  plots_1D(piRecw_Eta_Cut, nfiles);
  
  TCanvas *c22 = new TCanvas("c22"); //pi mom eff.
  plots_1D(piEff_P);
  
  TCanvas *c22a = new TCanvas("c22a"); 
  plots_1D(piTruthw_P_Uncut, nfiles);
  
  TCanvas *c22b = new TCanvas("c22b"); 
  plots_1D(piRecw_P_Cut, nfiles);

  //--------------------------------------------------------------------------------------------------------------------------------------------

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
  c7b->Print(fmt::format("{}.pdf", output_prefix).c_str());
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

