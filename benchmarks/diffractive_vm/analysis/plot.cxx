#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <TLorentzRotation.h>
#include <TVector2.h>
#include <TVector3.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

TH1D* makeHist(const char* name, const char* title, const char* xtit, const char* ytit,
               const int nBins, const double lower, const double higher, EColor color = kBlack)
{

  TH1D* temp = new TH1D(name, title, nBins, lower, higher);

  temp->SetMarkerSize(1.0);
  temp->SetMarkerStyle(20);
  temp->SetMarkerColor(color);
  temp->SetLineColor(color);
  temp->SetStats(kFALSE);

  temp->GetXaxis()->SetTitle(xtit);
  temp->GetXaxis()->SetTitleSize(0.05);
  temp->GetXaxis()->SetTitleFont(42);
  temp->GetXaxis()->SetTitleOffset(1.25);
  temp->GetXaxis()->SetLabelSize(0.05);
  temp->GetXaxis()->SetLabelOffset(0.01);
  temp->GetXaxis()->SetLabelFont(42);
  temp->GetXaxis()->SetLabelColor(kBlack);
  temp->GetXaxis()->CenterTitle();

  temp->GetYaxis()->SetTitle(ytit);
  temp->GetYaxis()->SetTitleSize(0.05);
  temp->GetYaxis()->SetTitleFont(42);
  temp->GetYaxis()->SetTitleOffset(1.4);
  temp->GetYaxis()->SetLabelSize(0.05);
  temp->GetYaxis()->SetLabelOffset(0.01);
  temp->GetYaxis()->SetLabelFont(42);
  temp->GetYaxis()->SetLabelColor(kBlack);
  temp->GetYaxis()->CenterTitle();

  return temp;
}

void fixedFontHist1D(TH1* h, Float_t xoffset = 1.5, Float_t yoffset = 2.3)
{
  h->SetLabelFont(43, "X");
  h->SetLabelFont(43, "Y");
  // h->SetLabelOffset(0.01);
  h->SetLabelSize(16);
  h->SetTitleFont(43);
  h->SetTitleSize(20);
  h->SetLabelSize(15, "Y");
  h->SetTitleFont(43, "Y");
  h->SetTitleSize(20, "Y");
  h->SetTitleOffset(xoffset, "X");
  h->SetTitleOffset(yoffset, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

int plot(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string vm_name       = config["vm_name"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const std::string test_tag      = config["test_tag"];
  const int         ebeam         = config["ebeam"];
  const int         pbeam         = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DIS electron analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix: {}\n", output_prefix);
  fmt::print(" - test tag: {}\n", test_tag);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);

  std::string output_name_dir = fmt::format("{}.root", output_prefix);
  cout << "Output file = " << output_name_dir << endl;
  TFile* file = new TFile(output_name_dir.c_str());

  TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
  TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
  TH1D* h_t_trk_REC = (TH1D*) file->Get("h_t_trk_REC");
  TH1D* h_t_combo_REC = (TH1D*) file->Get("h_t_combo_REC");
  TH2D* h_t_res = (TH2D*) file->Get("h_t_res");
  TH2D* h_trk_t_res = (TH2D*) file->Get("h_trk_t_res");
  TH2D* h_VM_res = (TH2D*) file->Get("h_VM_res");
  TH2D* h_energy_res = (TH2D*) file->Get("h_energy_res");
  TH1D* h_Q2_e = (TH1D*) file->Get("h_Q2_e");
  TH1D* h_y_e = (TH1D*) file->Get("h_y_e");  
  TH1D* h_energy_MC = (TH1D*) file->Get("h_energy_MC");
  TH1D* h_Q2REC_e = (TH1D*) file->Get("h_Q2REC_e");
  TH1D* h_yREC_e = (TH1D*) file->Get("h_yREC_e");
  TH1D* h_energy_REC = (TH1D*) file->Get("h_energy_REC");
  //Epz
  TH1D* h_Epz_REC = (TH1D*) file->Get("h_Epz_REC");
  TH1D* h_trk_Epz_REC = (TH1D*) file->Get("h_trk_Epz_REC");
  //cluster
  TH1D* h_EoverP_REC = (TH1D*) file->Get("h_EoverP_REC");  
  TH2D* h_emHits_position_REC = (TH2D*) file->Get("h_emHits_position_REC");//default cluster positio
  TH1D* h_energy_calibration_REC = (TH1D*) file->Get("h_energy_calibration_REC");
  TH1D* h_ClusOverHit_REC = (TH1D*) file->Get("h_ClusOverHit_REC");

  TString vm_label;
  TString daug_label;
  if (vm_name == "phi") {
    vm_label   = "#phi";
    daug_label = "K^{+}K^{-}";
  } else if (vm_name == "jpsi") {
    vm_label   = "J/#psi";
    daug_label = "e^{+}e^{-}";
  } else {
    throw std::runtime_error(fmt::format("Unknown vm_name = \"{}\"", vm_name));
  }

  {

    TCanvas c("c1", "c1", 1, 1, 600, 600);
    gPad->SetLogy(1);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0,
                           0.18, kBlack);
    base1->GetYaxis()->SetRangeUser(8e-2, 8e5);
    base1->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base1, 1., 1.2);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize() * 1.5);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize() * 1.5);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize() * 1.5);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize() * 1.5);
    base1->GetXaxis()->SetNdivisions(4, 4, 0);
    base1->GetYaxis()->SetNdivisions(5, 5, 0);
    base1->Draw();

    h_t_MC->Draw("same");

    h_t_REC->SetMarkerStyle(20);
    h_t_REC->Draw("PEsame");

    h_t_trk_REC->SetFillColorAlpha(kBlue, 0.4);
    h_t_trk_REC->SetFillStyle(1001);
    h_t_trk_REC->SetMarkerStyle(24);
    h_t_trk_REC->SetMarkerColor(kBlue);
    // h_t_trk_REC->Draw("PE3same");

    h_t_combo_REC->SetFillColorAlpha(kRed, 0.4);
    h_t_combo_REC->SetFillStyle(1001);
    h_t_combo_REC->SetMarkerStyle(24);
    h_t_combo_REC->SetMarkerColor(kRed);
    // h_t_combo_REC->Draw("PE3same");

    TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
    r42->SetNDC();
    r42->SetTextSize(22);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);
    r42->Draw("same");

    TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    r43->Draw("same");

    TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.95");
    r44->SetNDC();
    r44->SetTextSize(20);
    r44->SetTextFont(43);
    r44->SetTextColor(kBlack);
    r44->Draw("same");

    TLatex* r44_0 = new TLatex(0.18, 0.79,
                               "  |y_{" + vm_label + "}|<3.5, |M_{inv} #minus M_{" + vm_label +
                                   "}| < 0.02 GeV");
    r44_0->SetNDC();
    r44_0->SetTextSize(20);
    r44_0->SetTextFont(43);
    r44_0->SetTextColor(kBlack);
    r44_0->Draw("same");

    TLatex* r44_2 = new TLatex(0.18, 0.18, "" + vm_label + " #rightarrow " + daug_label);
    r44_2->SetNDC();
    r44_2->SetTextSize(30);
    r44_2->SetTextFont(43);
    r44_2->SetTextColor(kBlack);
    r44_2->Draw("same");

    TLegend* w7 = new TLegend(0.48, 0.68, 0.93, 0.76);
    w7->SetLineColor(kWhite);
    w7->SetFillColor(0);
    w7->SetTextSize(17);
    w7->SetTextFont(45);
    // if(filename=="MCclusterEnergy")
    //{
    //	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
    //	w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true EEMC E ", "P");
    //	w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only ", "P");

    //}
    // else if(filename=="MCvmAndelectron")
    //{
    //	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
    //	w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true e' ", "P");
    //	w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" MC w. RECO e' ", "P");
    //}
    // else{
    w7->AddEntry(h_t_MC, "Sartre " + vm_label + " MC ", "L");
    w7->AddEntry(h_t_REC, "Sartre " + vm_label + " RECO w. EEMC ", "P");
    // w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only", "P");
    // w7->AddEntry(h_t_combo_REC, "Sartre "+vm_label+" RECO best", "P");

    //}

    w7->Draw("same");

    c.Print(fmt::format("{}_benchmark-{}-dsigmadt.pdf", output_prefix, vm_name).c_str());
  }

  {
    TCanvas c("c2", "c2", 1, 1, 800, 800);
    gPad->SetTicks();
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.13);

    TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", " #delta t/t (resolution) ", 100, 0,
                           0.2, kBlack);
    base1->GetYaxis()->SetRangeUser(1e-2, 1000);
    base1->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base1, 1.2, 1.6);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize() * 1.5);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize() * 1.5);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize() * 1.8);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize() * 1.7);
    base1->GetXaxis()->SetNdivisions(5, 5, 0);
    base1->GetYaxis()->SetNdivisions(5, 5, 0);
    base1->Draw();
    TH2D* h_res;
    h_res = (TH2D*)h_t_res;
    // h_res=(TH2D*) h_trk_t_res;

    TH1D* h_res_1D = new TH1D("h_res_1D", "", 100, 0, 0.2);
    for (int ibin = 0; ibin < h_res->GetNbinsX(); ibin++) {
      TH1D*  tmp        = h_res->ProjectionY("tmp", ibin + 1, ibin + 1);
      double sigma      = tmp->GetStdDev();
      double sigmaerror = tmp->GetStdDevError();
      h_res_1D->SetBinContent(ibin + 1, sigma);
      h_res_1D->SetBinError(ibin + 1, sigmaerror);
    }

    h_res_1D->SetMarkerSize(1.6);
    h_res_1D->SetMarkerColor(kBlack);
    h_res_1D->SetLineColor(kBlack);
    h_res_1D->SetMarkerStyle(20);

    h_res_1D->Fit("pol0", "RMS0", "", 0.015, 0.025); // first  dip
    h_res_1D->Fit("pol0", "RMS0", "", 0.050, 0.060); // second dip
    h_res_1D->Fit("pol0", "RMS0", "", 0.100, 0.150); // third  dip
    h_res_1D->Draw("Psame");

    TLatex* r42 = new TLatex(0.15, 0.91, "eAu 18x110 GeV");
    r42->SetNDC();
    r42->SetTextSize(25);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);
    r42->Draw("same");

    TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    r43->Draw("same");

    TLatex* r44_2 = new TLatex(0.15, 0.18, vm_label + " #rightarrow " + daug_label);
    r44_2->SetNDC();
    r44_2->SetTextSize(30);
    r44_2->SetTextFont(43);
    r44_2->SetTextColor(kBlack);
    r44_2->Draw("same");

    TPad* drawPad = new TPad("pad_etalab_11", "pad_etalab_11", 0.16, 0.53, 0.47, 0.83);
    drawPad->SetLeftMargin(0.18);
    drawPad->SetRightMargin(0);
    drawPad->SetTopMargin(0.0);
    drawPad->SetBottomMargin(0.18);
    drawPad->Draw("same");
    drawPad->SetTicks();
    drawPad->SetLogz(1);
    drawPad->cd();
    TH1D* base2 = makeHist("base2", "", "MC", " resolution ", 100, 0, 0.2, kBlack);
    base2->GetYaxis()->SetRangeUser(-10, 2);
    base2->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base2, 3, 3);
    base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize() * 1);
    base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize() * 1);
    base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize() * 1);
    base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize() * 1);
    base2->GetXaxis()->SetNdivisions(5, 5, 0);
    base2->GetYaxis()->SetNdivisions(5, 5, 0);
    base2->Draw();

    h_res->Draw("colzsame");

    c.Print(fmt::format("{}_benchmark-{}-t-resolution.pdf", output_prefix, vm_name).c_str());
  }

  {
    TCanvas c("c1", "c1", 1, 1, 1600, 800);
    c.Divide(4, 2, 0.01, 0.01);
    c.cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_Q2_e->GetXaxis()->SetTitleSize(0.8 * h_Q2_e->GetXaxis()->GetTitleSize());
    h_Q2_e->GetXaxis()->SetLabelSize(0.8 * h_Q2_e->GetXaxis()->GetLabelSize());
    h_Q2_e->GetYaxis()->SetTitleSize(0.8 * h_Q2_e->GetYaxis()->GetTitleSize());
    h_Q2_e->GetYaxis()->SetLabelSize(0.8 * h_Q2_e->GetYaxis()->GetLabelSize());
    h_Q2_e->GetXaxis()->SetTitleOffset(1.6 * h_Q2_e->GetXaxis()->GetTitleOffset());
    h_Q2_e->GetYaxis()->SetTitleOffset(2.0 * h_Q2_e->GetYaxis()->GetTitleOffset());
    h_Q2_e->GetYaxis()->SetTitle("counts");
    h_Q2_e->Draw();
    h_Q2REC_e->SetMarkerStyle(24);
    h_Q2REC_e->Draw("PEsame");
    TLegend* w7 = new TLegend(0.28, 0.7, 0.53, 0.86);
    w7->SetLineColor(kWhite);
    w7->SetFillColor(0);
    w7->SetTextSize(17);
    w7->SetTextFont(45);
    w7->AddEntry(h_energy_MC, "MC ", "L");
    w7->AddEntry(h_energy_REC, "RECO", "P");
    w7->Draw("same");

    c.cd(2);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_y_e->GetXaxis()->SetTitleSize(0.8 * h_y_e->GetXaxis()->GetTitleSize());
    h_y_e->GetXaxis()->SetLabelSize(0.8 * h_y_e->GetXaxis()->GetLabelSize());
    h_y_e->GetYaxis()->SetTitleSize(0.8 * h_y_e->GetYaxis()->GetTitleSize());
    h_y_e->GetYaxis()->SetLabelSize(0.8 * h_y_e->GetYaxis()->GetLabelSize());
    h_y_e->GetXaxis()->SetTitleOffset(1.6 * h_y_e->GetXaxis()->GetTitleOffset());
    h_y_e->GetYaxis()->SetTitleOffset(2.0 * h_y_e->GetYaxis()->GetTitleOffset());
    h_y_e->GetYaxis()->SetTitle("counts");
    h_y_e->Draw();
    h_yREC_e->SetMarkerStyle(24);
    h_yREC_e->Draw("PEsame");
    w7->Draw("same");

    c.cd(3);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_energy_MC->GetXaxis()->SetTitleSize(0.8 * h_energy_MC->GetXaxis()->GetTitleSize());
    h_energy_MC->GetXaxis()->SetLabelSize(0.8 * h_energy_MC->GetXaxis()->GetLabelSize());
    h_energy_MC->GetYaxis()->SetTitleSize(0.8 * h_energy_MC->GetYaxis()->GetTitleSize());
    h_energy_MC->GetYaxis()->SetLabelSize(0.8 * h_energy_MC->GetYaxis()->GetLabelSize());
    h_energy_MC->GetXaxis()->SetTitleOffset(1.6 * h_energy_MC->GetXaxis()->GetTitleOffset());
    h_energy_MC->GetYaxis()->SetTitleOffset(2.5 * h_energy_MC->GetYaxis()->GetTitleOffset());
    h_energy_MC->GetYaxis()->SetTitle("counts");
    h_energy_MC->Draw();
    h_energy_REC->SetMarkerStyle(24);
    h_energy_REC->Draw("PEsame");
    w7->Draw("same");
    // TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
    // w7->SetLineColor(kWhite);
    // w7->SetFillColor(0);
    // w7->SetTextSize(17);
    // w7->SetTextFont(45);
    // w7->AddEntry(h_energy_MC, "MC ", "L");
    // w7->AddEntry(h_energy_REC, "RECO", "P");
    // w7->Draw("same");

    c.cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_Epz_REC->GetXaxis()->SetTitleSize(0.8 * h_Epz_REC->GetXaxis()->GetTitleSize());
    h_Epz_REC->GetXaxis()->SetLabelSize(0.8 * h_Epz_REC->GetXaxis()->GetLabelSize());
    h_Epz_REC->GetYaxis()->SetTitleSize(0.8 * h_Epz_REC->GetYaxis()->GetTitleSize());
    h_Epz_REC->GetYaxis()->SetLabelSize(0.8 * h_Epz_REC->GetYaxis()->GetLabelSize());
    h_Epz_REC->GetXaxis()->SetTitleOffset(1.6 * h_Epz_REC->GetXaxis()->GetTitleOffset());
    h_Epz_REC->GetYaxis()->SetTitleOffset(2.5 * h_Epz_REC->GetYaxis()->GetTitleOffset());
    h_Epz_REC->GetYaxis()->SetTitle("counts");
    h_Epz_REC->SetMarkerStyle(24);
    h_Epz_REC->Draw("PEsame");

    c.cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_energy_calibration_REC->GetXaxis()->SetTitleSize(
        0.8 * h_energy_calibration_REC->GetXaxis()->GetTitleSize());
    h_energy_calibration_REC->GetXaxis()->SetLabelSize(
        0.8 * h_energy_calibration_REC->GetXaxis()->GetLabelSize());
    h_energy_calibration_REC->GetYaxis()->SetTitleSize(
        0.8 * h_energy_calibration_REC->GetYaxis()->GetTitleSize());
    h_energy_calibration_REC->GetYaxis()->SetLabelSize(
        0.8 * h_energy_calibration_REC->GetYaxis()->GetLabelSize());
    h_energy_calibration_REC->GetXaxis()->SetTitleOffset(
        1.6 * h_energy_calibration_REC->GetXaxis()->GetTitleOffset());
    h_energy_calibration_REC->GetYaxis()->SetTitleOffset(
        2.5 * h_energy_calibration_REC->GetYaxis()->GetTitleOffset());
    h_energy_calibration_REC->GetYaxis()->SetTitle("counts");
    h_energy_calibration_REC->GetXaxis()->SetTitle("E_{reco} / E_{mc}");
    h_energy_calibration_REC->SetMarkerStyle(24);
    h_energy_calibration_REC->Draw("PEsame");

    c.cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_EoverP_REC->GetXaxis()->SetTitleSize(0.8 * h_EoverP_REC->GetXaxis()->GetTitleSize());
    h_EoverP_REC->GetXaxis()->SetLabelSize(0.8 * h_EoverP_REC->GetXaxis()->GetLabelSize());
    h_EoverP_REC->GetYaxis()->SetTitleSize(0.8 * h_EoverP_REC->GetYaxis()->GetTitleSize());
    h_EoverP_REC->GetYaxis()->SetLabelSize(0.8 * h_EoverP_REC->GetYaxis()->GetLabelSize());
    h_EoverP_REC->GetXaxis()->SetTitleOffset(1.6 * h_EoverP_REC->GetXaxis()->GetTitleOffset());
    h_EoverP_REC->GetYaxis()->SetTitleOffset(2.5 * h_EoverP_REC->GetYaxis()->GetTitleOffset());
    h_EoverP_REC->GetYaxis()->SetTitle("counts");
    h_EoverP_REC->SetMarkerStyle(24);
    h_EoverP_REC->Draw("PEsame");

    c.cd(7);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_ClusOverHit_REC->GetXaxis()->SetTitleSize(0.8 *
                                                h_ClusOverHit_REC->GetXaxis()->GetTitleSize());
    h_ClusOverHit_REC->GetXaxis()->SetLabelSize(0.8 *
                                                h_ClusOverHit_REC->GetXaxis()->GetLabelSize());
    h_ClusOverHit_REC->GetYaxis()->SetTitleSize(0.8 *
                                                h_ClusOverHit_REC->GetYaxis()->GetTitleSize());
    h_ClusOverHit_REC->GetYaxis()->SetLabelSize(0.8 *
                                                h_ClusOverHit_REC->GetYaxis()->GetLabelSize());
    h_ClusOverHit_REC->GetXaxis()->SetTitleOffset(1.6 *
                                                  h_ClusOverHit_REC->GetXaxis()->GetTitleOffset());
    h_ClusOverHit_REC->GetYaxis()->SetTitleOffset(2.5 *
                                                  h_ClusOverHit_REC->GetYaxis()->GetTitleOffset());
    h_ClusOverHit_REC->GetYaxis()->SetTitle("counts");
    h_ClusOverHit_REC->SetMarkerStyle(24);
    h_ClusOverHit_REC->Draw("PEsame");

    c.cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_emHits_position_REC->GetXaxis()->SetTitleSize(
        0.8 * h_emHits_position_REC->GetXaxis()->GetTitleSize());
    h_emHits_position_REC->GetXaxis()->SetLabelSize(
        0.8 * h_emHits_position_REC->GetXaxis()->GetLabelSize());
    h_emHits_position_REC->GetYaxis()->SetTitleSize(
        0.8 * h_emHits_position_REC->GetYaxis()->GetTitleSize());
    h_emHits_position_REC->GetYaxis()->SetLabelSize(
        0.8 * h_emHits_position_REC->GetYaxis()->GetLabelSize());
    h_emHits_position_REC->GetXaxis()->SetTitleOffset(
        1.6 * h_emHits_position_REC->GetXaxis()->GetTitleOffset());
    h_emHits_position_REC->GetYaxis()->SetTitleOffset(
        2. * h_emHits_position_REC->GetYaxis()->GetTitleOffset());
    h_emHits_position_REC->GetYaxis()->CenterTitle();
    h_emHits_position_REC->GetYaxis()->SetTitle("y(mm)");
    h_emHits_position_REC->Draw("colzsame");

    c.Print(fmt::format("{}_benchmark-{}-DIS-kinematics.pdf", output_prefix, vm_name).c_str());
  }

  return 0;
}
