#include "dvmp.h"
#include "plot.h"

#include <benchmark.h>
#include <mt.h>
#include <util.h>

#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <fmt/color.h>
#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

// Run VM invariant-mass-based benchmarks on an input reconstruction file for
// a desired vector meson (e.g. jpsi) and a desired decay particle (e.g. muon)
// Output figures are written to our output prefix (which includes the output
// file prefix), and labeled with our detector name.
// FIXME: MC does not trace back into particle history. Need to fix that
int vm_invar(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string vm_name       = config["vm_name"];
  const std::string decay_name    = config["decay"];
  const std::string detector      = config["detector"];
  std::string       output_prefix = config["output_prefix"];
  const std::string test_tag      = config["test_tag"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running VM invariant mass analysis...\n");
  fmt::print(" - Vector meson: {}\n", vm_name);
  fmt::print(" - Decay particle: {}\n", decay_name);
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix: {}\n", output_prefix);
  
  std::string VarName[4] = {"y", "Q2", "x", "t"};
  double width_target[4] = {.4, .09, .35, .07};
  //============================== test definition ==============================
  std::vector<eic::util::Test> Tests;
  for(int i = 0 ; i < 4 ; i++){
    eic::util::Test resolution_test_tmp{
      {{"name", fmt::format("{}_{}_resolution", test_tag, VarName[i])},
       {"title",
        fmt::format("{} Resolution for {} -> {} events with {}", VarName[i], vm_name, decay_name, detector)},
       {"description", fmt::format("{} resolution: relative difference with Gaussian fit", VarName[i])},
       {"quantity", "resolution"},
       {"target", fmt::format("{}", width_target[i])}}};
    Tests.push_back(resolution_test_tmp);
  }
  
  //============================== general settings ==============================
  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT(kNumThreads);
  TH1::SetDefaultSumw2();
  
  // The particles we are looking for. E.g. J/psi decaying into e+e-
  const double vm_mass    = util::get_pdg_mass(vm_name);
  const double decay_mass = util::get_pdg_mass(decay_name);

  // Ensure our output prefix always ends on a dot, a slash or a dash
  if (output_prefix.back() != '.' && output_prefix.back() != '/' && output_prefix.back() != '-') {
    output_prefix += "-";
  }

  // Open our input file file as a dataframe
  ROOT::RDataFrame d{"events", rec_file};
  
  //============================== redef functions ==============================
  auto momenta_sort_sim = [vm_name, decay_name](const std::vector<dd4pod::Geant4ParticleData>& parts){
    return util::momenta_sort_sim(parts, vm_name, decay_name);
  };
  auto momenta_sort_rec = [vm_name, decay_name](const std::vector<eic::ReconstructedParticleData>& parts){
    return util::momenta_sort_rec(parts, vm_name, decay_name);
  };

  //============================== Define analysis flow ==============================
  auto d_im = d.Define("p_rec_sorted", momenta_sort_rec, {"DummyReconstructedParticles"})
                  .Define("p_sim_sorted", momenta_sort_sim, {"mcparticles2"})
                  .Define("N", "p_rec_sorted.size()")
                  .Define("invariant_quantities_rec", util::calc_inv_quant, {"p_rec_sorted"})
                  .Define("invariant_quantities_sim", util::calc_inv_quant, {"p_sim_sorted"})
                  .Define("y_rec", util::get_y, {"invariant_quantities_rec"})
                  .Define("Q2_rec", util::get_Q2, {"invariant_quantities_rec"})
                  .Define("x_rec", util::get_x, {"invariant_quantities_rec"})
                  .Define("t_rec", util::get_t, {"invariant_quantities_rec"})
                  .Define("y_sim", util::get_y, {"invariant_quantities_sim"})
                  .Define("Q2_sim", util::get_Q2, {"invariant_quantities_sim"})
                  .Define("x_sim", util::get_x, {"invariant_quantities_sim"})
                  .Define("t_sim", util::get_t, {"invariant_quantities_sim"})
                  .Define("y_dif", "y_rec - y_sim")
                  .Define("Q2_dif", "Q2_rec - Q2_sim")
                  .Define("x_dif", "x_rec - x_sim")
                  .Define("t_dif", "t_rec - t_sim")
                  .Define("y_rdf", "(y_rec - y_sim)/y_sim")
                  .Define("Q2_rdf", "(Q2_rec - Q2_sim)/Q2_sim")
                  .Define("x_rdf", "(x_rec - x_sim)/x_sim")
                  .Define("t_rdf", "(t_rec - t_sim)/t_sim");
  
  //==============================hist def==============================                
  //ranges
  double range_l[4][4] = {{0., 0., -0.3, -1.5}, { 0.,  0., -1., -0.3}, {0.0, 0.0, -0.01, -1.}, {-1., -1., -0.1, -0.5}};
  double range_h[4][4] = {{1., 1.,  0.3,  1.5}, {15., 15.,  1.,  0.3}, {0.1, 0.1,  0.01,  1.}, { 0.,  0.,  0.1,  0.5}};
  
  //strings used
  std::string VarCate[5] = {"sim", "rec", "dif", "rdf", "svr"};
  std::string histName[4][5];
  std::string histTitles[4][5];
  std::string RawhistName[4][4];
  
  TH1D* h_Var1D[4][4];
  TH2D* h_Var2D[4];
  for(int i = 0 ; i < 4 ; i++){
    for(int j = 0 ; j < 4 ; j++){
      //construct histName
      histName[i][j] = "h_" + VarName[i] + "_" + VarCate[j];
      //construct histTitles
      histTitles[i][j] = ";";
      if(j > 1) histTitles[i][j] = histTitles[i][j] + "#Delta";
      if(i==1){
        histTitles[i][j] = histTitles[i][j] + "Q^{2}";
      }else{
        histTitles[i][j] = histTitles[i][j] + VarName[i];
      }
      if(j==3){
        histTitles[i][j] = histTitles[i][j] + "/";
        if(i==1){
          histTitles[i][j] = histTitles[i][j] + "Q^{2}";
        }else{
          histTitles[i][j] = histTitles[i][j] + VarName[i];
        }
      }
      histTitles[i][j] = histTitles[i][j] + ";#";
      //construct RawhistName
      RawhistName[i][j] = VarName[i] + "_" + VarCate[j];
      //get histograms
      auto h_tmp = d_im.Histo1D({fmt::format("{}_tmp", histName[i][j]).c_str(), histTitles[i][j].c_str(), 50, range_l[i][j], range_h[i][j]}, RawhistName[i][j].c_str());
      auto& htmp = *h_tmp;
      h_Var1D[i][j] = (TH1D*)htmp.Clone(histName[i][j].c_str());
    }
    //2d histogram
    histName[i][4] = "h_" + VarName[i] + "_" + VarCate[4];
    if(i==1){
      histTitles[i][4] = ";Q^{2}_{sim};Q^{2}_{rec};#";
    }else{
     histTitles[i][4] = fmt::format(";{}_{sim};{}_{rec};#", VarName[i], VarName[i]); 
    }
    auto h_tmp = d_im.Histo2D({fmt::format("{}_tmp", histName[i][4]).c_str(), histTitles[i][4].c_str(), 50, range_l[i][0], range_h[i][0], 50, range_l[i][0], range_h[i][0]}, RawhistName[i][0].c_str(), RawhistName[i][1].c_str());
    auto& htmp = *h_tmp;
    h_Var2D[i] = (TH2D*)htmp.Clone(histName[i][4].c_str());
  }
  double nEvents = h_Var1D[0][0]->Integral(0, -1);
  
  //==============================fit and plot==============================
  TFitResultPtr myFitPtr[4][2];
  TF1* myf[4][2];
  TText* tptr[4][4];
  TPaveText* t[4][4];
  for(int i = 0 ; i < 4 ; i++){
    TCanvas* ctmp = new TCanvas("ctmp", "ctmp", 1200, 900);
    ctmp->Divide(2, 2, 0.001, 0.001);
    //for pad 1: sim rec overlay
    h_Var1D[i][0]->SetLineColor(plot::kMpBlue);
    h_Var1D[i][0]->SetLineWidth(2);
    h_Var1D[i][0]->GetXaxis()->CenterTitle();
    h_Var1D[i][1]->SetLineColor(plot::kMpOrange);
    h_Var1D[i][1]->SetLineWidth(1);
    
    //for pad 2: sim~rec 2d
    //place holder
    
    //for pad 3: rec - sim
    h_Var1D[i][2]->SetLineColor(plot::kMpGrey);
    h_Var1D[i][2]->SetLineWidth(1);
    h_Var1D[i][2]->GetXaxis()->CenterTitle();
    
    //for pad 4: (rec - sim)/sim
    h_Var1D[i][3]->SetLineColor(plot::kMpGrey);
    h_Var1D[i][3]->SetLineWidth(1);
    h_Var1D[i][3]->GetXaxis()->CenterTitle();
    
    //initialize myf
    for(int j = 0 ; j < 2 ; j++){
      myf[i][j] = new TF1(Form("myf_%d_%d", i, j), "[2]*TMath::Gaus(x, [0], [1], 0)", range_l[i][j+2], range_h[i][j+2]);
      myf[i][j]->SetParameters(0., 0.25, nEvents/10.);
      myf[i][j]->SetParLimits(0, -0.5, 0.5);
      myf[i][j]->SetParLimits(1, 0., 1.0);
      myf[i][j]->SetParLimits(2, 0., nEvents*10.);
      myf[i][j]->SetNpx(1000);
      myf[i][j]->SetLineColor(plot::kMpRed);
      myf[i][j]->SetLineStyle(7);
    }
    
    //fit and plot
    for(int j = 0 ; j < 4 ; j++){
      ctmp->cd(j+1);
      t[i][j] = new TPaveText(.6, .8417, .9, .925, "NB NDC");
      t[i][j]->SetFillColorAlpha(kWhite, 0);
      t[i][j]->SetTextFont(43);
      t[i][j]->SetTextSize(25);
      switch(j){
        case 0://sim rec overlay
          h_Var1D[i][0]->Draw("hist");
          h_Var1D[i][1]->Draw("hist same");
          tptr[i][j] = t[i][j]->AddText("simulation");
          tptr[i][j]->SetTextColor(plot::kMpBlue);
          tptr[i][j] = t[i][j]->AddText("reconstructed");
          tptr[i][j]->SetTextColor(plot::kMpOrange);
          break;
        case 1:
          break;//2d
        case 2://dx
          h_Var1D[i][2]->Draw("hist");
          myFitPtr[i][0] = h_Var1D[i][2]->Fit(myf[i][0], "S 0", "", range_l[i][j], range_h[i][j]);
          myf[i][0]->Draw("same");
          if(i==1){
            tptr[i][j] = t[i][j]->AddText("#DeltaQ^{2}");
          }else{
            tptr[i][j] = t[i][j]->AddText(fmt::format("#Delta{}", VarName[i]).c_str());
          }
          tptr[i][j]->SetTextColor(1);
          break;
        case 3://dx/x
          h_Var1D[i][3]->Draw("hist");
          myFitPtr[i][1] = h_Var1D[i][3]->Fit(myf[i][1], "S 0", "", range_l[i][j], range_h[i][j]);
          myf[i][1]->Draw("same");
          if(i==1){
            tptr[i][j] = t[i][j]->AddText("#DeltaQ^{2}/Q^{2}");
          }else{
            tptr[i][j] = t[i][j]->AddText(fmt::format("#Delta{}/{}", VarName[i], VarName[i]).c_str());
          }
          tptr[i][j]->SetTextColor(1);
          break;
        default:
          break;
      }
      plot::draw_label(10, 100, detector);
      t[i][j]->Draw();
    }
    ctmp->Print(fmt::format("{}{}.png", output_prefix, VarName[i]).c_str());
    delete ctmp;
  }
  
  //============================== pass the test results ==============================
  for(int i = 0 ; i < 4 ; i++){
    double width = myf[i][1]->GetParameter(1);
    if(myFitPtr[i][1]->Status()!=0){
      Tests[i].error(-1);
    }else if(width > width_target[i]){
      Tests[i].fail(width);
    }else{
      Tests[i].pass(width);
    }
  }
  
  // write out our test data
  eic::util::write_test(Tests, fmt::format("{}invar.json", output_prefix));
  
  return 0;
}
