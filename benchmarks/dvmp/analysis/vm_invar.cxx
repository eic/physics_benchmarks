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

  // create our test definition
  std::vector<eic::util::Test> Tests;
  eic::util::Test y_resolution_test{
      {{"name", fmt::format("{}_y_resolution", test_tag)},
       {"title",
        fmt::format("y Resolution for {} -> {} events with {}", vm_name, decay_name, detector)},
       {"description", "y resolution: relative difference with Gaussian fit"},
       {"quantity", "resolution"},
       {"target", ".4"}}};
  Tests.push_back(y_resolution_test);
  
  eic::util::Test Q2_resolution_test{
      {{"name", fmt::format("{}_Q2_resolution", test_tag)},
       {"title",
        fmt::format("Q^2 Resolution for {} -> {} events with {}", vm_name, decay_name, detector)},
       {"description", "Q^2 resolution: relative difference with Gaussian fit"},
       {"quantity", "resolution"},
       {"target", ".09"}}};
  Tests.push_back(Q2_resolution_test);
  
  eic::util::Test x_resolution_test{
      {{"name", fmt::format("{}_x_resolution", test_tag)},
       {"title",
        fmt::format("x Resolution for {} -> {} events with {}", vm_name, decay_name, detector)},
       {"description", "x resolution: relative difference with Gaussian fit"},
       {"quantity", "resolution"},
       {"target", ".35"}}};
  Tests.push_back(x_resolution_test);
  
  eic::util::Test t_resolution_test{
      {{"name", fmt::format("{}_t_resolution", test_tag)},
       {"title",
        fmt::format("t Resolution for {} -> {} events with {}", vm_name, decay_name, detector)},
       {"description", "t resolution: relative difference with Gaussian fit"},
       {"quantity", "resolution"},
       {"target", ".07"}}};
  Tests.push_back(t_resolution_test);

  double width_target[4] = {0.4, 0.09, 0.35, 0.07};
  TH1::SetDefaultSumw2();
  
  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT(kNumThreads);

  // The particles we are looking for. E.g. J/psi decaying into e+e-
  const double vm_mass    = util::get_pdg_mass(vm_name);
  const double decay_mass = util::get_pdg_mass(decay_name);

  // Ensure our output prefix always ends on a dot, a slash or a dash
  if (output_prefix.back() != '.' && output_prefix.back() != '/' && output_prefix.back() != '-') {
    output_prefix += "-";
  }

  // Open our input file file as a dataframe
  ROOT::RDataFrame d{"events", rec_file};

  // utility lambda functions to bind the vector meson and decay particle
  // types
  
  auto momenta_sort_sim = [vm_name, decay_name](const std::vector<dd4pod::Geant4ParticleData>& parts){
    return util::momenta_sort_sim(parts, vm_name, decay_name);
  };
  auto momenta_sort_rec = [vm_name, decay_name](const std::vector<eic::ReconstructedParticleData>& parts){
    return util::momenta_sort_rec(parts, vm_name, decay_name);
  };
  //====================================================================

  // Define analysis flow
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
                  .Define("y_dif", "(y_rec - y_sim)/y_sim")
                  .Define("Q2_dif", "(Q2_rec - Q2_sim)/Q2_sim")
                  .Define("x_dif", "(x_rec - x_sim)/x_sim")
                  .Define("t_dif", "(t_rec - t_sim)/t_sim");
                  
  //================================================================
  //Factorizeation
  double func_range[4] = {1.5, 0.3, 1., 2.};
  double hist_range_l[4] = {0., 0., 0., -1.};
  double hist_range_h[4] = {1., 15., 0.1, 0.};
  
  std::string VarName[4] = {"y", "Q2", "x", "t"};
  /*std::string histName[4];
  std::string histTitle[4];
  std::string RawHist[4];
  for(int i = 0 ; i < 4 ; i++){
    histName[i] = "h_" + VarName[i] + "_sim_test";
    if(i!=1){
    }else{
      histTitle[i] = ";" + VarName[i] + ";#";
      histTitle[i] = ";Q^{2};#";
    }
    RawHist[i] = VarName[i] + "_sim";
  }
  
  TH1D* h_sim[4];
  {
    int i = 0;
    cout<<"================"<<histName[i]<<"================"<<endl;
    //auto h_tmp = d_im.Histo1D({histName[i], ";y;#", 50, hist_range_l[i], hist_range_h[i]}, "y_sim");          //using string variable
    //auto h_tmp = d_im.Histo1D({"h_y_sim_test", ";y;#", 50, hist_range_l[i], hist_range_h[i]}, "y_sim");       //directly quote the string
  }*/
  //==================================================================

  // Define output histograms
  //auto h_nu_sim = d_im.Histo1D({"h_nu_sim", ";#nu/1000;#", 100, 0., 2.}, "nu_sim");
  auto h_Q2_sim = d_im.Histo1D({"h_Q2_sim", ";Q^{2};#", 50, 0., 15.}, "Q2_sim");
  auto h_x_sim  = d_im.Histo1D({"h_x_sim", ";x;#", 50, 0., 0.1}, "x_sim");
  auto h_y_sim  = d_im.Histo1D({"h_y_sim", ";y;#", 50, 0., 1.}, "y_sim");
  auto h_t_sim  = d_im.Histo1D({"h_t_sim", ";t;#", 50, -1., 0.}, "t_sim");
  
  //auto h_nu_rec = d_im.Histo1D({"h_nu_rec", ";#nu/1000;#", 100, 0., 2.}, "nu_rec");
  auto h_Q2_rec = d_im.Histo1D({"h_Q2_rec", ";Q^{2};#", 50, 0., 15.}, "Q2_rec");
  auto h_x_rec  = d_im.Histo1D({"h_x_rec", ";x;#", 50, 0., 0.1}, "x_rec");
  auto h_y_rec  = d_im.Histo1D({"h_y_rec", ";y;#", 50, 0., 1.}, "y_rec");
  auto h_t_rec  = d_im.Histo1D({"h_t_rec", ";t;#", 50, -1., 0.}, "t_rec");
  
  auto h_y_dif   = d_im.Histo1D({"h_y_dif",  ";#Deltay/y;#",     50, -1.5, 1.5}, "y_dif");
  auto h_Q2_dif  = d_im.Histo1D({"h_Q2_dif", ";#DeltaQ^{2}/Q^{2};#", 50, -0.3, 0.3}, "Q2_dif");
  auto h_x_dif   = d_im.Histo1D({"h_x_dif",  ";#Deltax/x;#",     50, -1., 1.}, "x_dif");
  auto h_t_dif   = d_im.Histo1D({"h_t_dif",  ";#Deltat/t;#",     50, -0.5, 0.5}, "t_dif");
  
  double nEvents = h_y_dif->Integral(0, -1);
  
  TH1D* hist_sim[4] = {&(*h_y_sim), &(*h_Q2_sim), &(*h_x_sim), &(*h_t_sim)};
  TH1D* hist_rec[4] = {&(*h_y_rec), &(*h_Q2_rec), &(*h_x_rec), &(*h_t_rec)};
  TH1D* hist_dif[4] = {&(*h_y_dif), &(*h_Q2_dif), &(*h_x_dif), &(*h_t_dif)};
  TFitResultPtr myFitPtr[4];
  TF1* myf[4];
  TText* tptr[4][3];
  TPaveText* t[4][3];
  
  for(int i = 0 ; i < 4 ; i++){
    TCanvas* ctmp = new TCanvas("ctmp", "ctmp", 1800, 600);
    ctmp->Divide(3, 1, 0.001, 0.001);
    //for pad 1
    hist_sim[i]->SetLineColor(plot::kMpBlue);
    hist_sim[i]->SetLineWidth(2);
    hist_sim[i]->GetXaxis()->CenterTitle();
    hist_rec[i]->SetLineColor(plot::kMpOrange);
    hist_rec[i]->SetLineWidth(1);
    //for pad 2
    hist_dif[i]->SetLineColor(plot::kMpGrey);
    hist_dif[i]->SetLineWidth(1);
    hist_dif[i]->GetXaxis()->CenterTitle();
    myf[i] = new TF1(Form("myf_%d", i), "[2]*TMath::Gaus(x, [0], [1], 0)", -func_range[i], func_range[i]);
    myf[i]->SetParameters(0., 0.25, nEvents/10.);
    myf[i]->SetParLimits(0, -0.5, 0.5);
    myf[i]->SetParLimits(1, 0., 1.0);
    myf[i]->SetParLimits(2, 0., nEvents*10.);
    myf[i]->SetNpx(1000);
    myf[i]->SetLineColor(plot::kMpRed);
    myf[i]->SetLineStyle(7);
    //for pad 3
    
    //factorized part
    for(int j = 0 ; j < 2 ; j++){
      ctmp->cd(j+1);
      t[i][j] = new TPaveText(.6, .8417, .9, .925, "NB NDC");
      t[i][j]->SetFillColorAlpha(kWhite, 0);
      t[i][j]->SetTextFont(43);
      t[i][j]->SetTextSize(25);
      switch(j){
        case 0:
          hist_sim[i]->Draw("hist");
          hist_rec[i]->Draw("hist same");
          tptr[i][j] = t[i][j]->AddText(VarName[i].c_str());
          break;
        case 1:
          hist_dif[i]->Draw("hist");
          myFitPtr[i] = hist_dif[i]->Fit(myf[i], "S 0", "", -func_range[i], func_range[i]);
          myf[i]->Draw("same");
          tptr[i][j] = t[i][j]->AddText(fmt::format("#Delta{}/{}", VarName[i], VarName[i]).c_str());
          break;
        case 2:
          break;
        default:
          break;
      }
      plot::draw_label(10, 100, detector);
      tptr[i][j]->SetTextColor(plot::kMpOrange);
      t[i][j]->Draw();
    }
    ctmp->Print(fmt::format("{}-{}.png", output_prefix, VarName[i]).c_str());
    delete ctmp;
  }
  
    
  // Plot our histograms.
  // TODO: to start I'm explicitly plotting the histograms, but want to
  // factorize out the plotting code moving forward.
  //Before factorizing==========================================================================================
    /*TFitResultPtr myFitPtr[4];
    TF1* myf[4];
    for(int i = 0 ; i < 4 ; i++){
        myf[i] = new TF1(Form("myf_%d", i), "[2]*TMath::Gaus(x, [0], [1], 0)", -func_range[i], func_range[i]);
        myf[i]->SetParameters(0., 0.25, nEvents/10.);
        myf[i]->SetParLimits(0, -0.5, 0.5);
        myf[i]->SetParLimits(1, 0., 1.0);
        myf[i]->SetParLimits(2, 0., nEvents*10.);
        myf[i]->SetNpx(1000);
        myf[i]->SetLineColor(2);
        myf[i]->SetLineStyle(7);
    }*/
    
    
    // Print canvas to output file
    /*
    TCanvas c{"canvas2", "canvas2", 1200, 900};
    c.Divide(2, 2, 0.0001, 0.0001);
    //============================================================================
    //pad 1 nu_dif
    c.cd(1);
    auto& hy_dif = *h_y_dif;
    // histogram style
    hy_dif.SetLineColor(plot::kMpOrange);
    hy_dif.SetLineWidth(1);
    // axes
    hy_dif.GetXaxis()->CenterTitle();
    // draw everything
    hy_dif.DrawClone("hist");
    myFitPtr[0] = hy_dif.Fit(myf[0], "S 0", "", -func_range[0], func_range[0]);
    myf[0]->Draw("same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector);
    TText* tptr1;
    auto   t1 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t1->SetFillColorAlpha(kWhite, 0);
    t1->SetTextFont(43);
    t1->SetTextSize(25);
    tptr1 = t1->AddText("#Deltay/y");
    tptr1->SetTextColor(plot::kMpOrange);
    t1->Draw();


    // pad 2 Q2_dif
    c.cd(2);
    auto& hQ2_dif = *h_Q2_dif;
    // histogram style
    hQ2_dif.SetLineColor(plot::kMpOrange);
    hQ2_dif.SetLineWidth(1);
    // axes
    hQ2_dif.GetXaxis()->CenterTitle();
    // draw everything
    hQ2_dif.DrawClone("hist");
    myFitPtr[1] = hQ2_dif.Fit(myf[1], "S 0", "", -func_range[1], func_range[1]);
    myf[1]->Draw("same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector);
    TText* tptr2;
    auto   t2 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t2->SetFillColorAlpha(kWhite, 0);
    t2->SetTextFont(43);
    t2->SetTextSize(25);
    tptr2 = t2->AddText("#DeltaQ^{2}/Q^{2}");
    tptr2->SetTextColor(plot::kMpOrange);
    t2->Draw();
    
    // pad 3 x_dif
    c.cd(3);
    auto& hx_dif = *h_x_dif;
    // histogram style
    hx_dif.SetLineColor(plot::kMpOrange);
    hx_dif.SetLineWidth(1);
    // axes
    hx_dif.GetXaxis()->CenterTitle();
    // draw everything
    hx_dif.DrawClone("hist");
    myFitPtr[2] = hx_dif.Fit(myf[2], "S 0", "", -func_range[2], func_range[2]);
    myf[2]->Draw("same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector);
    TText* tptr3;
    auto   t3 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t3->SetFillColorAlpha(kWhite, 0);
    t3->SetTextFont(43);
    t3->SetTextSize(25);
    tptr3 = t3->AddText("#Deltax/x");
    tptr3->SetTextColor(plot::kMpOrange);
    t3->Draw();
    
    // pad 4 t_dif
    c.cd(4);
    auto& ht_dif = *h_t_dif;
    // histogram style
    ht_dif.SetLineColor(plot::kMpOrange);
    ht_dif.SetLineWidth(1);
    // axes
    ht_dif.GetXaxis()->CenterTitle();
    // draw everything
    ht_dif.DrawClone("hist");
    myFitPtr[3] = ht_dif.Fit(myf[3], "S 0", "", -func_range[3], func_range[3]);
    myf[3]->Draw("same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector);
    TText* tptr4;
    auto   t4 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t4->SetFillColorAlpha(kWhite, 0);
    t4->SetTextFont(43);
    t4->SetTextSize(25);
    tptr4 = t4->AddText("#Deltat/t");
    tptr4->SetTextColor(plot::kMpOrange);
    t4->Draw();
    //============================================================================
    c.Print(fmt::format("{}InvariantQuantities.png", output_prefix).c_str());*/
    
  //Before factorizing==========================================================================================

  for(int i = 0 ; i < 4 ; i++){
    double width = myf[i]->GetParameter(1);
    if(myFitPtr[i]->Status()!=0){
      Tests[i].error(-1);
    }else if(width > width_target[i]){
      Tests[i].fail(width);
    }else{
      Tests[i].pass(width);
    }
  }
  
  // write out our test data
  eic::util::write_test(Tests, fmt::format("{}invar.json", output_prefix));

  // That's all!
  return 0;
}
