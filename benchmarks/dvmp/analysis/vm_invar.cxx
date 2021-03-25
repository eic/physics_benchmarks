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
                  .Define("y_diff", "(y_rec - y_sim)/y_sim")
                  .Define("Q2_diff", "(Q2_rec - Q2_sim)/Q2_sim")
                  .Define("x_diff", "(x_rec - x_sim)/x_sim")
                  .Define("t_diff", "(t_rec - t_sim)/t_sim");
                  
  //================================================================
  //Factorized
  TString VarName[4] = {"y", "Q2", "x", "t"};
  
  auto h_sim[4];
  auto h_rec[4];
  auto h_diff[4];
  
  double fun_range[4] = {1.5, 0.3, 1., 2.};
  double hist_range_l[4] = {0., 0., 0., -1.};
  double hist_range_h[4] = {1., 15., 0.1, 0.};
  
  h_sim[0] = (TH1D*)d_im.Histo1D({"h_" + VarName[0] + "_sim", ";Q^{2};#", 50, hist_range_l[0], hist_range_h[0]}, VarName[0] + "_sim");
  
  /*for(int i = 0 ; i < 4 ; i++){
    if(i==1){
      h_sim[i] = (TH1D*)d_im.Histo1D({"h_" + VarName[i] + "_sim", ";Q^{2};#", 50, hist_range_l[i], hist_range_h[i]}, VarName[i] + "_sim");
    }else{
      h_sim[i] = (TH1D*)d_im.Histo1D({"h_" + VarName[i] + "_sim", ";" + VarName[i] + ";#", 50, hist_range_l[i], hist_range_h[i]}, VarName[i] + "_sim");
    }
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
  
  
  
  auto h_y_diff   = d_im.Histo1D({"h_y_diff",  ";#Deltay/y;#",     50, -1.5, 1.5}, "y_diff");
  auto h_Q2_diff  = d_im.Histo1D({"h_Q2_diff", ";#DeltaQ^{2}/Q^{2};#", 50, -0.3, 0.3}, "Q2_diff");
  auto h_x_diff   = d_im.Histo1D({"h_x_diff",  ";#Deltax/x;#",     50, -1., 1.}, "x_diff");
  auto h_t_diff   = d_im.Histo1D({"h_t_diff",  ";#Deltat/t;#",     50, -0.5, 0.5}, "t_diff");
  
  double nEvents = h_y_diff->Integral(0, -1);
  // Plot our histograms.
  // TODO: to start I'm explicitly plotting the histograms, but want to
  // factorize out the plotting code moving forward.
    TFitResultPtr myFitPtr[4];
    TF1* myf[4];
    for(int i = 0 ; i < 4 ; i++){
        myf[i] = new TF1(Form("myf_%d", i), "[2]*TMath::Gaus(x, [0], [1], 0)", -fun_range[i], fun_range[i]);
        myf[i]->SetParameters(0., 0.25, nEvents/10.);
        myf[i]->SetParLimits(0, -0.5, 0.5);
        myf[i]->SetParLimits(1, 0., 1.0);
        /*if(i==3){
          myf[i]->SetParameter(1, 0.1);
          //myf[i]->SetParameter(2, nEvents);
        }*/
        myf[i]->SetParLimits(2, 0., nEvents*10.);
        myf[i]->SetNpx(1000);
        myf[i]->SetLineColor(2);
        myf[i]->SetLineStyle(7);
    }
    
    
    // Print canvas to output file
    TCanvas c{"canvas2", "canvas2", 1200, 900};
    c.Divide(2, 2, 0.0001, 0.0001);
    //============================================================================
    //pad 1 nu_diff
    c.cd(1);
    auto& hy_diff = *h_y_diff;
    // histogram style
    hy_diff.SetLineColor(plot::kMpOrange);
    hy_diff.SetLineWidth(1);
    // axes
    hy_diff.GetXaxis()->CenterTitle();
    // draw everything
    hy_diff.DrawClone("hist");
    myFitPtr[0] = hy_diff.Fit(myf[0], "S 0", "", -fun_range[0], fun_range[0]);
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


    // pad 2 Q2_diff
    c.cd(2);
    auto& hQ2_diff = *h_Q2_diff;
    // histogram style
    hQ2_diff.SetLineColor(plot::kMpOrange);
    hQ2_diff.SetLineWidth(1);
    // axes
    hQ2_diff.GetXaxis()->CenterTitle();
    // draw everything
    hQ2_diff.DrawClone("hist");
    myFitPtr[1] = hQ2_diff.Fit(myf[1], "S 0", "", -fun_range[1], fun_range[1]);
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
    
    // pad 3 x_diff
    c.cd(3);
    auto& hx_diff = *h_x_diff;
    // histogram style
    hx_diff.SetLineColor(plot::kMpOrange);
    hx_diff.SetLineWidth(1);
    // axes
    hx_diff.GetXaxis()->CenterTitle();
    // draw everything
    hx_diff.DrawClone("hist");
    myFitPtr[2] = hx_diff.Fit(myf[2], "S 0", "", -fun_range[2], fun_range[2]);
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
    
    // pad 4 t_diff
    c.cd(4);
    auto& ht_diff = *h_t_diff;
    // histogram style
    ht_diff.SetLineColor(plot::kMpOrange);
    ht_diff.SetLineWidth(1);
    // axes
    ht_diff.GetXaxis()->CenterTitle();
    // draw everything
    ht_diff.DrawClone("hist");
    myFitPtr[3] = ht_diff.Fit(myf[3], "S 0", "", -fun_range[3], fun_range[3]);
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
    c.Print(fmt::format("{}InvariantQuantities.png", output_prefix).c_str());

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
