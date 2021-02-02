#include "benchmark.hh"
#include "mt.h"
#include "plot.h"
#include "util.h"

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
// TODO: I think it would be better to pass small json configuration file to
//       the test, instead of this ever-expanding list of function arguments.
// FIXME: MC does not trace back into particle history. Need to fix that
int vm_invar(const std::string& config_name) {
  // read our configuration
  std::ifstream config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file = config["rec_file"];
  const std::string vm_name = config["vm_name"];
  const std::string decay_name = config["decay"];
  const std::string detector = config["detector"];
  std::string output_prefix = config["output_prefix"];
  const std::string test_tag = config["test_tag"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running VM invariant mass analysis...\n");
  fmt::print(" - Vector meson: {}\n", vm_name);
  fmt::print(" - Decay particle: {}\n", decay_name);
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - output prefix: {}\n", output_prefix);

  // create our test definition
  // test_tag
  eic::util::Test vm_mass_resolution_test{
      {{"name",
        fmt::format("{}_{}_{}_mass_resolution", test_tag, vm_name, decay_name)},
       {"title",
        fmt::format("{} -> {} Invariant Mass Resolution", vm_name, decay_name)},
       {"description", "Invariant Mass Resolution calculated from raw "
                       "tracking data using a Gaussian fit."},
       {"quantity", "resolution"},
       {"target", ".1"}}};

  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT(kNumThreads);

  // The particles we are looking for. E.g. J/psi decaying into e+e-
  const double vm_mass = util::get_pdg_mass(vm_name);
  const double decay_mass = util::get_pdg_mass(decay_name);

  // Ensure our output prefix always ends on a dot, a slash or a dash
  if (output_prefix.back() != '.' && output_prefix.back() != '/' &&
      output_prefix.back() != '-') {
    output_prefix += "-";
  }

  // Open our input file file as a dataframe
  ROOT::RDataFrame d{"events", rec_file};

  // utility lambda functions to bind the vector meson and decay particle
  // types
  auto momenta_from_tracking =
      [decay_mass](const std::vector<eic::TrackParametersData>& tracks) {
        return util::momenta_from_tracking(tracks, decay_mass);
      };
  auto calc_inv_quant_rec = 
      [vm_mass](const std::vector<ROOT::Math::PxPyPzMVector>& parts) {
        return util::calc_inv_quant_rec(parts, vm_mass);
      };

  //====================================================================
    
  // Define analysis flow
  auto d_im =
      d.Define("p_rec", momenta_from_tracking, {"outputTrackParameters"})
          .Define("N", "p_rec.size()")
          .Define("p_sim", util::momenta_from_simulation, {"mcparticles2"})
          //================================================================
          .Define("invariant_quantities_rec", calc_inv_quant_rec, {"p_rec"})
          .Define("invariant_quantities_sim", util::calc_inv_quant_simu, {"p_sim"})
          .Define("nu_rec" , util::get_nu, {"invariant_quantities_rec"})
          .Define("Q2_rec" , util::get_Q2, {"invariant_quantities_rec"})
          .Define("x_rec" ,  util::get_x, {"invariant_quantities_rec"})
          .Define("t_rec",   util::get_t, {"invariant_quantities_rec"})
          .Define("nu_sim" , util::get_nu, {"invariant_quantities_sim"})
          .Define("Q2_sim" , util::get_Q2, {"invariant_quantities_sim"})
          .Define("x_sim" ,  util::get_x, {"invariant_quantities_sim"})
          .Define("t_sim",   util::get_t, {"invariant_quantities_sim"});
          //================================================================

  // Define output histograms

  auto h_nu_rec = d_im.Histo1D(
      {"h_nu_rec", ";#nu/1000;#", 100, 0., 2.}, "nu_rec");
  auto h_Q2_rec = d_im.Histo1D(
      {"h_Q2_rec", ";Q^{2};#", 100, 0., 15.}, "Q2_rec");
  auto h_x_rec = d_im.Histo1D(
      {"h_x_rec", ";x;#", 100, 0., 0.1}, "x_rec");
  auto h_t_rec = d_im.Histo1D(
      {"h_t_rec", ";t;#", 100, -1., 0.}, "t_rec");

  
  auto h_nu_sim = d_im.Histo1D(
      {"h_nu_sim", ";#nu/1000;#", 100, 0., 2.}, "nu_sim");
  auto h_Q2_sim = d_im.Histo1D(
      {"h_Q2_sim", ";Q^{2};#", 100, 0., 15.}, "Q2_sim");
  auto h_x_sim = d_im.Histo1D(
      {"h_x_sim", ";x;#", 100, 0., 0.1}, "x_sim");
  auto h_t_sim = d_im.Histo1D(
      {"h_t_sim", ";t;#", 100, -1., 0.}, "t_sim");


  // Plot our histograms.
  // TODO: to start I'm explicitly plotting the histograms, but want to
  // factorize out the plotting code moving forward.
  {
    
    // Print canvas to output file
    
    TCanvas c{"canvas2", "canvas2", 1200, 1200};
    c.Divide(2, 2, 0.0001, 0.0001);
    //pad 1 nu
    c.cd(1);
    //gPad->SetLogx(false);
    //gPad->SetLogy(false);
    auto& hnu_rec = *h_nu_rec;
    auto& hnu_sim = *h_nu_sim;
    // histogram style
    hnu_rec.SetLineColor(plot::kMpOrange);
    hnu_rec.SetLineWidth(2);
    hnu_sim.SetLineColor(plot::kMpBlue);
    hnu_sim.SetLineWidth(2);
    // axes
    hnu_rec.GetXaxis()->CenterTitle();
    //hnu.GetXaxis()->SetTitle("#times1000");
    // draw everything
    hnu_sim.DrawClone("hist");
    hnu_rec.DrawClone("hist same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector, vm_name, "#nu");
    TText* tptr1;
    auto t1 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t1->SetFillColorAlpha(kWhite, 0);
    t1->SetTextFont(43);
    t1->SetTextSize(25);
    tptr1 = t1->AddText("simulated");
    tptr1->SetTextColor(plot::kMpBlue);
    tptr1 = t1->AddText("reconstructed");
    tptr1->SetTextColor(plot::kMpOrange);
    t1->Draw();
    
    //pad 2 Q2
    c.cd(2);
    //gPad->SetLogx(false);
    //gPad->SetLogy(false);
    auto& hQ2_rec = *h_Q2_rec;
    auto& hQ2_sim = *h_Q2_sim;
    // histogram style
    hQ2_rec.SetLineColor(plot::kMpOrange);
    hQ2_rec.SetLineWidth(2);
    hQ2_sim.SetLineColor(plot::kMpBlue);
    hQ2_sim.SetLineWidth(2);
    // axes
    hQ2_rec.GetXaxis()->CenterTitle();
    //hnu.GetXaxis()->SetTitle("#times1000");
    // draw everything
    hQ2_sim.DrawClone("hist");
    hQ2_rec.DrawClone("hist same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector, vm_name, "Q^{2}");
    TText* tptr2;
    auto t2 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t2->SetFillColorAlpha(kWhite, 0);
    t2->SetTextFont(43);
    t2->SetTextSize(25);
    tptr2 = t2->AddText("simulated");
    tptr2->SetTextColor(plot::kMpBlue);
    tptr2 = t2->AddText("reconstructed");
    tptr2->SetTextColor(plot::kMpOrange);
    t2->Draw();

    
    //pad 3 x
    c.cd(3);
    //gPad->SetLogx(false);
    //gPad->SetLogy(false);
    auto& hx_rec = *h_x_rec;
    auto& hx_sim = *h_x_sim;
    // histogram style
    hx_rec.SetLineColor(plot::kMpOrange);
    hx_rec.SetLineWidth(2);
    hx_sim.SetLineColor(plot::kMpBlue);
    hx_sim.SetLineWidth(2);
    // axes
    hx_rec.GetXaxis()->CenterTitle();
    //hnu.GetXaxis()->SetTitle("#times1000");
    // draw everything
    hx_sim.DrawClone("hist");
    hx_rec.DrawClone("hist same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector, vm_name, "x");
    TText* tptr3;
    auto t3 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t3->SetFillColorAlpha(kWhite, 0);
    t3->SetTextFont(43);
    t3->SetTextSize(25);
    tptr3 = t3->AddText("simulated");
    tptr3->SetTextColor(plot::kMpBlue);
    tptr3 = t3->AddText("reconstructed");
    tptr3->SetTextColor(plot::kMpOrange);
    t3->Draw();
    
    //pad 4 t
    c.cd(4);
    //gPad->SetLogx(false);
    //gPad->SetLogy(false);
    auto& ht_rec = *h_t_rec;
    auto& ht_sim = *h_t_sim;
    // histogram style
    ht_rec.SetLineColor(plot::kMpOrange);
    ht_rec.SetLineWidth(2);
    ht_sim.SetLineColor(plot::kMpBlue);
    ht_sim.SetLineWidth(2);
    // axes
    ht_rec.GetXaxis()->CenterTitle();
    //hnu.GetXaxis()->SetTitle("#times1000");
    // draw everything
    ht_sim.DrawClone("hist");
    ht_rec.DrawClone("hist same");
    // FIXME hardcoded beam configuration
    plot::draw_label(10, 100, detector, vm_name, "t");
    TText* tptr4;
    auto t4 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t4->SetFillColorAlpha(kWhite, 0);
    t4->SetTextFont(43);
    t4->SetTextSize(25);
    tptr4 = t4->AddText("simulated");
    tptr4->SetTextColor(plot::kMpBlue);
    tptr4 = t4->AddText("reconstructed");
    tptr4->SetTextColor(plot::kMpOrange);
    t4->Draw();

    
    c.Print(fmt::format("{}InvariantQuantities.png", output_prefix).c_str());
  }

  // TODO we're not actually doing an IM fit yet, so for now just return an
  // error for the test result
  vm_mass_resolution_test.error(-1);

  // write out our test data
  eic::util::write_test(vm_mass_resolution_test,
                           fmt::format("{}vm_invar.json", output_prefix));

  // That's all!
  return 0;
}
