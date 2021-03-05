#include "dis.h"
#include "plot.h"

#include <benchmark.h>
#include <mt.h>
#include <util.h>

#include "ROOT/RDataFrame.hxx"
#include <TFitResult.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <eicd/ReconstructedParticleData.h>
#include <fmt/color.h>
#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

// Beam energies
double eBeamEnergy;
double pBeamEnergy;

// Gaussian variables
double gausMean;
double gausSTD;

// Get a vector of 4-momenta from the reconstructed data.
inline auto
momenta_from_reconstruction_pid(const std::vector<eic::ReconstructedParticleData>& parts)
{
  std::vector<std::pair<ROOT::Math::PxPyPzEVector, int>> momenta(parts.size());
  // transform our reconstructed particle data into 4-momenta
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    return std::make_pair(ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy},
                          part.pid);
  });
  return momenta;
}

// Get a vector of 4-momenta from the simulated data.
inline auto momenta_from_simulation_pid(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<std::pair<ROOT::Math::PxPyPzEVector, int>> momenta(parts.size());
  // transform our simulation particle data into 4-momenta
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    double energy = sqrt(part.psx * part.psx + part.psy * part.psy + part.psz * part.psz +
                         part.mass * part.mass);
    return std::make_pair(ROOT::Math::PxPyPzEVector{part.psx, part.psy, part.psz, energy},
                          part.pdgID);
  });
  return momenta;
}

// Finds particle with highest momentum within the event and returns the 4 Vector
// Currentky there is a cut so that only electrons are retained, but that can be changed
// by chaning the PID
inline auto find_scat(const std::vector<std::pair<ROOT::Math::PxPyPzEVector, int>>& mom)
{
  ROOT::Math::PxPyPzEVector elec_cand = {0, 0, 0, 0};
  int                       partID    = 11; // For electrons
  for (auto mom_cand : mom) {
    if (mom_cand.first.E() > elec_cand.E() && mom_cand.second == partID) {
      elec_cand = mom_cand.first;
    }
  }
  return elec_cand;
}

// Q2 calculation from 4 Vector
// Recall, the electron ebeam is defined as coming from -z and the proton beam from +z
inline auto Q2(ROOT::Math::PxPyPzEVector& mom)
{
  ROOT::Math::PxPyPzEVector beamMom = {0, 0, -eBeamEnergy, eBeamEnergy};
  return -(mom - beamMom).M2();
}

// Multiplies a double by a gaussian distributed number
inline auto randomize(double& inDouble)
{
  TRandom3 rand(0);
  return rand.Gaus(gausMean, gausSTD) * inDouble;
}

int dis_electrons(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string detector      = config["detector"];
  std::string       output_prefix = config["output_prefix"];
  const std::string test_tag      = config["test_tag"];
  eBeamEnergy                     = config["ebeam"];
  pBeamEnergy                     = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DIS electron analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix: {}\n", output_prefix);
  fmt::print(" - test tag: {}\n", test_tag);

  // create our test definition
  // test_tag
  eic::util::Test dis_Q2_resolution{
      {{"name", fmt::format("{}_Q2_resolution", test_tag)},
       {"title", "DIS Q2 resolution"},
       {"description",
        fmt::format("DIS Q2 resolution with {}, estimated using a Gaussian fit.", detector)},
       {"quantity", "resolution (in %)"},
       {"target", "0.1"}}};

  // Run this in multi-threaded mode if desired
  // ROOT::EnableImplicitMT(kNumThreads);

  // PIDs for reference
  // electron = 11
  // pi_0     = 111
  // pi_plus  = 211
  // pi_minus = -211
  // k_0      = 311
  // k_plus   = 321
  // k_minus  = -321
  // proton   = 2212
  // neutron  = 2112

  // Gausian variable declarations
  gausMean = 1.0;
  gausSTD  = 0.2;

  const double electron_mass = util::get_pdg_mass("electron");

  // Ensure our output prefix always ends on a dot, a slash or a dash
  // Necessary when generating output plots
  if (output_prefix.back() != '.' && output_prefix.back() != '/' && output_prefix.back() != '-') {
    output_prefix += "-";
  }

  ROOT::RDataFrame d("events", rec_file);

  // utility lambda functions to bind the reconstructed particle type
  // (as we have no PID yet)
  auto momenta_from_tracking =
      [electron_mass](const std::vector<eic::TrackParametersData>& tracks) {
        return util::momenta_from_tracking(tracks, electron_mass);
      };

  /*
    //Old dataframe
    auto d0 = d.Define("p_rec", momenta_from_tracking, {"outputTrackParameters"})
                  .Define("N", "p_rec.size()")
                  .Define("p_sim", util::momenta_from_simulation, {"mcparticles2"})
                  .Define("mom_sim", util::mom, {"p_sim"})
                  .Define("mom_rec", util::mom, {"p_rec"});
  */

  auto d0 = d.Define("p_recon", momenta_from_reconstruction_pid, {"DummyReconstructedParticles"})
                .Define("elec_recon", find_scat, {"p_recon"})
                .Define("elec_recon_Q2", Q2, {"elec_recon"})
                .Define("elec_recon_Q2_rand", randomize, {"elec_recon_Q2"})
                .Define("p_sim", momenta_from_simulation_pid, {"mcparticles2"})
                .Define("elec_sim", find_scat, {"p_sim"})
                .Define("elec_sim_Q2", Q2, {"elec_sim"})
                .Define("dQ2", "elec_recon_Q2_rand - elec_sim_Q2")
                .Define("dQ2_relative", "(elec_recon_Q2_rand - elec_sim_Q2)/elec_sim_Q2");

  // Testing script
  /*
    //auto dis = d0.Display({"pt_recon", "elec_sim_Q2", "elec_recon_Q2"});
    //dis -> Print();
    std:vector<string> nameStr = {"elec_recon_Q2_rand", "elec_sim_Q2", "dQ2", "dQ2_relative"};
    for (auto name : nameStr)
    {
      printf("%s Min(%f) Mean(%f) Max(%f)\n", name.c_str(), *d0.Min<double>(name),
  *d0.Mean<double>(name), *d0.Max<double>(name) );
    }
  /**/
  // Momentum
  // auto h_mom_sim = d0.Histo1D({"h_mom_sim", "; GeV; counts", 100, 0, 50}, "mom_sim");
  // auto h_mom_rec = d0.Histo1D({"h_mom_rec", "; GeV; counts", 100, 0, 50}, "mom_rec");

  // Q2
  auto h_Q2_sim = d0.Histo1D({"h_Q2_sim", "; GeV^{2}; Counts", 50, 0, 25}, "elec_sim_Q2");
  auto h_Q2_rec = d0.Histo1D({"h_Q2_rec", "; GeV^{2}; Counts", 50, 0, 25}, "elec_recon_Q2_rand");
  auto h_Q2_res =
      d0.Histo1D({"h_Q2_res", "; Q^{2} Resolution (GeV^{2}); Counts", 50, -120, 120}, "dQ2");
  auto h_Q2_rel_res = d0.Histo1D(
      {"h_Q2_rel_res", "; Q^{2} Relative Resolution ; Counts", 50, -2, 2}, "dQ2_relative");

  TFitResultPtr f1  = h_Q2_res->Fit("gaus", "S");
  const double* res = f1->GetParams();
  if (res[2] <= 0.1) {
    dis_Q2_resolution.pass(res[2]);
  } else {
    dis_Q2_resolution.fail(res[2]);
  }
  f1->Print("V");

  auto c = new TCanvas();

  // Plot our histograms.
  // TODO: to start I'm explicitly plotting the histograms, but want to
  // factorize out the plotting code moving forward.
  {
    TCanvas c{"canvas", "canvas", 1200, 1200};
    c.cd();
    // gPad->SetLogx(false);
    gPad->SetLogy(true);
    // auto& h1 = *h_mom_sim;
    // auto& h2 = *h_mom_rec;
    auto& h1 = *h_Q2_sim;
    auto& h2 = *h_Q2_rec;
    // histogram style
    h1.SetLineColor(plot::kMpBlue);
    h1.SetLineWidth(2);
    h2.SetLineColor(plot::kMpOrange);
    h2.SetLineWidth(2);
    // axes
    h1.GetXaxis()->CenterTitle();
    h1.GetYaxis()->CenterTitle();
    // draw everything
    h1.DrawClone("hist");
    h2.DrawClone("hist same");
    // FIXME hardcoded beam configuration
    plot::draw_label(eBeamEnergy, pBeamEnergy, detector);
    TText* tptr1;
    auto   t1 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t1->SetFillColorAlpha(kWhite, 0);
    t1->SetTextFont(43);
    t1->SetTextSize(25);
    tptr1 = t1->AddText("simulated");
    tptr1->SetTextColor(plot::kMpBlue);
    tptr1 = t1->AddText("reconstructed");
    tptr1->SetTextColor(plot::kMpOrange);
    t1->Draw();

    c.Print(fmt::format("{}momentum.png", output_prefix).c_str());
  }
  // Plot Q2 and Relatitive Q2
  {
    TCanvas c{"canvas", "canvas", 1200, 1200};
    c.cd();
    // gPad->SetLogx(false);
    gPad->SetLogy(true);
    // auto& h1 = *h_mom_sim;
    // auto& h2 = *h_mom_rec;
    auto& h1 = *h_Q2_res;
    auto& h2 = *h_Q2_rel_res;
    // histogram style
    h1.SetLineColor(plot::kMpBlue);
    h1.SetLineWidth(2);
    h2.SetLineColor(plot::kMpOrange);
    h2.SetLineWidth(2);
    // axes
    h1.GetXaxis()->CenterTitle();
    h1.GetYaxis()->CenterTitle();
    h2.GetXaxis()->CenterTitle();
    h2.GetYaxis()->CenterTitle();
    // draw everything
    h1.DrawClone("hist");
    plot::draw_label(eBeamEnergy, pBeamEnergy, detector);
    TText* tptr1;
    auto   t1 = new TPaveText(.6, .8417, .9, .925, "NB NDC");
    t1->SetFillColorAlpha(kWhite, 0);
    t1->SetTextFont(43);
    t1->SetTextSize(25);
    tptr1 = t1->AddText("Q^{2} Resolution");
    tptr1->SetTextColor(plot::kMpBlue);
    t1->Draw();
    c.Print(fmt::format("{}Q2_resolution.png", output_prefix).c_str());
    t1->Clear();

    h2.DrawClone("hist");
    plot::draw_label(eBeamEnergy, pBeamEnergy, detector);
    t1->SetFillColorAlpha(kWhite, 0);
    t1->SetTextFont(43);
    t1->SetTextSize(25);
    tptr1 = t1->AddText("Q^{2} Relative Resolution");
    tptr1->SetTextColor(plot::kMpOrange);
    t1->Draw();
    c.Print(fmt::format("{}Q2_rel_resolution.png", output_prefix).c_str());
  }
  eic::util::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}
