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
#include <utility>
#include <vector>

// Reconstuction functions
bool mom_sort_recon(eic::ReconstructedParticleData& part1, eic::ReconstructedParticleData& part2)
{
  return (part1.p.x * part1.p.x + part1.p.y * part1.p.y + part1.p.z * part1.p.z <
          part2.p.x * part2.p.x + part2.p.y * part2.p.y + part2.p.z * part2.p.z);
}

inline auto momenta_from_recon(const std::vector<eic::ReconstructedParticleData>& parts)
{
  std::vector<eic::ReconstructedParticleData> sort_parts = parts;
  sort(sort_parts.begin(), sort_parts.end(), mom_sort_recon);
  return sort_parts;
}

inline auto Q2_from_recon(const std::vector<eic::ReconstructedParticleData>& parts)
{
  std::vector<double> Q2Vec(parts.size());
  double              beamEnergy = 18;
  std::transform(parts.begin(), parts.end(), Q2Vec.begin(), [beamEnergy](const auto& part) {
    double q2 = pow(beamEnergy - part.energy, 2.0) - part.p.x * part.p.x - part.p.y * part.p.y -
                pow(beamEnergy - part.p.z, 2.0);
    return -q2;
  });
  return Q2Vec;
}

// Simulation functions
bool mom_sort_sim(dd4pod::Geant4ParticleData& part1, dd4pod::Geant4ParticleData& part2)
{
  return (part1.psx * part1.psx + part1.psy * part1.psy + part1.psz * part1.psz <
          part2.psx * part2.psx + part2.psy * part2.psy + part2.psz * part2.psz);
}

inline auto momenta_from_sim(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<dd4pod::Geant4ParticleData> sort_parts = parts;
  sort(sort_parts.begin(), sort_parts.end(), mom_sort_sim);
  return sort_parts;
}

inline auto Q2_from_sim(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<double> Q2Vec(parts.size());
  double              beamEnergy = 18;
  std::transform(parts.begin(), parts.end(), Q2Vec.begin(), [beamEnergy](const auto& part) {
    double energy = sqrt(part.psx * part.psx + part.psy * part.psy + part.psz * part.psz +
                         part.mass * part.mass);
    double q2     = pow(beamEnergy - energy, 2.0) - part.psx * part.psx - part.psy * part.psy -
                pow(part.psz - beamEnergy, 2.0);
    return -q2;
  });
  return Q2Vec;
}

inline auto elec_PID_sim(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<dd4pod::Geant4ParticleData> electrons;
  for (auto part : parts) {
    if (part.pdgID == 11) {
      electrons.push_back(part);
    }
  }
  return electrons;
}

inline auto mass_from_recon(const std::vector<eic::ReconstructedParticleData>& parts)
{
  std::vector<double> mass(parts.size());
  std::transform(parts.begin(), parts.end(), mass.begin(),
                 [](const auto& part) { return part.mass; });
  return mass;
}

inline auto mass_random(const std::vector<double>& massVec)
{
  std::vector<double> mass(massVec.size());
  TRandom3            rand;
  std::transform(massVec.begin(), massVec.end(), mass.begin(),
                 [&rand](const auto& inmass) { return inmass * rand.Gaus(1.0, 0.2); });
  return mass;
}

inline auto mass_from_simulation(const std::vector<ROOT::Math::PxPyPzMVector>& momenta)
{
  std::vector<double> mass(momenta.size());
  std::transform(momenta.begin(), momenta.end(), mass.begin(),
                 [](const auto& mom) { return mom.mass(); });
  return mass;
}

inline auto getSimPID(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<int> pid(parts.size());
  std::transform(parts.begin(), parts.end(), pid.begin(),
                 [](const auto& part) { return part.pdgID; });
  return pid;
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

  // Particle number enumeration
  enum sidis_particle_ID {
    electron = 11,
    pi_0     = 111,
    pi_plus  = 211,
    pi_minus = -211,
    k_0      = 311,
    k_plus   = 321,
    k_minus  = -321,
    proton   = 2212,
    neutron  = 2112
  };

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

  auto d0 = d.Define("sorted_recon", momenta_from_recon, {"DummyReconstructedParticles"})
                .Define("Q2_recon", Q2_from_recon, {"sorted_recon"})
                .Define("elec_Q2_recon", "Q2_recon[0]")
                .Define("sorted_sim", momenta_from_sim, {"mcparticles2"})
                .Define("Q2_sim", Q2_from_sim, {"sorted_sim"})
                .Define("elec_Q2_sim", "Q2_sim[0]")
                .Define("electrons_sim", elec_PID_sim, {"sorted_sim"})
                .Define("Q2_sim_elec_pid", Q2_from_sim, {"sorted_sim"})
                .Define("elec_Q2_sim_pid", "Q2_sim_elec_pid[0]");
  // Testing script
  /*
    auto dis = d0.Display({"Q2_sim", "elec_Q2_recon", "elec_Q2_sim", "elec_Q2_sim_pid"});
    dis -> Print();
    cout << *d0.Max<double>("elec_Q2_recon") << " " << *d0.Min<double>("elec_Q2_recon") << endl;
  */
  // cout << *d0.Count() << endl;
  // Momentum
  // auto h_mom_sim = d0.Histo1D({"h_mom_sim", "; GeV; counts", 100, 0, 50}, "mom_sim");
  // auto h_mom_rec = d0.Histo1D({"h_mom_rec", "; GeV; counts", 100, 0, 50}, "mom_rec");
  // Q2
  auto h_Q2_sim = d0.Histo1D({"h_Q2_sim", "; GeV; counts", 100, -50, 50}, "elec_Q2_sim");
  auto h_Q2_rec = d0.Histo1D({"h_Q2_rec", "; GeV; counts", 100, -50, 50}, "elec_Q2_recon");

  TH1D* h_Q2_res      = (TH1D*)h_Q2_sim->Clone();
  TH1D* h_Q2_rec_copy = (TH1D*)h_Q2_rec->Clone();
  h_Q2_res->Scale(1.0 / h_Q2_res->Integral());
  h_Q2_res->Add(h_Q2_rec_copy, -1.0 / h_Q2_rec_copy->Integral());

  TFitResultPtr f1 = h_Q2_res->Fit("gaus", "S");
  f1->Print("V");
  /*
    printf("chisq %f A %f mean %f sigma %f \n", f1 -> Chi2(),
                                                f1 -> GetParameter(0),
                                                f1 -> GetParameter(1),
                                                f1 -> GetParameter(2));
  */

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
    h2.Scale(2.0);
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
    plot::draw_label(18, 275, detector);
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
  eic::util::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}
