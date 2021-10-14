#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include "TFile.h"

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"
#include "eicd/InclusiveKinematicsData.h"
#include "eicd/ReconstructedParticleData.h"

//particles properties
  auto momenta_from_reconstruction(const std::vector<eic::ReconstructedParticleData>& parts) {
    std::vector<ROOT::Math::PxPyPzEVector> momenta{parts.size()};
    std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
      return ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy};
    });
    return momenta;
  }
  auto getPt(const std::vector<ROOT::Math::PxPyPzEVector>& mom) {
    std::vector<double> PtVec(mom.size() );
    ROOT::Math::PxPyPzEVector beamMom = {0, 0, -18, 18};
    std::transform(mom.begin(), mom.end(), PtVec.begin(), [beamMom](const auto& part) {
      return part.Pt();
    });
    return PtVec;
  }

int diffractive_phi_analysis(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const std::string test_tag      = config["test_tag"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running diffractive phi analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix: {}\n", output_prefix);
  fmt::print(" - test tag: {}\n", test_tag);

  // disable test for now - Kong.Tu
  // create our test definition
  // test_tag
  // common_bench::Test dis_Q2_resolution{
  //     {{"name", fmt::format("{}_Q2_resolution", test_tag)},
  //      {"title", "DIS Q2 resolution"},
  //      {"description",
  //       fmt::format("DIS Q2 resolution with {}, estimated using a Gaussian fit.", detector)},
  //      {"quantity", "resolution (in %)"},
  //      {"target", "0.1"}}};

  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT(kNumThreads);
  ROOT::RDataFrame d("events", rec_file);

  auto combinatorial_diff_ratio = [] (
      const ROOT::VecOps::RVec<float>& v1,
      const ROOT::VecOps::RVec<float>& v2
  ) {
    std::vector<float> v;
    for (auto& i1: v1) {
      for (auto& i2: v2) {
        if (i1 != 0) {
          v.push_back((i1-i2)/i1);
        }
      }
    }
    return v;
  };

  //event kinematics
  auto d0 = d.Define("Q2_sim", "InclusiveKinematicsTruth.Q2")
             .Define("Q2_rec", "InclusiveKinematicsElectron.Q2")
             .Define("Q2_res", combinatorial_diff_ratio, {"Q2_sim", "Q2_rec"})
             .Define("x_sim", "InclusiveKinematicsTruth.x")
             .Define("x_rec", "InclusiveKinematicsElectron.x")
             .Define("x_res", combinatorial_diff_ratio, {"x_sim", "x_rec"})
             ;

  //Q2
  auto h_Q2_sim = d0.Histo1D({"h_Q2_sim", "; GeV^2; counts", 100, -5, 25}, "Q2_sim");
  auto h_Q2_rec = d0.Histo1D({"h_Q2_rec", "; GeV^2; counts", 100, -5, 25}, "Q2_rec");
  auto h_Q2_res = d0.Histo1D({"h_Q2_res", ";      ; counts", 100, -1,  1}, "Q2_res");
  //x
  auto h_x_sim = d0.Histo1D({"h_x_sim", "; ; counts", 100, 0, +1}, "x_sim");
  auto h_x_rec = d0.Histo1D({"h_x_rec", "; ; counts", 100, 0, +1}, "x_rec");
  auto h_x_res = d0.Histo1D({"h_x_res", "; ; counts", 100, -1, 1}, "x_res");

  
  auto d1 = d.Define("p", momenta_from_reconstruction, {"ReconstructedChargedParticles"}).Define("Pt", getPt, {"p"});
  auto h_Pt_rec = d1.Histo1D({"h_Pt_rec", "; GeV; counts", 100, 0, 25}, "Pt");


  TString output_name_dir = output_prefix.c_str();
  TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");
  h_Q2_sim->Write();
  h_Q2_rec->Write();
  h_Q2_res->Write();

  h_x_sim->Write();
  h_x_rec->Write();
  h_x_res->Write();

  h_Pt_rec->Write();

  output->Write();
  output->Close();

  // common_bench::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}
