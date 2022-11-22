#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "TCanvas.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"
#include "fmt/color.h"

R__LOAD_LIBRARY(libedm4eic.so)

#include "edm4eic/TrackParametersCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/InclusiveKinematicsCollection.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

auto p_track = [](std::vector<edm4eic::TrackParametersData> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(std::abs(1.0/(in[i].qOverP)));
  }
  return result;
};

auto momentum = [](std::vector<ROOT::Math::PxPyPzMVector> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(in[i].E());
  }
  return result;
};

auto theta = [](std::vector<ROOT::Math::PxPyPzMVector> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(in[i].Theta()*180/M_PI);
  }
  return result;
};

auto recfourvec = [](ROOT::VecOps::RVec<edm4eic::ReconstructedParticleData> const& in) {
  std::vector<ROOT::Math::PxPyPzMVector> result;
  ROOT::Math::PxPyPzMVector lv;
  for (size_t i = 0; i < in.size(); ++i) {
    lv.SetCoordinates(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].mass);
    result.push_back(lv);
  }
  return result;
};

auto delta_p = [](const std::vector<double>& tracks, const std::vector<double>& thrown) {
  std::vector<double> res;
  for (const auto& p1 : thrown) {
    for (const auto& p2 : tracks) {
      res.push_back(p1 - p2);
    }
  }
  return res;
};


void demo(const char* fname = "rec_dvcs.root"){

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green), "Running u_omega analysis...\n");

  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df("events", fname);

  using ROOT::Math::PxPyPzMVector;
  PxPyPzMVector p_ebeam = {0,0,-10, 0.000511};
  PxPyPzMVector p_pbeam = {0,0,275,  0.938 };

  auto q_vec = [=](PxPyPzMVector const& p) {
    return p_ebeam - p;
  };

  auto df0 = df.Define("n_parts", "ReconstructedParticles.size()")
               .Define("isQ2gt1", "InclusiveKinematicsTruth.Q2 > 1.0")
               .Define("n_Q2gt1", "isQ2gt1.size()");


  auto h_n_parts = df0.Histo1D({"h_n_parts", "; h_n_parts n", 10, 0, 10}, "n_parts");
  auto h_Q2      = df0.Histo1D({"h_Q2", "; Q^{2} [GeV^{2}/c^{2}]", 100, 0, 30}, "InclusiveKinematicsTruth.Q2");
  auto n_Q2gt1   = df0.Mean("n_Q2gt1");
  auto n_parts   = df0.Mean("n_parts");

  // ---------------------------
  // Do evaluation

  auto c = new TCanvas();
  h_Q2->DrawCopy();
  c->SaveAs("results/u_omega/Q2.png");
  c->SaveAs("results/u_omega/Q2.pdf");
  fmt::print("{} u_omega events Q2>1\n",*n_Q2gt1);
  fmt::print("{} tracks per event\n",*n_parts);

  c = new TCanvas();
  h_n_parts->DrawCopy();
  c->SaveAs("results/u_omega/n_parts.png");
  c->SaveAs("results/u_omega/n_parts.pdf");

}
