#include "pleaseIncludeMe.h"
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
  

  auto d1 = d.Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles"})
             .Define("vm", vector_sum, {"p1","p2"}).Define("Pt2",getPt2,{"vm"}).Define("Mass",getMass,{"vm"});

  auto h_Pt_rec = d1.Histo1D({"h_Pt_rec", "; GeV; counts", 100, 0, 25}, "Pt2");
  auto h_Mass_rec = d1.Histo1D({"h_Mass_rec", "; GeV; counts", 100, 0, 4}, "Mass");

  TString output_name_dir = output_prefix.c_str();
  TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");
  h_Q2_sim->Write();
  h_Q2_rec->Write();
  h_Q2_res->Write();

  h_x_sim->Write();
  h_x_rec->Write();
  h_x_res->Write();

  h_Pt_rec->Write();
  h_Mass_rec->Write();

  output->Write();
  output->Close();

  TFile* output = new TFile(output_prefix.c_str()+"_output.root","RECREATE");
  h_Q2_sim.Write();
  
  // common_bench::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}
