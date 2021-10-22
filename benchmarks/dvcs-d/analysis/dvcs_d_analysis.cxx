#include "pleaseIncludeMe.h"
int dvcs_d_analysis(const std::string& config_name)
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


  // Run this in multi-threaded mode if desired
  ROOT::EnableImplicitMT(kNumThreads);
  ROOT::RDataFrame d("events", rec_file);

  /*
  Block 1
  - default defines examples, filling Q2/x sim,rec res.
  */

  //event kinematics
  auto d0 = d.Define("Q2_sim", "InclusiveKinematicsTruth.Q2")
             .Define("Q2_rec", "InclusiveKinematicsElectron.Q2")
             .Define("Q2_res", combinatorial_diff_ratio, {"Q2_sim", "Q2_rec"})
             .Define("x_sim", "InclusiveKinematicsTruth.x")
             .Define("x_rec", "InclusiveKinematicsElectron.x")
             .Define("x_res", combinatorial_diff_ratio, {"x_sim", "x_rec"});

  //Q2
  auto h_Q2_sim = d0.Histo1D({"h_Q2_sim", "; GeV^2; counts", 100, -5, 25}, "Q2_sim");
  auto h_Q2_rec = d0.Histo1D({"h_Q2_rec", "; GeV^2; counts", 100, -5, 25}, "Q2_rec");
  auto h_Q2_res = d0.Histo1D({"h_Q2_res", ";      ; counts", 100, -1,  1}, "Q2_res");
  //x
  auto h_x_sim = d0.Histo1D({"h_x_sim", "; ; counts", 100, 0, +1}, "x_sim");
  auto h_x_rec = d0.Histo1D({"h_x_rec", "; ; counts", 100, 0, +1}, "x_rec");
  auto h_x_res = d0.Histo1D({"h_x_res", "; ; counts", 100, -1, 1}, "x_res");
  
  /*
  Block 2
  - Kong's examples on filtering events on the REC level
  - Reconstruct gamma, pt, eta 
  - scattered electron finder.
  */

  //y,Q2 cuts 
  auto kineCut = [](const ROOT::VecOps::RVec<float>& qsq, const ROOT::VecOps::RVec<float>& y_rec) { 
    if(qsq.size()<1||y_rec.size()<1) return 0;
    if(qsq[0] > 1. && qsq[0] < 10. && y_rec[0] < 0.95 && y_rec[0] > 0.01) return 1;
    else return 0;
  };

  auto d1 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("scatID_value","InclusiveKinematicsElectron.scatID.value")
             .Define("scatID_source","InclusiveKinematicsElectron.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("scatElec",findScatElec,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"}).Define("etaElec",getEta,{"scatElec"})
             .Define("gammaREC",findGamma,{"ReconstructedParticles"}).Define("MassREC",getMass,{"gammaREC"}).Define("gamma_rec_pt",getPt,{"gammaREC"}).Define("gamma_rec_eta",getEta,{"gammaREC"})
             .Define("protonREC",findScatProton,{"ReconstructedFFParticles"}).Define("proton_rec_eta",getEta,{"protonREC"}).Define("proton_rec_phi",getPhi,{"protonREC"})
             .Define("proton_rec_pt",getPt,{"protonREC"})
             .Define("t_REC",giveme_t_REC,{"protonREC","mcparticles"}).Define("t_REC_A",giveme_t,{"gammaREC","scatElec"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Q2_elec = d1.Histo1D({"h_Q2_elec", "; GeV^2; counts", 100, -5, 25}, "Q2_elec");
  auto h_y_elec = d1.Histo1D({"h_y_elec", "; ; counts", 100, 0, 1}, "y_elec");
  auto h_Eta_scatElec_REC = d1.Histo1D({"h_Eta_scatElec_REC",";eta; counts",100,-11,9}, "etaElec");
  auto h_Mass_gamma_REC = d1.Histo1D({"h_Mass_gamma_REC", "; GeV; counts", 1000, 0, 4}, "MassREC");
  auto h_Pt_gamma_REC = d1.Histo1D({"h_Pt_gamma_REC", "; GeV; counts", 50, 0, 2}, "gamma_rec_pt");
  auto h_Eta_gamma_REC = d1.Histo1D({"h_Eta_gamma_REC", "; ; counts", 100, -11, 9}, "gamma_rec_eta");
  auto h_Pt_proton_REC = d1.Histo1D({"h_Pt_proton_REC", "; GeV; counts", 50, 0, 5}, "proton_rec_pt");
  auto h_Eta_proton_REC = d1.Histo1D({"h_Eta_proton_REC", "; ; counts", 50, -11,9}, "proton_rec_eta");
  auto h_EtaPhi_proton_REC = d1.Histo2D({"h_EtaPhi_proton_REC",";eta;phi",100,-10,10,100,-PI,PI},"proton_rec_eta","proton_rec_phi");
  auto h_t_REC = d1.Histo1D({"h_t_REC", "; ; counts", 50, 0, 5}, "t_REC");
  auto h_t_REC_A = d1.Histo1D({"h_t_REC_A", "; ; counts", 50, 0, 5}, "t_REC_A");

  /*
  Block 3
  - Kong's examples on gen particles distributions, including 
  - e', VM, and VM daughters
  */

  auto d2 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("scatElecMC",findScatElecMC, {"mcparticles"}).Define("etaElecMC",getEta,{"scatElecMC"})
             .Define("gammaMC",findGammaMC,{"mcparticles"}).Define("MassMC",getMass,{"gammaMC"}).Define("gamma_mc_pt",getPt,{"gammaMC"}).Define("gamma_mc_eta",getEta,{"gammaMC"})
             .Define("protonMC",findScatProtonMC,{"mcparticles"}).Define("proton_mc_pt",getPt,{"protonMC"},Define("proton_mc_eta"),getEta,{"protonMC"})
             .Define("t_MC",giveme_t_MC,{"mcparticles"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Eta_scatElec_MC = d2.Histo1D({"h_Eta_scatElec_MC",";eta; counts",100,-11,9}, "etaElecMC");
  auto h_Mass_MC = d2.Histo1D({"h_Mass_MC",";Mass; counts",100,0,4}, "MassMC");
  auto h_Pt_gamma_MC = d2.Histo1D({"h_Pt_gamma_MC", "; GeV; counts", 50, 0, 2}, "gamma_mc_pt");
  auto h_Eta_gamma_MC = d2.Histo1D({"h_Eta_gamma_MC", "; ; counts", 100, -11, 9}, "gamma_mc_eta");
  auto h_Pt_proton_MC = d2.Histo1D({"h_Pt_proton_MC", "; GeV; counts", 50, 0, 5}, "proton_mc_pt");
  auto h_Eta_proton_MC = d2.Histo1D({"h_Eta_proton_MC", "; ; counts", 100, -11, 9}, "proton_mc_eta");
  auto h_t_MC = d2.Histo1D({"h_t_MC", "; ; counts", 50, 0, 2}, "t_MC");
 
  TString output_name_dir = output_prefix.c_str();
  TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");
  
  //Block 1
  
  h_Q2_sim->Write();
  h_Q2_rec->Write();
  h_Q2_res->Write();
  h_x_sim->Write();
  h_x_rec->Write();
  h_x_res->Write();

  //Block 2

  h_Q2_elec->Write();
  h_y_elec->Write();
  h_Eta_scatElec_REC->Write();
  h_Mass_gamma_REC->Write();
  h_Pt_gamma_REC->Write();
  h_Eta_gamma_REC->Write();
  h_Pt_proton_REC->Write();
  h_Eta_proton_REC->Write();
  h_EtaPhi_proton_REC->Write();
  h_t_REC->Write();
  h_t_REC_A->Write();

   //Block 3

  h_Eta_scatElec_MC->Write();
  h_Mass_MC->Write();
  h_Pt_gamma_MC->Write();
  h_Eta_gamma_MC->Write();
  h_Pt_proton_MC->Write();
  h_Eta_proton_MC->Write();
  h_t_MC->Write();

  output->Write();
  output->Close();

  // common_bench::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}
