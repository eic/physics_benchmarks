#include "pleaseIncludeMe.h"
int test_phi_analysis(const std::string& config_name, const int vm_type=1, const int mc_type=0)
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

  //CHOOSE which VM, 0 = rho, 1 = phi, 2 = j/psi
  which_vm = vm_type;
  which_mc = mc_type;

  //y,Q2 cuts 
  auto kineCut = [](const ROOT::VecOps::RVec<float>& qsq, const ROOT::VecOps::RVec<float>& y_rec) { 
    if(qsq.size()<1||y_rec.size()<1) return 0;
    if(qsq[0] > 1. && qsq[0] < 10. && y_rec[0] < 0.95 && y_rec[0] > 0.01) return 1;
    else return 0;
  };
  /*
  Block 4 
  - resolution on t
  - t resolution
  */

  auto d3 = d.Define("Q2_elec", "InclusiveKinematicsTruth.Q2")
             .Define("y_elec", "InclusiveKinematicsTruth.y")
             .Define("scatID_value","InclusiveKinematicsTruth.scatID.value")
             .Define("scatID_source","InclusiveKinematicsTruth.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"})
             .Define("scatElec",findScatElecTest,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"})
             .Define("vm", vector_sum, {"p1","p2"})
             .Define("t_rec", giveme_t_A, {"vm","scatElec","mcparticles"})
             .Define("scatElecMC",findScatElecMC, {"mcparticles"})
             .Define("VMMC",findVMMC,{"mcparticles"})
             .Define("t_MC",giveme_t_MC_E,{"VMMC","scatElecMC","mcparticles"})
             .Define("t_res",giveme_resolution,{"t_MC","t_rec"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_t_rec = d3.Histo1D({"h_t_rec", "; GeV^{2}; counts",100,0,0.2}, "t_rec");
  auto h_t_MC = d3.Histo1D({"h_t_MC",";t; counts",100,0,0.2}, "t_MC");
  auto h_t_res = d3.Histo1D({"h_t_res",";res; counts",100,-1,1},"t_res");
  auto h_t_res_2D = d3.Histo2D({"h_t_res_2D",";-t;res",100,0,0.2,100,-1,1},"t_MC","t_res");

  /*
  Block 5
  - single particle (vm, dvcs photon) efficiency, fake, resolution
  */

  auto d4 = d.Define("Q2_elec", "InclusiveKinematicsTruth.Q2")
             .Define("y_elec", "InclusiveKinematicsTruth.y")
             .Define("scatID_value","InclusiveKinematicsTruth.scatID.value")
             .Define("scatID_source","InclusiveKinematicsTruth.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("scatElec",findScatElecTest,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"}).Define("e_rec_pt",getPt,{"scatElec"})
             .Define("scatElecMC",findScatElecMC, {"mcparticles"}).Define("e_mc_pt",getPt,{"scatElecMC"})
             .Define("e_res_pt",resolution_MC_match_REC_electron,{"scatElecMC","scatElec"})
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"})
             .Define("vm", vector_sum, {"p1","p2"}).Define("vm_rec_pt",getPtVM,{"vm"})
             .Define("VMMC",findVMMC,{"mcparticles"}).Define("vm_mc_pt",getPtVM,{"VMMC"})
             .Define("vm_mc_match_rec",findVM_MC_match_REC,{"VMMC","vm"})
             .Define("vm_mcMrec_pt",getPtVM,{"vm_mc_match_rec"})
             .Define("vm_rec_not_match_mc",findVM_REC_NOT_match_MC,{"vm","VMMC"})
             .Define("vm_recNMmc_pt",getPtVM,{"vm_rec_not_match_mc"})
             .Define("vm_res_pt",resolution_MC_match_REC,{"VMMC","vm"})
             .Define("protonMC",findScatProtonMC,{"mcparticles"}).Define("proton_mc_p",getP,{"protonMC"})
             .Define("VMMC_daugPlus",findVM_DaugPlus_MC,{"mcparticles"}).Define("ptVMMC_daugPlus",getPt,{"VMMC_daugPlus"})
             .Define("VMMC_daugPlusMatchREC",findVM_MC_match_REC,{"VMMC_daugPlus","p1"}).Define("ptVMMC_daugPlus_match",getPt,{"VMMC_daugPlusMatchREC"})
             .Define("p1_rec_pt",getPt,{"p1"}).Define("p1_res_pt",resolution_MC_match_REC,{"VMMC_daugPlus","p1"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Pt_VM_MC_total = d4.Histo1D({"h_Pt_VM_MC_total", "; GeV; counts", 50, 0, 2}, "vm_mc_pt");
  auto h_Pt_VM_MC_match = d4.Histo1D({"h_Pt_VM_MC_match", "; GeV; counts", 50, 0, 2}, "vm_mcMrec_pt");
  auto h_Pt_VM_REC_total = d4.Histo1D({"h_Pt_VM_REC_total", "; GeV; counts", 50, 0, 2}, "vm_rec_pt");
  auto h_Pt_VM_REC_not_match = d4.Histo1D({"h_Pt_VM_REC_not_match", "; GeV; counts", 50, 0, 2}, "vm_recNMmc_pt");
  auto h_Pt_VM_MC_res = d4.Histo2D({"h_Pt_VM_MC_res",";pt;res",50,0,2,300,-0.15,0.15},"vm_mc_pt","vm_res_pt");
  auto h_P_proton_MC = d4.Histo1D({"h_P_proton_MC", "; GeV; counts", 200, 0, 200}, "proton_mc_p");
  auto h_Pt_track_REC = d4.Histo1D({"h_Pt_track_REC", "; GeV; counts", 50, 0, 2}, "p1_rec_pt");
  auto h_Pt_track_MC = d4.Histo1D({"h_Pt_track_MC", "; GeV; counts", 50, 0, 2}, "ptVMMC_daugPlus");
  auto h_Pt_track_MC_match = d4.Histo1D({"h_Pt_track_MC_match", "; GeV; counts", 50, 0, 2}, "ptVMMC_daugPlus_match");
  auto h_Pt_track_MC_res = d4.Histo2D({"h_Pt_track_MC_res",";pt;res",50,0,2,300,-0.15,0.15},"ptVMMC_daugPlus","p1_res_pt");
  auto h_Pt_e_REC = d4.Histo1D({"h_Pt_e_REC", "; GeV; counts", 50, 0, 2}, "e_rec_pt");
  auto h_Pt_e_MC = d4.Histo1D({"h_Pt_e_MC", "; GeV; counts", 50, 0, 2}, "e_mc_pt");
  auto h_Pt_e_MC_res = d4.Histo2D({"h_Pt_e_MC_res",";pt;res",50,0,2,300,-0.15,0.15},"e_mc_pt","e_res_pt");

 
  TString output_name_dir = output_prefix.c_str();
  TFile* output = new TFile(output_name_dir+"_output.root","RECREATE");
  
  //Block 4
  h_t_rec->Write();
  h_t_MC->Write();
  h_t_res->Write();
  h_t_res_2D->Write();

  //Block 5
  h_Pt_VM_MC_total->Write();
  h_Pt_VM_MC_match->Write();
  h_Pt_VM_REC_total->Write();
  h_Pt_VM_REC_not_match->Write();
  h_Pt_VM_MC_res->Write();
  h_P_proton_MC->Write();
  h_Pt_track_REC->Write();
  h_Pt_track_MC->Write();
  h_Pt_track_MC_match->Write();
  h_Pt_track_MC_res->Write();
  h_Pt_e_REC->Write();
  h_Pt_e_MC->Write();
  h_Pt_e_MC_res->Write();

  output->Write();
  output->Close();

  return 0;
}
