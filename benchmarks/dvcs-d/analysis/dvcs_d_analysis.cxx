#include "pleaseIncludeMe.h"
int dvcs_d_analysis(const std::string& config_name, const int vm_type=1)
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
  
  /*
  Block 2
  - Kong's examples on filtering events on the REC level
  - Reconstruct VM (e.g., phi) thru decays, inv_mass, pt, eta 
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
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles"})
             .Define("scatID_value","InclusiveKinematicsElectron.scatID.value")
             .Define("scatID_source","InclusiveKinematicsElectron.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("scatElec",findScatElec,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"}).Define("etaElec",getEta,{"scatElec"})
             .Define("vm", vector_sum, {"p1","p2"}).Define("Mass",getMass,{"vm"}).Define("vm_rec_pt", getPtVM, {"vm"}).Define("vm_rec_eta", getEtaVM, {"vm"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Q2_elec = d1.Histo1D({"h_Q2_elec", "; GeV^2; counts", 100, -5, 25}, "Q2_elec");
  auto h_y_elec = d1.Histo1D({"h_y_elec", "; ; counts", 100, 0, 1}, "y_elec");
  auto h_Eta_scatElec_REC = d1.Histo1D({"h_Eta_scatElec_REC",";eta; counts",100,-11,9}, "etaElec");
  auto h_Mass_REC = d1.Histo1D({"h_Mass_REC", "; GeV; counts", 1000, 0, 4}, "Mass");
  auto h_Pt_VM_REC = d1.Histo1D({"h_Pt_VM_REC", "; GeV; counts", 50, 0, 2}, "vm_rec_pt");
  auto h_Eta_VM_REC = d1.Histo1D({"h_Eta_VM_REC", "; ; counts", 100, -11, 9}, "vm_rec_eta");

  /*
  Block 3
  - Kong's examples on gen particles distributions, including 
  - e', VM, and VM daughters
  */

  auto d2 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("scatElecMC",findScatElecMC, {"mcparticles"}).Define("etaElecMC",getEta,{"scatElecMC"})
             .Define("VMMC",findVMMC,{"mcparticles"}).Define("MassMC",getMass,{"VMMC"}).Define("vm_mc_pt",getPtVM,{"VMMC"}).Define("vm_mc_eta",getEtaVM,{"VMMC"})
             .Define("VMMC_daugPlus",findVM_DaugPlus_MC,{"mcparticles"}).Define("etaVMMC_daugPlus",getEta,{"VMMC_daugPlus"}).Define("ptVMMC_daugPlus",getPt,{"VMMC_daugPlus"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Eta_scatElec_MC = d2.Histo1D({"h_Eta_scatElec_MC",";eta; counts",100,-11,9}, "etaElecMC");
  auto h_Mass_MC = d2.Histo1D({"h_Mass_MC",";Mass; counts",100,0,4}, "MassMC");
  auto h_Pt_VM_MC = d2.Histo1D({"h_Pt_VM_MC", "; GeV; counts", 50, 0, 2}, "vm_mc_pt");
  auto h_Eta_VM_MC = d2.Histo1D({"h_Eta_VM_MC", "; ; counts", 100, -11, 9}, "vm_mc_eta");
  auto h_Pt_VMdaugPlus_MC = d2.Histo1D({"h_Pt_VMdaugPlus_MC",";pt; counts",100,0,9}, "ptVMMC_daugPlus");
  auto h_Eta_VMdaugPlus_MC = d2.Histo1D({"h_Eta_VMdaugPlus_MC",";eta; counts",100,-11,9}, "etaVMMC_daugPlus");

  /*
  Block 4 
  - resolution on t
  - t resolution
  */

  auto d3 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles"})
             .Define("scatID_value","InclusiveKinematicsElectron.scatID.value")
             .Define("scatID_source","InclusiveKinematicsElectron.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("scatElec",findScatElec,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"})
             .Define("vm", vector_sum, {"p1","p2"})
             .Define("t_rec", giveme_t, {"vm","scatElec"})
             .Define("scatElecMC",findScatElecMC, {"mcparticles"})
             .Define("VMMC",findVMMC,{"mcparticles"})
             .Define("t_MC",giveme_t,{"VMMC","scatElecMC"})
             .Define("t_res",giveme_resolution,{"t_MC","t_rec"})     
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_t_rec = d3.Histo1D({"h_t_rec", "; GeV^{2}; counts",50,0,2}, "t_rec");
  auto h_t_MC = d3.Histo1D({"h_t_MC",";t; counts",50,0,2}, "t_MC");
  auto h_t_res = d3.Histo1D({"h_t_res",";res; counts",100,-1,1},"t_res");
  auto h_t_res_2D = d3.Histo2D({"h_t_res_2D",";-t;res",50,0,2,100,-1,1},"t_MC","t_res");

  /*
  Block 5
  - single particle (vm, dvcs photon) efficiency, fake, resolution
  */

  auto d4 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles"})
             .Define("vm", vector_sum, {"p1","p2"}).Define("vm_rec_pt",getPtVM,{"vm"})
             .Define("VMMC",findVMMC,{"mcparticles"}).Define("vm_mc_pt",getPtVM,{"VMMC"})
             .Define("vm_mc_match_rec",findVM_MC_match_REC,{"VMMC","vm"})
             .Define("vm_mcMrec_pt",getPtVM,{"vm_mc_match_rec"})
             .Define("vm_rec_not_match_mc",findVM_REC_NOT_match_MC,{"vm","VMMC"})
             .Define("vm_recNMmc_pt",getPtVM,{"vm_rec_not_match_mc"})
             .Define("vm_res_pt",resolution_MC_match_REC,{"VMMC","vm"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_Pt_VM_MC_total = d4.Histo1D({"h_Pt_VM_MC_total", "; GeV; counts", 50, 0, 2}, "vm_mc_pt");
  auto h_Pt_VM_MC_match = d4.Histo1D({"h_Pt_VM_MC_match", "; GeV; counts", 50, 0, 2}, "vm_mcMrec_pt");
  auto h_Pt_VM_REC_total = d4.Histo1D({"h_Pt_VM_REC_total", "; GeV; counts", 50, 0, 2}, "vm_rec_pt");
  auto h_Pt_VM_REC_not_match = d4.Histo1D({"h_Pt_VM_REC_not_match", "; GeV; counts", 50, 0, 2}, "vm_recNMmc_pt");
  auto h_Pt_VM_MC_res = d4.Histo2D({"h_Pt_VM_MC_res",";pt;res",50,0,2,100,-1,1},"vm_mc_pt","vm_res_pt");

  /*
  Block 6
  - t rec with veto conditions applied, including FF detectors
  */

  auto eventVetoCut = [](const std::vector<eic::ReconstructedParticleData>& FF,
    const std::vector<eic::ReconstructedParticleData>& parts)
  { 
    bool keepThisEvent_ = true;
    if(FF.size()>0) keepThisEvent_ = false;
    int mult=0;
    for(auto&i1 : parts){
      TLorentzVector p(i1.p.x,i1.p.y,i1.p.z,i1.mass);
      if(p.Pt()>0.05&&fabs(p.Eta())<3.5) mult++;
    }
    if(mult>3) keepThisEvent_ = false;

    return keepThisEvent_;
  };

  auto d5 = d.Define("Q2_elec", "InclusiveKinematicsElectron.Q2")
             .Define("y_elec", "InclusiveKinematicsElectron.y")
             .Define("p1", momenta_from_reconstruction_plus, {"ReconstructedChargedParticles"})
             .Define("p2", momenta_from_reconstruction_minus, {"ReconstructedChargedParticles"})
             .Define("scatID_value","InclusiveKinematicsElectron.scatID.value")
             .Define("scatID_source","InclusiveKinematicsElectron.scatID.source")
             .Define("scatID_cand_value",scatID_cand_value, {"scatID_value"})
             .Define("scatID_cand_source",scatID_cand_value, {"scatID_source"})
             .Define("scatElec",findScatElec,{"ReconstructedChargedParticles","scatID_cand_value","scatID_cand_source"})
             .Define("vm", vector_sum, {"p1","p2"})
             .Define("t_rec", giveme_t, {"vm","scatElec"})
             .Filter(eventVetoCut,{"ReconstructedFFParticles","ReconstructedChargedParticles"})
             .Filter(kineCut,{"Q2_elec","y_elec"});

  auto h_t_rec_veto = d5.Histo1D({"h_t_rec_veto", "; GeV^{2}; counts",50,0,2}, "t_rec");

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
  h_Mass_REC->Write();
  h_Pt_VM_REC->Write();
  h_Eta_VM_REC->Write();

  //Block 3
  h_Eta_scatElec_MC->Write();
  h_Mass_MC->Write();
  h_Pt_VM_MC->Write();
  h_Eta_VM_MC->Write();
  h_Pt_VMdaugPlus_MC->Write();
  h_Eta_VMdaugPlus_MC->Write();

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

  //Block 6
  h_t_rec_veto->Write();

  output->Write();
  output->Close();

  // common_bench::write_test({dis_Q2_resolution}, fmt::format("{}dis_electrons.json", output_prefix));

  return 0;
}